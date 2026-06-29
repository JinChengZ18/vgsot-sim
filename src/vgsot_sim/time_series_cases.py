from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional

import numpy as np
from tqdm import tqdm

from .configs import (
    PhysicalConstantsConfig,
    SotOnlyConstantCurrentConfig,
    SotSwitchingNoVcmaConfig,
    TerminalVoltageControlConfig,
)
from .dynamic_switching import switching
from .dynamic_switching_vector import switching_vector
from .electronic import electronic
from .initialize import init
from .tmr import tmr
from .thermal import self_heating_step
from .material_temperature import ms_of_T, ki_of_T


@dataclass
class SimResult:
    time_s: np.ndarray
    mz: np.ndarray
    r_mtj: np.ndarray
    v_mtj: np.ndarray
    i_sot: np.ndarray
    switch_energy_j: float
    theta: Optional[np.ndarray] = None
    phi: Optional[np.ndarray] = None
    v1: Optional[np.ndarray] = None
    v2: Optional[np.ndarray] = None
    v3: Optional[np.ndarray] = None
    # Self-heating diagnostics. When `enable_self_heating=True` in the
    # stepper, these hold the per-step instantaneous temperature and the
    # corresponding T-corrected material parameters fed into LLG.
    T_K: Optional[np.ndarray] = None
    Ms_T: Optional[np.ndarray] = None
    Ki_T: Optional[np.ndarray] = None


@dataclass
class SweepResult:
    time_s: np.ndarray
    mz_curves: Dict[str, np.ndarray]
    r_mtj_curves: Dict[str, np.ndarray]
    pulse_curves: Dict[str, np.ndarray]
    switch_energy_j: Dict[str, float]
    pulse_ylabel: str

    @property
    def curves(self) -> Dict[str, np.ndarray]:
        return self.mz_curves


def _maybe_tqdm(iterable: Iterable[int], show_progress: bool, desc: str = ""):
    return tqdm(iterable, desc=desc) if show_progress else iterable


def _compute_switch_energy_j(
    time_s: np.ndarray,
    v_mtj: np.ndarray,
    r_mtj: np.ndarray,
    i_sot: np.ndarray,
    constants: PhysicalConstantsConfig,
) -> float:
    safe_r = np.where(np.abs(r_mtj) > 1e-30, r_mtj, np.nan)
    p_mtj = np.where(np.isfinite(safe_r), (v_mtj ** 2) / safe_r, 0.0)
    p_sot = (i_sot ** 2) * constants.R_W
    p_total = p_mtj + p_sot
    active = (np.abs(v_mtj) > 1e-15) | (np.abs(i_sot) > 1e-18)
    if not np.any(active):
        return 0.0
    # `trapezoid` is the post-NumPy-2.0 name; fall back to the deprecated
    # `trapz` for older interpreters.
    _trap = getattr(np, "trapezoid", np.trapz)
    return float(_trap(p_total[active], time_s[active]))


def run_piecewise_terminal_voltage(
    cfg: TerminalVoltageControlConfig,
    *,
    show_progress: bool = True,
    enable_self_heating: bool = False,
    T_ambient_K: float = 300.0,
    demag_mode: str = "ellipsoid",
    integrator: str = "cayley",
    sigma_SH=None,
    rng=None,
) -> SimResult:
    constants = cfg.constants
    sim_mid2_step = cfg.sim_mid2_step if cfg.sim_mid2_step is not None else cfg.sim_mid1_step

    n_steps = cfg.sim_end_step + 1
    time_idx = np.arange(cfg.sim_start_step - 1, cfg.sim_end_step + 1)
    time_s = time_idx * constants.t_step

    mz_arr = np.zeros(n_steps)
    theta_arr = np.zeros(n_steps)
    phi_arr = np.zeros(n_steps)
    r_arr = np.zeros(n_steps)
    v_mtj_arr = np.zeros(n_steps)
    i_sot_arr = np.zeros(n_steps)
    v1_arr = np.zeros(n_steps)
    v2_arr = np.zeros(n_steps)
    v3_arr = np.zeros(n_steps)
    if enable_self_heating:
        T_arr  = np.full(n_steps, T_ambient_K, dtype=float)
        Ms_arr = np.full(n_steps, constants.Ms, dtype=float)
        Ki_arr = np.full(n_steps, constants.Ki, dtype=float)
    else:
        T_arr = Ms_arr = Ki_arr = None

    r0, theta, mz, phi = init(cfg.pap, constants, rng=rng)
    r_arr[0], theta_arr[0], mz_arr[0], phi_arr[0] = r0, theta, mz, phi
    T_now = T_ambient_K

    def stage_params(i: int):
        if i < cfg.sim_mid1_step:
            return cfg.v_stage1, cfg.estt_stage1, cfg.esot_stage1
        if i < sim_mid2_step:
            return cfg.v_stage2, cfg.estt_stage2, cfg.esot_stage2
        return cfg.v_stage3, cfg.estt_stage3, cfg.esot_stage3

    loop = _maybe_tqdm(range(cfg.sim_start_step - 1, cfg.sim_end_step), show_progress, desc="terminal_voltage")
    for i in loop:
        (v1, v2, v3), estt, esot = stage_params(i)
        r_mtj = r_arr[i]
        i_sot, v_mtj = electronic(v1, v2, v3, r_mtj, constants)

        if enable_self_heating:
            V_SOT_effective = i_sot * constants.R_W
            T_now = self_heating_step(
                T_now, v_mtj, V_SOT_effective, r_mtj, constants,
                dt=constants.t_step, T_0=T_ambient_K,
            )
            Ms_T = ms_of_T(T_now, constants)
            Ki_T = ki_of_T(T_now, constants)
            T_arr[i + 1] = T_now
            Ms_arr[i + 1] = Ms_T
            Ki_arr[i + 1] = Ki_T
        else:
            Ms_T = None
            Ki_T = None

        if integrator == "cayley":
            # Cartesian + Cayley-rotation step (norm-preserving, accepts an
            # explicit sigma_SH 3-vector). Mirrors run_piecewise_direct_excitation.
            m_vec = np.array([
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                np.cos(theta),
            ])
            sig_SH = np.array([-1.0, 0.0, 0.0]) if sigma_SH is None else np.asarray(sigma_SH, dtype=float)
            m_vec_new = switching_vector(
                m_vec, v_mtj, i_sot, r_mtj, estt, esot,
                VNV=cfg.vnv, NON=cfg.non, R_SOT_FL_DL=cfg.r_sot_fl_dl,
                sigma_SH=sig_SH,
                constants=constants,
                Ki_T=Ki_T, Ms_T=Ms_T, T=(T_now if enable_self_heating else None), demag_mode=demag_mode,
                rng=rng,
            )
            mz = float(m_vec_new[2])
            theta_tmp = float(np.arccos(np.clip(m_vec_new[2], -1.0, 1.0)))
            phi_tmp = float(np.arctan2(m_vec_new[1], m_vec_new[0]))
        else:
            mz, phi_tmp, theta_tmp = switching(
                v_mtj,
                i_sot,
                r_mtj,
                theta,
                phi,
                estt,
                esot,
                VNV=cfg.vnv,
                NON=cfg.non,
                R_SOT_FL_DL=cfg.r_sot_fl_dl,
                constants=constants,
                Ki_T=Ki_T, Ms_T=Ms_T, T=(T_now if enable_self_heating else None), demag_mode=demag_mode,
                rng=rng,
            )
        phi, theta = phi_tmp, theta_tmp
        r_next = tmr(v_mtj, mz, constants)

        v1_arr[i] = v1
        v2_arr[i] = v2
        v3_arr[i] = v3
        i_sot_arr[i] = i_sot
        v_mtj_arr[i] = v_mtj
        mz_arr[i + 1] = mz
        theta_arr[i + 1] = theta
        phi_arr[i + 1] = phi
        r_arr[i + 1] = r_next

    if n_steps > 1:
        v1_arr[-1] = v1_arr[-2]
        v2_arr[-1] = v2_arr[-2]
        v3_arr[-1] = v3_arr[-2]
        i_sot_arr[-1] = i_sot_arr[-2]
        v_mtj_arr[-1] = v_mtj_arr[-2]

    return SimResult(
        time_s=time_s,
        mz=mz_arr,
        r_mtj=r_arr,
        v_mtj=v_mtj_arr,
        i_sot=i_sot_arr,
        switch_energy_j=_compute_switch_energy_j(time_s, v_mtj_arr, r_arr, i_sot_arr, constants),
        theta=theta_arr,
        phi=phi_arr,
        v1=v1_arr,
        v2=v2_arr,
        v3=v3_arr,
        T_K=T_arr, Ms_T=Ms_arr, Ki_T=Ki_arr,
    )


def run_piecewise_direct_excitation(
    *,
    sim_start_step: int,
    sim_mid1_step: int,
    sim_mid2_step: Optional[int],
    sim_end_step: int,
    pap: int,
    v_mtj_stage1: float,
    v_mtj_stage2: float,
    v_mtj_stage3: float,
    i_sot_stage1: float,
    i_sot_stage2: float,
    i_sot_stage3: float,
    estt_stage1: int,
    esot_stage1: int,
    estt_stage2: int,
    esot_stage2: int,
    estt_stage3: int,
    esot_stage3: int,
    vnv: int,
    non: int,
    r_sot_fl_dl: float,
    show_progress: bool = True,
    constants: PhysicalConstantsConfig | None = None,
    enable_self_heating: bool = False,
    T_ambient_K: float = 300.0,
    demag_mode: str = "ellipsoid",
    integrator: str = "cayley",
    sigma_SH=None,
    rng=None,
) -> SimResult:
    constants = constants or PhysicalConstantsConfig()
    if sim_mid2_step is None:
        sim_mid2_step = sim_end_step

    n_steps = sim_end_step + 1
    time_idx = np.arange(sim_start_step - 1, sim_end_step + 1)
    time_s = time_idx * constants.t_step

    mz_arr = np.zeros(n_steps)
    theta_arr = np.zeros(n_steps)
    phi_arr = np.zeros(n_steps)
    r_arr = np.zeros(n_steps)
    v_arr = np.zeros(n_steps)
    i_sot_arr = np.zeros(n_steps)
    if enable_self_heating:
        T_arr  = np.full(n_steps, T_ambient_K, dtype=float)
        Ms_arr = np.full(n_steps, constants.Ms, dtype=float)
        Ki_arr = np.full(n_steps, constants.Ki, dtype=float)
    else:
        T_arr = Ms_arr = Ki_arr = None

    r0, theta, mz, phi = init(pap, constants, rng=rng)
    r_arr[0], theta_arr[0], mz_arr[0], phi_arr[0] = r0, theta, mz, phi
    T_now = T_ambient_K

    def stage_inputs(i: int):
        if i < sim_mid1_step:
            return v_mtj_stage1, i_sot_stage1, estt_stage1, esot_stage1
        if i < sim_mid2_step:
            return v_mtj_stage2, i_sot_stage2, estt_stage2, esot_stage2
        return v_mtj_stage3, i_sot_stage3, estt_stage3, esot_stage3

    loop = _maybe_tqdm(range(sim_start_step - 1, sim_end_step), show_progress, desc="direct_excitation")
    for i in loop:
        v_mtj, i_sot, estt, esot = stage_inputs(i)
        r_mtj = r_arr[i]

        if enable_self_heating:
            # Thermal RC update: V_SOT across the channel is the I_SOT
            # times its (configured) resistance. The MgO heat-sink path is
            # already inside `self_heating_step`.
            V_SOT_effective = i_sot * constants.R_W
            T_now = self_heating_step(
                T_now, v_mtj, V_SOT_effective, r_mtj, constants,
                dt=constants.t_step, T_0=T_ambient_K,
            )
            Ms_T = ms_of_T(T_now, constants)
            Ki_T = ki_of_T(T_now, constants)
            T_arr[i + 1] = T_now
            Ms_arr[i + 1] = Ms_T
            Ki_arr[i + 1] = Ki_T
        else:
            Ms_T = None
            Ki_T = None

        if integrator == "cayley":
            # Cartesian + Cayley-rotation step (norm-preserving, accepts
            # an explicit sigma_SH 3-vector). Convert (θ, φ) ↔ m at the
            # interface so the existing scalar arrays stay populated.
            m_vec = np.array([
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                np.cos(theta),
            ])
            sig_SH = np.array([-1.0, 0.0, 0.0]) if sigma_SH is None else np.asarray(sigma_SH, dtype=float)
            m_vec_new = switching_vector(
                m_vec, v_mtj, i_sot, r_mtj, estt, esot,
                VNV=vnv, NON=non, R_SOT_FL_DL=r_sot_fl_dl,
                sigma_SH=sig_SH,
                constants=constants,
                Ki_T=Ki_T, Ms_T=Ms_T, T=(T_now if enable_self_heating else None), demag_mode=demag_mode,
                rng=rng,
            )
            mz = float(m_vec_new[2])
            theta_tmp = float(np.arccos(np.clip(m_vec_new[2], -1.0, 1.0)))
            phi_tmp = float(np.arctan2(m_vec_new[1], m_vec_new[0]))
        else:
            mz, phi_tmp, theta_tmp = switching(
                v_mtj,
                i_sot,
                r_mtj,
                theta,
                phi,
                estt,
                esot,
                VNV=vnv,
                NON=non,
                R_SOT_FL_DL=r_sot_fl_dl,
                constants=constants,
                Ki_T=Ki_T, Ms_T=Ms_T, T=(T_now if enable_self_heating else None), demag_mode=demag_mode,
                rng=rng,
            )
        phi, theta = phi_tmp, theta_tmp
        r_next = tmr(v_mtj, mz, constants)

        v_arr[i] = v_mtj
        i_sot_arr[i] = i_sot
        mz_arr[i + 1] = mz
        theta_arr[i + 1] = theta
        phi_arr[i + 1] = phi
        r_arr[i + 1] = r_next

    if n_steps > 1:
        v_arr[-1] = v_arr[-2]
        i_sot_arr[-1] = i_sot_arr[-2]

    return SimResult(
        time_s=time_s,
        mz=mz_arr,
        r_mtj=r_arr,
        v_mtj=v_arr,
        i_sot=i_sot_arr,
        switch_energy_j=_compute_switch_energy_j(time_s, v_arr, r_arr, i_sot_arr, constants),
        theta=theta_arr,
        phi=phi_arr,
        T_K=T_arr, Ms_T=Ms_arr, Ki_T=Ki_arr,
    )


def terminal_voltage_control(cfg: TerminalVoltageControlConfig | None = None, *, show_progress: bool = True) -> SimResult:
    cfg = cfg or TerminalVoltageControlConfig()
    return run_piecewise_terminal_voltage(cfg, show_progress=show_progress)


def sot_only_constant_current(cfg: SotOnlyConstantCurrentConfig | None = None, *, show_progress: bool = True) -> SimResult:
    cfg = cfg or SotOnlyConstantCurrentConfig()
    return run_piecewise_direct_excitation(
        sim_start_step=cfg.sim_start_step,
        sim_mid1_step=cfg.sim_mid1_step,
        sim_mid2_step=cfg.sim_end_step,
        sim_end_step=cfg.sim_end_step,
        pap=cfg.pap,
        v_mtj_stage1=cfg.v_mtj_stage1,
        v_mtj_stage2=cfg.v_mtj_stage2,
        v_mtj_stage3=cfg.v_mtj_stage2,
        i_sot_stage1=cfg.i_sot_stage1,
        i_sot_stage2=cfg.i_sot_stage2,
        i_sot_stage3=cfg.i_sot_stage2,
        estt_stage1=0,
        esot_stage1=1,
        estt_stage2=0,
        esot_stage2=1,
        estt_stage3=0,
        esot_stage3=1,
        vnv=cfg.vnv,
        non=cfg.non,
        r_sot_fl_dl=cfg.r_sot_fl_dl,
        show_progress=show_progress,
        constants=cfg.constants,
    )


def sot_switching_no_vcma(cfg: SotSwitchingNoVcmaConfig | None = None, *, show_progress: bool = True) -> SweepResult:
    cfg = cfg or SotSwitchingNoVcmaConfig()
    time_s = None
    mz_curves, r_mtj_curves, pulse_curves, switch_energy_j = {}, {}, {}, {}
    for i_sot in cfg.i_sot_list:
        res = run_piecewise_direct_excitation(
            sim_start_step=cfg.sim_start_step,
            sim_mid1_step=cfg.sim_mid1_step,
            sim_mid2_step=cfg.sim_end_step,
            sim_end_step=cfg.sim_end_step,
            pap=cfg.pap,
            v_mtj_stage1=cfg.v_mtj,
            v_mtj_stage2=cfg.v_mtj,
            v_mtj_stage3=cfg.v_mtj,
            i_sot_stage1=i_sot,
            i_sot_stage2=cfg.i_sot_relax,
            i_sot_stage3=cfg.i_sot_relax,
            estt_stage1=0,
            esot_stage1=1,
            estt_stage2=0,
            esot_stage2=1,
            estt_stage3=0,
            esot_stage3=1,
            vnv=cfg.vnv,
            non=cfg.non,
            r_sot_fl_dl=cfg.r_sot_fl_dl,
            show_progress=show_progress,
            constants=cfg.constants,
        )
        label = f"I_SOT={i_sot*1e6:.1f}uA"
        if time_s is None:
            time_s = res.time_s
        mz_curves[label] = res.mz
        r_mtj_curves[label] = res.r_mtj
        pulse_curves[label] = res.i_sot * 1e6
        switch_energy_j[label] = res.switch_energy_j
    return SweepResult(
        time_s=time_s,
        mz_curves=mz_curves,
        r_mtj_curves=r_mtj_curves,
        pulse_curves=pulse_curves,
        switch_energy_j=switch_energy_j,
        pulse_ylabel=r"$I_{\mathrm{SOT}}$ ($\mu$A)",
    )


TIME_SERIES_CASES = (
    "terminal_voltage_control",
    "sot_only_constant_current",
    "sot_switching_no_vcma",
)
