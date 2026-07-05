from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from tqdm import trange

from .configs import PhysicalConstantsConfig, SerSotNoVcmaThermalConfig
from .time_series_cases import run_piecewise_direct_excitation


@dataclass
class SerResult:
    """Switching-error-rate result.

    The legacy field `ser` stays for backward compatibility; new code should
    prefer the `psw` property which gives the switching *success* probability
    (P_sw = 1 − SER). Chapter 2 of the thesis reports P_sw throughout, so any
    figure rendered against the experiment is more directly read off `psw`.
    """
    x: np.ndarray
    ser: np.ndarray
    x_label: str

    @property
    def psw(self) -> np.ndarray:
        """Switching success probability  P_sw = 1 − SER."""
        return 1.0 - np.asarray(self.ser)


def _trial_seed(seed: int, i_sot: float, trial_idx: int) -> int:
    """Mix (master seed, I_SOT, trial index) into a deterministic 31-bit seed.

    Identical formula used by both the legacy `np.random.seed` path and the
    modern `np.random.default_rng` path, so the **stream of seeds** is the
    same in both modes — only the underlying generator changes.
    """
    return (int(abs(seed)) * 0x9E3779B1
            + int(round(i_sot * 1e9)) * 0x85EBCA6B
            + int(trial_idx) * 0xC2B2AE35) & 0x7FFFFFFF


def _switching_error_rate_single_isot(
    i_sot: float,
    cfg: SerSotNoVcmaThermalConfig,
    *,
    show_progress: bool = True,
    enable_self_heating: bool = False,
    T_ambient_K: float = 300.0,
    seed: int | None = None,
    rng_mode: str = "legacy",
    integrator: str | None = None,
) -> float:
    """Inner Monte-Carlo: count trials whose final m_z misses `target_mz`
    by more than `failure_tol`.

    `rng_mode` controls how per-trial randomness is reproduced:
      - "legacy"      : sets the global `np.random.seed(...)` once per trial
                        (back-compat behaviour, default).
      - "generator"   : builds a fresh `np.random.default_rng(trial_seed)` per
                        trial and forwards it through to the LLG stepper via
                        the new `rng=` plumbing. Cleaner — does not touch the
                        global numpy state, friendlier to multiprocessing.

    `integrator` lets callers force the LLG step type ("euler_spherical",
    "cayley", or the 2nd-order "cayley_midpoint"); falls back to the stepper
    default when None.
    """
    failures = 0
    loop = trange(cfg.trials, desc=f"MC isot={i_sot:.3e}A", disable=not show_progress)
    common_kw = dict(
        sim_start_step=cfg.sim_start_step,
        sim_mid1_step=cfg.sim_mid1_step,
        sim_mid2_step=cfg.sim_end_step,
        sim_end_step=cfg.sim_end_step,
        pap=cfg.pap,
        v_mtj_stage1=cfg.v_mtj,
        v_mtj_stage2=cfg.v_mtj,
        v_mtj_stage3=cfg.v_mtj,
        i_sot_stage1=i_sot,
        i_sot_stage2=0.0,
        i_sot_stage3=0.0,
        estt_stage1=0, esot_stage1=1,
        estt_stage2=0, esot_stage2=1,
        estt_stage3=0, esot_stage3=1,
        vnv=cfg.vnv, non=cfg.non, r_sot_fl_dl=cfg.r_sot_fl_dl,
        show_progress=False,
        constants=cfg.constants,
        enable_self_heating=enable_self_heating,
        T_ambient_K=T_ambient_K,
    )
    if integrator is not None:
        common_kw["integrator"] = integrator

    for trial_idx in loop:
        rng = None
        if seed is not None:
            ts = _trial_seed(seed, i_sot, trial_idx)
            if rng_mode == "generator":
                rng = np.random.default_rng(ts)
            elif rng_mode == "legacy":
                np.random.seed(ts)
            else:
                raise ValueError(
                    f"rng_mode must be 'legacy' or 'generator', got {rng_mode!r}"
                )
        res = run_piecewise_direct_excitation(rng=rng, **common_kw)
        final_mz = float(res.mz[cfg.sim_end_step])
        if abs(final_mz - cfg.target_mz) > cfg.failure_tol:
            failures += 1

    return failures / float(cfg.trials)


def ser_sot_no_vcma_thermal(
    cfg: SerSotNoVcmaThermalConfig | None = None,
    *,
    show_progress: bool = True,
    enable_self_heating: bool = False,
    T_ambient_K: float = 300.0,
    seed: int | None = None,
    rng_mode: str = "legacy",
    integrator: str | None = None,
) -> SerResult:
    """Monte-Carlo SER vs I_SOT.

    Toggles preserved for back-compat:
      `seed`        — int master seed (None → no seeding).
      `rng_mode`    — "legacy" (default, np.random.seed) or "generator"
                       (np.random.default_rng, no global state).
      `integrator`  — None (stepper default) or "euler_spherical" / "cayley" /
                       "cayley_midpoint" (true 2nd-order implicit midpoint).
    """
    cfg = cfg or SerSotNoVcmaThermalConfig()
    ser = np.array(
        [_switching_error_rate_single_isot(
            i, cfg, show_progress=show_progress,
            enable_self_heating=enable_self_heating, T_ambient_K=T_ambient_K,
            seed=seed, rng_mode=rng_mode, integrator=integrator,
        ) for i in cfg.i_sot_list],
        dtype=float,
    )
    return SerResult(x=np.array(cfg.i_sot_list, dtype=float), ser=ser, x_label=r"$I_{\mathrm{SOT}}$ (A)")


@dataclass
class VariabilitySweepResult:
    """D2D variability MC: wafer-averaged Sigmoid slope vs CV(Δ).

    Wraps `vgsot_sim.analysis.variability` into a self-contained vgsot-sim
    case so the whole §2.3.5 figure can be regenerated from one call.

    Attributes
    ----------
    cv_sweep : np.ndarray
        CV(Δ) values scanned (dimensionless).
    beta_eff : np.ndarray
        Wafer-averaged Sigmoid slope at each CV (V⁻¹).
    vth_eff : np.ndarray
        Wafer-averaged Sigmoid midpoint at each CV (V).
    F_func : np.ndarray
        D2D transfer function F(CV) = β_eff(CV) / β_eff(CV=0).
    eta_c : float
        C2C calibration factor β_meas / β_NB^fit at CV=0.
    beta_combined : np.ndarray
        η_c · F(CV) · β_NB^fit — the joint D2D + C2C prediction.
    """
    cv_sweep: np.ndarray
    beta_eff: np.ndarray
    vth_eff: np.ndarray
    F_func: np.ndarray
    eta_c: float
    beta_combined: np.ndarray
    Delta: float
    Vc0: float
    tw_ns: float


def variability_sweep(
    *,
    Delta: float,
    Vc0: float,
    tw_ns: float = 0.75,
    tau0_ns: float = 1.0,
    beta_meas: float = 44.6,
    cv_sweep=None,
    V_dense=None,
    N_samples: int = 20_000,
    rng=None,
) -> VariabilitySweepResult:
    """Run the Brinkman-decomposed D2D MC sweep and return the joint
    prediction band β^eff(CV_Δ) = η_c · F(CV_Δ) · β_NB^fit.

    Defaults reproduce the chapter §2.3.5 numbers (Device A, P→AP, 0.75 ns).
    """
    from .analysis.variability import transfer_function_F

    if cv_sweep is None:
        cv_sweep = np.array([0.00, 0.03, 0.05, 0.07, 0.0771, 0.10, 0.13, 0.15,
                             0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60])
        cv_sweep = np.sort(np.unique(np.round(cv_sweep, 4)))
    cv_sweep = np.asarray(cv_sweep, dtype=float)
    if rng is None:
        rng = np.random.default_rng(seed=42)

    beta_eff, vth_eff = transfer_function_F(
        cv_sweep, Delta, Vc0,
        tw_ns=tw_ns, tau0_ns=tau0_ns,
        V_dense=V_dense, N_samples=N_samples, rng=rng,
    )
    beta_NB_fit = beta_eff[0]               # CV=0 reference
    F_func = beta_eff / beta_NB_fit
    eta_c  = float(beta_meas / beta_NB_fit)
    beta_combined = eta_c * F_func * beta_NB_fit

    return VariabilitySweepResult(
        cv_sweep=cv_sweep, beta_eff=beta_eff, vth_eff=vth_eff,
        F_func=F_func, eta_c=eta_c, beta_combined=beta_combined,
        Delta=Delta, Vc0=Vc0, tw_ns=tw_ns,
    )


SER_CASES = (
    "ser_sot_no_vcma_thermal",
    "variability_sweep",
)
