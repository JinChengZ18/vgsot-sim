from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np
from tqdm import trange

from . import constants as C
from .configs import SerOptimizedVgsotConfig, SerSotNoVcmaThermalConfig
from .time_series_cases import run_piecewise_direct_excitation, run_two_pulse_optimized


@dataclass
class SerResult:
    x: np.ndarray
    ser: np.ndarray
    x_label: str


@dataclass
class SerOptimizedResult:
    t1_s: np.ndarray
    ser: np.ndarray
    mz_at_t1_avg: np.ndarray


def _switching_error_rate_single_isot(
    i_sot: float,
    cfg: SerSotNoVcmaThermalConfig,
    *,
    show_progress: bool = True,
) -> float:
    failures = 0
    loop = trange(cfg.trials, desc=f"MC isot={i_sot:.3e}A", disable=not show_progress)
    for _ in loop:
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
            i_sot_stage2=0.0,
            i_sot_stage3=0.0,
            estt_stage1=0,
            esot_stage1=1,
            estt_stage2=0,
            esot_stage2=1,
            estt_stage3=0,
            esot_stage3=1,
            vnv=cfg.vnv,
            non=cfg.non,
            r_sot_fl_dl=cfg.r_sot_fl_dl,
            show_progress=False,
        )
        final_mz = float(res.mz[cfg.sim_end_step])
        if abs(final_mz - cfg.target_mz) > cfg.failure_tol:
            failures += 1

    return failures / float(cfg.trials)


def ser_sot_no_vcma_thermal(
    cfg: SerSotNoVcmaThermalConfig | None = None,
    *,
    show_progress: bool = True,
) -> SerResult:
    cfg = cfg or SerSotNoVcmaThermalConfig()
    ser = np.array(
        [_switching_error_rate_single_isot(i, cfg, show_progress=show_progress) for i in cfg.i_sot_list],
        dtype=float,
    )
    return SerResult(x=np.array(cfg.i_sot_list, dtype=float), ser=ser, x_label=r"$I_{\mathrm{SOT}}$ (A)")


def ser_optimized_vgsot(
    cfg: SerOptimizedVgsotConfig | None = None,
    *,
    show_progress: bool = True,
) -> SerOptimizedResult:
    cfg = cfg or SerOptimizedVgsotConfig()

    i_sot = cfg.i_sot
    if i_sot is None:
        i_sot = (2 * C.e * C.u0 * C.Ms * C.tf * C.A2 * (-50 * 1000 / (4 * C.pi))) / (C.h_bar * C.theta_SH)

    sim_end_step = int(cfg.sim_total_time_s / C.t_step)

    ser_list: List[float] = []
    mz_avg_list: List[float] = []

    for t1_s in cfg.t1_list_s:
        t2_s = cfg.total_pulse_s - t1_s
        failures = 0
        mz_sum = 0.0

        loop = trange(cfg.iterations_num, desc=f"MC t1={t1_s:.3e}s", disable=not show_progress)
        for _ in loop:
            res = run_two_pulse_optimized(
                t1_s=t1_s,
                t2_s=t2_s,
                v_mtj_1=cfg.v_mtj_1,
                v_mtj_2=cfg.v_mtj_2,
                i_sot_1=i_sot,
                i_sot_2=0.0,
                sim_total_time_s=cfg.sim_total_time_s,
                pap=cfg.pap,
                non=cfg.non,
                vnv=cfg.vnv,
                r_sot_fl_dl=cfg.r_sot_fl_dl,
                show_progress=False,
            )

            idx_t1 = min(max(int(t1_s / C.t_step), 0), len(res.mz) - 1)
            mz_sum += float(res.mz[idx_t1])

            final_mz = float(res.mz[sim_end_step])
            if abs(final_mz - cfg.target_final_mz) > cfg.failure_tol:
                failures += 1

        ser_list.append(failures / float(cfg.iterations_num))
        mz_avg_list.append(mz_sum / float(cfg.iterations_num))

    return SerOptimizedResult(
        t1_s=np.array(cfg.t1_list_s, dtype=float),
        ser=np.array(ser_list, dtype=float),
        mz_at_t1_avg=np.array(mz_avg_list, dtype=float),
    )


SER_CASES = (
    "ser_sot_no_vcma_thermal",
    "ser_optimized_vgsot",
)