from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple


@dataclass
class TerminalVoltageControlConfig:
    sim_start_step: int = 1
    sim_mid1_step: int = 2000
    sim_mid2_step: Optional[int] = None
    sim_end_step: int = 5000
    pap: int = 1

    v_stage1: Tuple[float, float, float] = (1.0, 0.0, 0.1)
    v_stage2: Tuple[float, float, float] = (-1.0, 0.0, 0.0)
    v_stage3: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    estt_stage1: int = 0
    esot_stage1: int = 1
    estt_stage2: int = 0
    esot_stage2: int = 1
    estt_stage3: int = 1
    esot_stage3: int = 1

    vnv: int = 1
    non: int = 1
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 5e-10


@dataclass
class SotOnlyConstantCurrentConfig:
    sim_start_step: int = 1
    sim_mid1_step: int = 2000
    sim_end_step: int = 5000
    pap: int = 1

    i_sot_stage1: float = -95e-6
    i_sot_stage2: float = 0.0
    v_mtj_stage1: float = 0.0
    v_mtj_stage2: float = 0.0

    vnv: int = 1
    non: int = 1
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 5e-10


@dataclass
class SotSwitchingNoVcmaConfig:
    i_sot_list: Sequence[float] = (-85e-6, -90e-6, -95e-6, -100e-6)
    sim_start_step: int = 1
    sim_mid1_step: int = 2000
    sim_end_step: int = 5000
    pap: int = 1

    non: int = 1
    v_mtj: float = 0.0
    i_sot_relax: float = 0.0
    vnv: int = 0
    r_sot_fl_dl: float = 0.0
    tick_spacing_s: float = 5e-10


@dataclass
class SerSotNoVcmaThermalConfig:
    i_sot_list: Sequence[float] = (-100e-6, -98e-6, -96e-6, -94e-6, -92e-6, -90e-6)
    trials: int = 1000
    sim_start_step: int = 1
    sim_mid1_step: int = 2000
    sim_end_step: int = 5000
    pap: int = 1

    non: int = 1
    v_mtj: float = 0.0
    vnv: int = 0
    r_sot_fl_dl: float = 0.0
    target_mz: float = -1.0
    failure_tol: float = 1e-1


@dataclass
class VcmaAssistedSwitchingIsotSweepConfig:
    v_mtj: float = 1.2
    i_sot_list: Sequence[float] = (-90e-6, -30e-6, -18e-6, -16e-6)
    sim_start_step: int = 1
    sim_end_step: int = 25000
    pap: int = 1

    non: int = 1
    vnv: int = 1
    r_sot_fl_dl: float = 0.0
    tick_spacing_s: float = 5e-9


@dataclass
class VcmaAssistedSwitchingVmtjSweepConfig:
    i_sot: Optional[float] = None
    v_mtj_list: Sequence[float] = (1.3189, 1.3191, 1.333, 1.3489, 1.4937)
    sim_start_step: int = 1
    sim_end_step: int = 25000
    pap: int = 1

    non: int = 1
    vnv: int = 1
    r_sot_fl_dl: float = 0.0
    tick_spacing_s: float = 5e-9


@dataclass
class OptimizedVgsotSwitchingConfig:
    v_mtj_1: float = 1.4937
    v_mtj_2: float = -1.0
    i_sot: Optional[float] = None
    t_pairs_s: Sequence[Tuple[float, float]] = (
        (25e-9, 0.0),
        (1.4e-9, 1.6e-9),
        (1.8e-9, 1.2e-9),
        (2.2e-9, 0.8e-9),
    )
    sim_total_time_s: float = 25e-9
    pap: int = 1

    non: int = 1
    vnv: int = 1
    r_sot_fl_dl: float = 0.0
    tick_spacing_s: float = 5e-9


@dataclass
class SerOptimizedVgsotConfig:
    iterations_num: int = 300
    v_mtj_1: float = 1.4937
    v_mtj_2: float = -1.0
    i_sot: Optional[float] = None
    t1_list_s: Sequence[float] = (1.3e-9, 1.4e-9, 1.5e-9, 1.6e-9, 1.7e-9, 1.8e-9, 1.9e-9)
    total_pulse_s: float = 3e-9
    sim_total_time_s: float = 25e-9
    pap: int = 1

    non: int = 1
    vnv: int = 1
    r_sot_fl_dl: float = 0.0
    target_final_mz: float = 1.0
    failure_tol: float = 1e-1



