from __future__ import annotations

__version__ = "0.1.0"

from . import anisotropy, constants, dynamic_switching, electronic, initialize, stochastic, tmr
from .configs import (
    OptimizedVgsotSwitchingConfig,
    SerOptimizedVgsotConfig,
    SerSotNoVcmaThermalConfig,
    SotOnlyConstantCurrentConfig,
    SotSwitchingNoVcmaConfig,
    TerminalVoltageControlConfig,
    VcmaAssistedSwitchingIsotSweepConfig,
    VcmaAssistedSwitchingVmtjSweepConfig,
)
from .cases import (
    ALL_CASES,
    SerOptimizedResult,
    SerResult,
    SimResult,
    SweepResult,
    optimized_vgsot_switching,
    run_piecewise_direct_excitation,
    run_piecewise_terminal_voltage,
    run_two_pulse_optimized,
    ser_optimized_vgsot,
    ser_sot_no_vcma_thermal,
    sot_only_constant_current,
    sot_switching_no_vcma,
    terminal_voltage_control,
    vcma_assisted_switching_isot_sweep,
    vcma_assisted_switching_vmtj_sweep,
)
from .result_io import (
    build_stem,
    config_to_params,
    ensure_result_dir,
    save_grouped_timeseries_csv,
    save_single_plot,
    save_three_panel_plot,
    save_timeseries_csv,
    save_two_panel_plot,
    save_xy_csv,
)

__all__ = [
    "__version__",
    "constants", "stochastic", "initialize", "electronic", "anisotropy", "dynamic_switching", "tmr",
    "TerminalVoltageControlConfig", "SotOnlyConstantCurrentConfig", "SotSwitchingNoVcmaConfig", "SerSotNoVcmaThermalConfig",
    "VcmaAssistedSwitchingIsotSweepConfig", "VcmaAssistedSwitchingVmtjSweepConfig", "OptimizedVgsotSwitchingConfig", "SerOptimizedVgsotConfig",
    "SimResult", "SweepResult", "SerResult", "SerOptimizedResult",
    "terminal_voltage_control", "sot_only_constant_current", "sot_switching_no_vcma", "ser_sot_no_vcma_thermal",
    "vcma_assisted_switching_isot_sweep", "vcma_assisted_switching_vmtj_sweep", "optimized_vgsot_switching", "ser_optimized_vgsot",
    "run_piecewise_terminal_voltage", "run_piecewise_direct_excitation", "run_two_pulse_optimized",
    "ensure_result_dir", "build_stem", "config_to_params", "save_grouped_timeseries_csv", "save_timeseries_csv", "save_xy_csv",
    "save_single_plot", "save_two_panel_plot", "save_three_panel_plot", "ALL_CASES",
]
