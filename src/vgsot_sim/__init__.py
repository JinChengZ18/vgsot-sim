from __future__ import annotations

__version__ = "0.1.1"

from . import anisotropy, constants, dynamic_switching, electronic, initialize, rtn, stochastic, tmr
from .configs import (
    PhysicalConstantsConfig,
    SerSotNoVcmaThermalConfig,
    SotOnlyConstantCurrentConfig,
    SotSwitchingNoVcmaConfig,
    TerminalVoltageControlConfig,
)
from .cases import (
    ALL_CASES,
    SerResult,
    SimResult,
    SweepResult,
    run_piecewise_direct_excitation,
    run_piecewise_terminal_voltage,
    ser_sot_no_vcma_thermal,
    sot_only_constant_current,
    sot_switching_no_vcma,
    terminal_voltage_control,
)
from .constants import get_constants, reset_constants, set_constants
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
    # Physics modules
    "constants", "stochastic", "initialize", "electronic", "anisotropy", "dynamic_switching", "tmr",
    "rtn",
    # Config dataclasses (one per case + the shared physics constants)
    "PhysicalConstantsConfig",
    "TerminalVoltageControlConfig",
    "SotOnlyConstantCurrentConfig",
    "SotSwitchingNoVcmaConfig",
    "SerSotNoVcmaThermalConfig",
    # Result containers
    "SimResult", "SweepResult", "SerResult",
    # High-level cases
    "terminal_voltage_control",
    "sot_only_constant_current",
    "sot_switching_no_vcma",
    "ser_sot_no_vcma_thermal",
    # Low-level kernels
    "run_piecewise_terminal_voltage",
    "run_piecewise_direct_excitation",
    # IO helpers
    "ensure_result_dir", "build_stem", "config_to_params",
    "save_grouped_timeseries_csv", "save_timeseries_csv", "save_xy_csv",
    "save_single_plot", "save_two_panel_plot", "save_three_panel_plot",
    # Misc
    "ALL_CASES",
    "get_constants", "set_constants", "reset_constants",
]
