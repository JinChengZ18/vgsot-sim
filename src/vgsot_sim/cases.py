from __future__ import annotations

from .ser_cases import SER_CASES, SerOptimizedResult, SerResult, ser_optimized_vgsot, ser_sot_no_vcma_thermal
from .time_series_cases import (
    TIME_SERIES_CASES,
    SimResult,
    SweepResult,
    optimized_vgsot_switching,
    run_piecewise_direct_excitation,
    run_piecewise_terminal_voltage,
    run_two_pulse_optimized,
    sot_only_constant_current,
    sot_switching_no_vcma,
    terminal_voltage_control,
    vcma_assisted_switching_isot_sweep,
    vcma_assisted_switching_vmtj_sweep,
)

ALL_CASES = TIME_SERIES_CASES + SER_CASES

__all__ = [
    "SimResult",
    "SweepResult",
    "SerResult",
    "SerOptimizedResult",
    "run_piecewise_terminal_voltage",
    "run_piecewise_direct_excitation",
    "run_two_pulse_optimized",
    "terminal_voltage_control",
    "sot_only_constant_current",
    "sot_switching_no_vcma",
    "vcma_assisted_switching_isot_sweep",
    "vcma_assisted_switching_vmtj_sweep",
    "optimized_vgsot_switching",
    "ser_sot_no_vcma_thermal",
    "ser_optimized_vgsot",
    "ALL_CASES",
]

