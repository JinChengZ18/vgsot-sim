from __future__ import annotations

from .ser_cases import SER_CASES, SerResult, ser_sot_no_vcma_thermal
from .time_series_cases import (
    TIME_SERIES_CASES,
    SimResult,
    SweepResult,
    run_piecewise_direct_excitation,
    run_piecewise_terminal_voltage,
    sot_only_constant_current,
    sot_switching_no_vcma,
    terminal_voltage_control,
)

ALL_CASES = TIME_SERIES_CASES + SER_CASES

__all__ = [
    "SimResult",
    "SweepResult",
    "SerResult",
    "run_piecewise_terminal_voltage",
    "run_piecewise_direct_excitation",
    "terminal_voltage_control",
    "sot_only_constant_current",
    "sot_switching_no_vcma",
    "ser_sot_no_vcma_thermal",
    "ALL_CASES",
]
