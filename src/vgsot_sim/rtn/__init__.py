"""Random-telegraph-noise (RTN) reservoir-computing primitive for vgsot-sim.

A continuous-time two-state Markov model of a free-running, low-barrier
superparamagnetic sMTJ. Complements the pulse-switching LLG solver (it does not
replace it): the LLG kernels answer "did this write pulse switch the cell?",
while this module answers "how does the cell hop between states in time under a
sustained bias?" — the device primitive reservoir computing feeds on.
"""
from __future__ import annotations

from .telegraph import (
    TelegraphArray,
    TelegraphParams,
    neel_brown_rate,
    relaxation_time,
    simulate_trace,
    stationary_mean,
    tau_max,
    up_down_rates,
)

__all__ = [
    "TelegraphParams",
    "TelegraphArray",
    "neel_brown_rate",
    "up_down_rates",
    "stationary_mean",
    "relaxation_time",
    "tau_max",
    "simulate_trace",
]
