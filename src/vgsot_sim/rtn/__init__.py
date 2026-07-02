"""Random-telegraph-noise (RTN) reservoir-computing primitive for vgsot-sim.

A continuous-time two-state Markov model of a free-running, low-barrier
superparamagnetic sMTJ. Complements the pulse-switching LLG solver (it does not
replace it): the LLG kernels answer "did this write pulse switch the cell?",
while this module answers "how does the cell hop between states in time under a
sustained bias?" — a candidate single-**node** primitive for reservoir computing.

This is one node, NOT a reservoir: there are no input weights, no node coupling,
and no trained readout here. ``V`` is a phenomenological effective bias/tilt valid
for ``|V| < Vc0``, not literally ``V_MTJ``. Building an actual reservoir and
calibrating ``V`` to a device drive are in ``scripts/10_rtn_reservoir/`` (see its README).
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
from .reservoir import (
    DelayReservoir,
    Reservoir,
    ReservoirConfig,
    information_processing_capacity,
    memory_capacity,
    narma10,
    narma10_task,
    parity_task,
    ridge_fit,
    ridge_predict,
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
    # reservoir layer (stage 3)
    "Reservoir",
    "ReservoirConfig",
    "DelayReservoir",
    "ridge_fit",
    "ridge_predict",
    "memory_capacity",
    "narma10",
    "narma10_task",
    "parity_task",
    "information_processing_capacity",
]
