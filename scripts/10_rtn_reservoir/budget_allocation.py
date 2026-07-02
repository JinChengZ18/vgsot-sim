"""D3 — device-budget allocation: given B = n x R physical sMTJs, how to split
logical nodes (n) vs devices-averaged-per-node (R) for maximum memory capacity?

Stochastic (binary-device) mode of the heterogeneous W_in reservoir; ties the
reservoir story to the thesis §2.3.2 sampling-budget analysis. Writes
budget_allocation.json next to this script.

Run:  PYTHONPATH=src python scripts/10_rtn_reservoir/budget_allocation.py
"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from vgsot_sim.rtn.reservoir import Reservoir, ReservoirConfig, memory_capacity

HERE = Path(__file__).resolve().parent
OUT = HERE / "budget_allocation.json"
cfg = ReservoirConfig()

res = {"budgets": []}
for B in (64, 256, 1024):
    row = {"B": B, "allocations": []}
    print(f"budget B = {B} devices:")
    n = 8
    while n <= B:
        R = B // n
        r = Reservoir(replace(cfg, n_nodes=n), seed=7)
        mc = memory_capacity(r, max_delay=30, n_samples=2000, washout=200,
                             mode="stochastic", n_replicas=R, seed=7)
        print(f"  n={n:5d} x R={R:4d}: MC = {mc.mc_total:5.2f}")
        row["allocations"].append({"n": n, "R": R, "mc": mc.mc_total})
        n *= 2
    best = max(row["allocations"], key=lambda a: a["mc"])
    row["best"] = best
    print(f"  -> best: n={best['n']}, R={best['R']}, MC={best['mc']:.2f}")
    res["budgets"].append(row)

OUT.write_text(json.dumps(res, indent=2))
print(f"wrote {OUT}")
