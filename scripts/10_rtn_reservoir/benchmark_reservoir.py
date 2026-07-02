"""Stage-3 reservoir benchmarks: memory capacity, NARMA-10, baselines, device-noise penalty.

Run: PYTHONPATH=src python scripts/10_rtn_reservoir/benchmark_reservoir.py
Writes reservoir_results.json next to this script. Fast (mean-field path is vectorised).
See this directory's README.
"""
from __future__ import annotations

import json
import os
from dataclasses import replace

import numpy as np

from vgsot_sim.rtn.reservoir import (
    DelayReservoir,
    Reservoir,
    ReservoirConfig,
    information_processing_capacity,
    memory_capacity,
    narma10_task,
    parity_task,
)

OUT = os.path.join(os.path.dirname(__file__), "reservoir_results.json")
cfg = ReservoirConfig()              # tuned defaults: Delta~U(0.5,4), a_in=0.5, dt=1 ns
res = {}

# --- A2 baseline: identical nodes, broadcast same input (the audit's collapse) ---
base = Reservoir(replace(cfg, delta_range=(2.0, 2.0), bias_spread=0.0), seed=1)
base.W_in[:] = 1.0
mc_base = memory_capacity(base, max_delay=40, seed=1)
res["baseline_broadcast_identical"] = {"n": base.n, "mc_total": mc_base.mc_total}
print(f"A2 baseline (identical + broadcast): MC = {mc_base.mc_total:.2f}")

# --- heterogeneous W_in reservoir: MC vs n_nodes ---
print("\nheterogeneous reservoir — memory capacity vs n:")
mc_vs_n = []
for n in (25, 50, 100, 200):
    r = Reservoir(replace(cfg, n_nodes=n), seed=3)
    mc = memory_capacity(r, max_delay=min(60, n + 5), n_samples=3000, seed=3)
    print(f"  n={n:3d}  MC = {mc.mc_total:.2f}")
    mc_vs_n.append({"n": n, "mc_total": mc.mc_total, "mc_k": mc.mc_k.tolist()})
res["mc_vs_n"] = mc_vs_n

# --- NARMA-10 ---
print("\nNARMA-10 (n=200):")
narma = []
r = Reservoir(replace(cfg, n_nodes=200), seed=4)
for alpha in (1e-8, 1e-6, 1e-4):
    tr = narma10_task(r, n_samples=4000, alpha=alpha, seed=4)
    print(f"  alpha={alpha:.0e}  NRMSE={tr.nrmse:.3f}  R2={tr.r2:.3f}")
    narma.append({"alpha": alpha, "nrmse": tr.nrmse, "r2": tr.r2})
res["narma10_n200"] = narma

# --- device stochasticity penalty: stochastic binary states vs replicas ---
print("\nstochastic (binary device) MC vs replicas (n=100):")
sto = []
r = Reservoir(replace(cfg, n_nodes=100), seed=5)
for R in (1, 4, 16, 64):
    mc = memory_capacity(r, max_delay=40, n_samples=2500, mode="stochastic", n_replicas=R, seed=5)
    print(f"  replicas={R:3d}  MC = {mc.mc_total:.2f}")
    sto.append({"replicas": R, "mc_total": mc.mc_total})
res["stochastic_mc_vs_replicas"] = sto

# --- IPC + parity, and the delay-reservoir comparison ---
print("\nIPC + parity-N (n=100):")
r = Reservoir(replace(cfg, n_nodes=100), seed=6)
ipc = information_processing_capacity(r, seed=6)
par = {N: parity_task(r, n_bits=N, seed=6) for N in (2, 3, 4)}
print(f"  filter-bank: IPC={ipc.total:.2f} (deg1={ipc.deg1:.2f} deg2={ipc.deg2:.2f})  "
      f"parity2/3/4={par[2]:.2f}/{par[3]:.2f}/{par[4]:.2f}")
res["filterbank_ipc_parity"] = {"ipc_total": ipc.total, "ipc_deg1": ipc.deg1,
                                "ipc_deg2": ipc.deg2, "parity": par}

dr = DelayReservoir(n_virtual=100, Delta=2.0, fb_gain=0.3, a_in=0.5, bias_spread=0.2, seed=7)
mc_dr = memory_capacity(dr, max_delay=60, n_samples=3000, seed=7)
ipc_dr = information_processing_capacity(dr, seed=7)
par_dr = {N: parity_task(dr, n_bits=N, seed=7) for N in (2, 3, 4)}
print(f"  delay (Appeltant): MC={mc_dr.mc_total:.2f}  IPC={ipc_dr.total:.2f} "
      f"(deg1={ipc_dr.deg1:.2f} deg2={ipc_dr.deg2:.2f})  "
      f"parity2/3/4={par_dr[2]:.2f}/{par_dr[3]:.2f}/{par_dr[4]:.2f}")
print("  -> delay trades linear memory for nonlinear capacity; does NOT beat the MC ceiling")
res["delay_reservoir"] = {"mc_total": mc_dr.mc_total, "ipc_total": ipc_dr.total,
                          "ipc_deg1": ipc_dr.deg1, "ipc_deg2": ipc_dr.deg2, "parity": par_dr}

with open(OUT, "w") as f:
    json.dump(res, f, indent=2)
print(f"\nWrote {OUT}")
