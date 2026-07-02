"""Fig 2.21 panels — minimal reservoir validation (serves thesis §2.4.6).

Emits three CLEAN panels (no letters/numbers, English-only) to article/ppt/panels/:

  ch02_21_a.png  per-delay memory capacity MC_k: heterogeneous W_in reservoir vs
                 the broadcast-identical baseline (the i.i.d. collapse);
  ch02_21_b.png  total MC vs node count (filter-bank saturation) + baseline;
  ch02_21_c.png  stochastic binary-device MC vs replicas per node, vs the
                 mean-field limit (the device-averaging cost).

Everything is recomputed live with the SAME fixed seeds as
benchmark_reservoir.py and asserted against the committed reservoir_results.json
(if present) so every plotted number is reproducible.

Run:  PYTHONPATH=src python scripts/10_rtn_reservoir/plot_reservoir_figs.py
"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from vgsot_sim.rtn.reservoir import Reservoir, ReservoirConfig, memory_capacity

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PANELS = REPO / "article" / "ppt" / "panels"
PANELS.mkdir(parents=True, exist_ok=True)

NAVY, TEAL, CRIMSON, GREY = "#1F5FA8", "#1A6B5A", "#A82038", "0.45"
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Arial", "Liberation Serif"],
    "font.size": 11,
    "mathtext.fontset": "stix",
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

cfg = ReservoirConfig()          # tuned defaults (Delta~U(0.5,4), a_in=0.5, dt=1 ns)
ref = None
res_path = HERE / "reservoir_results.json"
if res_path.exists():
    ref = json.loads(res_path.read_text())

# ── recompute with benchmark seeds ───────────────────────────────────────
# baseline: identical nodes, broadcast input (benchmark seed=1)
base = Reservoir(replace(cfg, delta_range=(2.0, 2.0), bias_spread=0.0), seed=1)
base.W_in[:] = 1.0
mc_base = memory_capacity(base, max_delay=40, seed=1)

# heterogeneous reservoir MC vs n (benchmark seed=3)
mc_by_n = {}
for n in (25, 50, 100, 200):
    r = Reservoir(replace(cfg, n_nodes=n), seed=3)
    mc_by_n[n] = memory_capacity(r, max_delay=min(60, n + 5), n_samples=3000, seed=3)

# stochastic device MC vs replicas (benchmark seed=5, n=100)
r_sto = Reservoir(replace(cfg, n_nodes=100), seed=5)
mc_sto = {R: memory_capacity(r_sto, max_delay=40, n_samples=2500, mode="stochastic",
                             n_replicas=R, seed=5).mc_total for R in (1, 4, 16, 64)}

# consistency asserts vs the committed benchmark JSON (same seeds -> identical)
if ref is not None:
    assert abs(mc_base.mc_total - ref["baseline_broadcast_identical"]["mc_total"]) < 1e-9
    for row in ref["mc_vs_n"]:
        assert abs(mc_by_n[row["n"]].mc_total - row["mc_total"]) < 1e-9
    for row in ref["stochastic_mc_vs_replicas"]:
        assert abs(mc_sto[row["replicas"]] - row["mc_total"]) < 1e-9
    print("consistency vs reservoir_results.json: OK")

# ── (a) MC_k vs delay ────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(4.6, 3.4))
het = mc_by_n[100]
ax.plot(het.delays, het.mc_k, "-o", ms=3.5, color=NAVY, lw=1.5)
ax.plot(mc_base.delays, mc_base.mc_k, "-s", ms=3.5, color=GREY, lw=1.3)
ax.text(11, 0.44, rf"heterogeneous $W_{{\rm in}}$ ($n=100$), MC$={het.mc_total:.1f}$",
        color=NAVY, fontsize=9.5)
ax.text(9, 0.10, rf"broadcast identical nodes, MC$={mc_base.mc_total:.2f}$",
        color=GREY, fontsize=9.5)
ax.set_xlim(0, 42)
ax.set_ylim(-0.03, 1.05)
ax.set_xlabel(r"delay $k$")
ax.set_ylabel(r"MC$_k=\mathrm{corr}^2(\hat u[t-k],\,u[t-k])$")
ax.set_title("Per-delay memory capacity", fontsize=11.5, pad=10)
fig.savefig(PANELS / "ch02_21_a.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ── (b) MC vs n ──────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(4.6, 3.4))
ns = sorted(mc_by_n)
totals = [mc_by_n[n].mc_total for n in ns]
ax.plot(ns, totals, "-o", ms=6, color=NAVY, lw=1.6)
for n, v in zip(ns, totals):
    ax.annotate(f"{v:.1f}", xy=(n, v), xytext=(0, 7), textcoords="offset points",
                ha="center", fontsize=9, color=NAVY)
ax.axhline(mc_base.mc_total, color=GREY, lw=1.2, ls="--")
ax.text(198, mc_base.mc_total + 0.25, "broadcast baseline", color=GREY,
        fontsize=9.5, ha="right")
ax.set_xlim(0, 215)
ax.set_ylim(0, 10)
ax.set_xlabel(r"number of nodes $n$")
ax.set_ylabel("total MC")
ax.set_title("Capacity vs node count (filter-bank limit)", fontsize=11.5, pad=10)
fig.savefig(PANELS / "ch02_21_b.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ── (c) stochastic devices: MC vs replicas ───────────────────────────────
fig, ax = plt.subplots(figsize=(4.6, 3.4))
reps = sorted(mc_sto)
vals = [mc_sto[R] for R in reps]
ax.semilogx(reps, vals, "-o", ms=6, color=CRIMSON, lw=1.6, base=2)
for R, v in zip(reps, vals):
    ax.annotate(f"{v:.2f}", xy=(R, v), xytext=(0, 7), textcoords="offset points",
                ha="center", fontsize=9, color=CRIMSON)
mf = mc_by_n[100].mc_total
ax.axhline(mf, color=NAVY, lw=1.2, ls="--")
ax.text(1.05, mf - 0.55, rf"mean-field limit (MC$={mf:.1f}$)", color=NAVY, fontsize=9.5)
ax.set_xticks(reps, [str(R) for R in reps])
ax.set_ylim(0, 9)
ax.set_xlabel(r"devices per node $R$ (binary states, averaged)")
ax.set_ylabel("total MC")
ax.set_title(r"Device-averaging cost ($n=100$)", fontsize=11.5, pad=10)
fig.savefig(PANELS / "ch02_21_c.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"wrote 3 panels -> {PANELS}")
