"""Fig 2.20 panels — sLLG bridge validation (serves thesis §2.4.5).

Emits three CLEAN panels (no letters/numbers, English-only) to article/ppt/panels/:

  ch02_20_a.png  low-barrier free-running m_z(t) trace (fresh short run, fixed seed)
                 with the hysteretic dwell thresholds;
  ch02_20_b.png  dwell-time survival (CCDF, semilog-y) from the committed
                 bridge_results.json vs the exponential law — the CV~1 signature;
  ch02_20_c.png  <m_z> vs dimensionless tilt: sLLG points, least-squares A*tanh fit,
                 ideal tanh for contrast (amplitude compression).

(b)/(c) consume scripts/10_rtn_reservoir/bridge_results.json (regenerate with
validate_bridge.py — fixed seeds, deterministic). Letters + export via
scripts/build_ppt_figs.py.

Run:  PYTHONPATH=src python scripts/10_rtn_reservoir/plot_bridge_figs.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.rtn import bridge

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PANELS = REPO / "article" / "ppt" / "panels"
PANELS.mkdir(parents=True, exist_ok=True)
RESULTS = HERE / "bridge_results.json"
if not RESULTS.exists():
    raise SystemExit("bridge_results.json missing - run validate_bridge.py first "
                     "(fixed seeds, ~20 min).")
res = json.loads(RESULTS.read_text())

NAVY, TEAL, CRIMSON = "#1F5FA8", "#1A6B5A", "#A82038"
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Arial", "Liberation Serif"],
    "font.size": 11,
    "mathtext.fontset": "stix",
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

DT = 4e-12
c0 = PhysicalConstantsConfig()
c0 = type(c0)(**{**c0.__dict__, "t_step": DT})

# ── (a) free-running m_z(t) trace, Delta = 2 ─────────────────────────────
c = bridge.low_barrier_constants(2.0, c0)
n_steps = 500_000                                  # 2 us
mz = bridge.free_run(c, n_steps, h_ex_z=0.0, seed=42)
t_ns = np.arange(n_steps) * DT * 1e9
fig, ax = plt.subplots(figsize=(4.6, 3.4))
ax.plot(t_ns[::5], mz[::5], color=NAVY, lw=0.4)
for thr in (0.5, -0.5):
    ax.axhline(thr, color=CRIMSON, lw=0.9, ls="--")
ax.text(t_ns[-1] * 0.995, 0.56, r"$+m_z^{\rm thr}$", ha="right", fontsize=9, color=CRIMSON)
ax.text(t_ns[-1] * 0.995, -0.72, r"$-m_z^{\rm thr}$", ha="right", fontsize=9, color=CRIMSON)
ax.set_xlim(0, t_ns[-1])
ax.set_ylim(-1.25, 1.25)
ax.set_xlabel(r"$t$ (ns)")
ax.set_ylabel(r"$m_z$")
ax.set_title(r"Free-running sLLG ($\Delta=2$, zero bias)", fontsize=11.5, pad=10)
fig.savefig(PANELS / "ch02_20_a.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ── (b) dwell-time survival vs exponential ───────────────────────────────
dwell = np.array(res["main_delta2"]["dwell_ns"], dtype=float)
tau_bar = res["main_delta2"]["tau_dwell_ns"]
cv = res["main_delta2"]["dwell_cv"]
d_sorted = np.sort(dwell)
surv = 1.0 - np.arange(1, len(d_sorted) + 1) / len(d_sorted)
surv = np.clip(surv, 1.0 / len(d_sorted) / 2, None)   # keep last point plottable
fig, ax = plt.subplots(figsize=(4.6, 3.4))
ax.semilogy(d_sorted, surv, drawstyle="steps-post", color=NAVY, lw=1.4)
tt = np.linspace(0, d_sorted.max() * 1.05, 200)
ax.semilogy(tt, np.exp(-tt / tau_bar), color=CRIMSON, lw=1.4, ls="--")
ax.text(0.97, 0.93, rf"$N={len(dwell)}$", transform=ax.transAxes, ha="right", fontsize=10)
ax.text(0.97, 0.84, rf"$\bar\tau_{{\rm dwell}}={tau_bar:.0f}$ ns",
        transform=ax.transAxes, ha="right", fontsize=10)
ax.text(0.97, 0.75, rf"$\mathrm{{CV}}={cv:.2f}$", transform=ax.transAxes, ha="right", fontsize=10)
ax.text(150, np.exp(-150 / tau_bar) * 1.6, r"$e^{-t/\bar\tau}$", color=CRIMSON, fontsize=10)
ax.set_xlabel(r"dwell time $t$ (ns)")
ax.set_ylabel(r"survival $P(\tau_{\rm dwell}>t)$")
ax.set_title(r"Dwell-time statistics ($\Delta=2$, 8 $\mu$s)", fontsize=11.5, pad=10)
fig.savefig(PANELS / "ch02_20_b.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ── (c) <m_z> vs tilt: A*tanh fit ────────────────────────────────────────
tilt = np.array(res["tanh_calib"]["tilt"], dtype=float)
mzm = np.array(res["tanh_calib"]["mz_mean"], dtype=float)
th = np.tanh(tilt)
A = float(np.sum(mzm * th) / np.sum(th * th))        # least-squares amplitude
xx = np.linspace(-2.6, 2.6, 300)
fig, ax = plt.subplots(figsize=(4.6, 3.4))
ax.plot(xx, np.tanh(xx), color="0.55", lw=1.2, ls="--")
ax.plot(xx, A * np.tanh(xx), color=CRIMSON, lw=1.6)
ax.plot(tilt, mzm, "o", ms=7, mec=NAVY, mfc="white", mew=1.6, zorder=3)
ax.text(1.15, 0.99, r"ideal $\tanh$", color="0.45", fontsize=10)
ax.text(1.5, 0.44, rf"$A\tanh(\cdot)$, $A={A:.2f}$", color=CRIMSON, fontsize=10)
ax.text(-2.45, 0.78, "sLLG time averages", color=NAVY, fontsize=10)
ax.axhline(0.0, color="0.85", lw=0.7, zorder=0)
ax.axvline(0.0, color="0.85", lw=0.7, zorder=0)
ax.set_xlabel(r"dimensionless tilt $\mu_0 M_s V_{\rm mag} h_z / (k_B T)$")
ax.set_ylabel(r"$\langle m_z\rangle$")
ax.set_ylim(-1.1, 1.1)
ax.set_title(r"Transfer function ($\Delta=2$)", fontsize=11.5, pad=10)
fig.savefig(PANELS / "ch02_20_c.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"wrote 3 panels -> {PANELS}   (fitted A = {A:.3f})")
