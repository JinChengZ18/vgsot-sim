"""Fig 2.19 panels — RTN node physics (serves thesis §2.4.2–2.4.4).

Emits three CLEAN panels (no panel letters, no figure numbers, English-only,
mathtext subscripts) to article/ppt/panels/:

  ch02_19_a.png  stationary mean <s> = tanh(Delta V / Vc0) for a Delta family,
                 with the physical domain |V| < Vc0 marked;
  ch02_19_b.png  relaxation time tau(V) (semilog) for the same family — the
                 memory/nonlinearity trade-off; tau -> tau0 asymptote;
  ch02_19_c.png  sample RTN traces at three bias points (Delta = 3.8).

Panel letters + export to article/figs/ are added by scripts/build_ppt_figs.py
(the PPT assembly flow); this script never writes into article/figs/.

Run:  PYTHONPATH=src python scripts/10_rtn_reservoir/plot_node_figs.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from vgsot_sim.rtn.telegraph import relaxation_time, simulate_trace, stationary_mean, TelegraphParams

REPO = Path(__file__).resolve().parents[2]
PANELS = REPO / "article" / "ppt" / "panels"
PANELS.mkdir(parents=True, exist_ok=True)

# Palette consistent with demo/plot_t_circuit.py
THU_DEEP, NAVY, TEAL, CRIMSON, AMBER = "#660874", "#1F5FA8", "#1A6B5A", "#A82038", "#C47A00"
CHARCOAL = "#2B2B2B"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Arial", "Liberation Serif"],
    "font.size": 11,
    "mathtext.fontset": "stix",
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

VC0 = 0.884
TAU0 = 1.0
DELTAS = [(2.0, NAVY), (3.8, TEAL), (5.15, CRIMSON)]
V = np.linspace(-1.15, 1.15, 601)


def _shade_domain(ax, ymin, ymax):
    """Grey out |V| > Vc0 (outside the two-state model's physical domain)."""
    for sgn in (-1, 1):
        ax.axvspan(sgn * VC0, sgn * 1.15, color="0.88", zorder=0)
        ax.axvline(sgn * VC0, color="0.45", lw=1.0, ls="--", zorder=1)
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(ymin, ymax)


# ── (a) stationary mean ──────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(4.6, 3.4))
_shade_domain(ax, -1.12, 1.12)
for d, c in DELTAS:
    ax.plot(V, stationary_mean(V, Delta=d, Vc0=VC0), color=c, lw=1.8)
# direct labels stacked in the empty lower-right quadrant, thin leaders to curves
anchors = {2.0: 0.42, 3.8: 0.30, 5.15: 0.20}      # leader target V per curve
for k, (d, c) in enumerate(DELTAS):
    xa = anchors[d]
    ya = float(stationary_mean(xa, Delta=d, Vc0=VC0))
    ax.annotate(rf"$\Delta={d:g}$", xy=(xa, ya), xytext=(0.62, -0.18 - 0.30 * k),
                color=c, fontsize=10,
                arrowprops=dict(arrowstyle="-", color=c, lw=0.7, shrinkB=2))
ax.text(VC0 + 0.10, -0.42, r"$+V_{c0}$", fontsize=9, color="0.35", rotation=90)
ax.text(-VC0 - 0.22, -0.42, r"$-V_{c0}$", fontsize=9, color="0.35", rotation=90)
ax.axhline(0.0, color="0.8", lw=0.7, zorder=0)
ax.set_xlabel(r"$V$ (V)")
ax.set_ylabel(r"$\langle s\rangle_{\infty}$")
ax.set_title("Bias-dependent stationary mean", fontsize=11.5, pad=10)
fig.savefig(PANELS / "ch02_19_a.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ── (b) relaxation time ──────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(4.6, 3.4))
for d, c in DELTAS:
    ax.semilogy(V, relaxation_time(V, tau0=TAU0, Delta=d, Vc0=VC0), color=c, lw=1.8)
_shade_domain(ax, 0.2, 400.0)
ax.axhline(TAU0, color="0.45", lw=1.0, ls=":")
ax.text(-1.10, TAU0 * 1.25, r"$\tau_0$", fontsize=10, color="0.35")
# centred labels in the clear band above each curve's peak
ax.text(0.0, np.exp(5.15) / 2.0 * 1.7, r"$\Delta=5.15$", color=CRIMSON,
        fontsize=10, ha="center")
ax.text(0.0, np.exp(3.8) / 2.0 * 1.55, r"$\Delta=3.8$", color=TEAL,
        fontsize=10, ha="center")
ax.text(0.0, np.exp(2.0) / 2.0 * 1.5, r"$\Delta=2$", color=NAVY,
        fontsize=10, ha="center")
ax.set_xlabel(r"$V$ (V)")
ax.set_ylabel(r"$\tau(V)$ (ns)")
ax.set_title("Relaxation time (fading memory)", fontsize=11.5, pad=10)
fig.savefig(PANELS / "ch02_19_b.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ── (c) sample RTN traces (Delta = 3.8) ──────────────────────────────────
fig, ax = plt.subplots(figsize=(4.6, 3.4))
params = TelegraphParams(tau0=TAU0, Delta=3.8, Vc0=VC0)
DT, NSTEP = 1.0, 600                          # 1 ns steps, 600 ns window
BIAS = [(0.25, CRIMSON), (0.10, TEAL), (0.00, NAVY)]
t = np.arange(NSTEP) * DT
for k, (v, c) in enumerate(BIAS):
    tr = simulate_trace(np.full(NSTEP, v), DT, n=1, params=params, seed=11 + k)[:, 0]
    off = 3.0 * (len(BIAS) - 1 - k)
    ax.step(t, tr + off, where="post", color=c, lw=1.0)
    ax.text(-12, off, rf"$V={v:.2f}$ V", ha="right", va="center", fontsize=10, color=c)
    for lev in (-1, 1):
        ax.axhline(lev + off, color="0.9", lw=0.5, zorder=0)
ax.set_xlim(-95, NSTEP * DT)
ax.set_yticks([])
ax.set_xlabel(r"$t$ (ns)")
ax.set_ylabel(r"$s(t)$  (offset)")
ax.set_title(r"Telegraph traces ($\Delta=3.8$)", fontsize=11.5)
ax.spines["left"].set_visible(False)
fig.savefig(PANELS / "ch02_19_c.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# consistency asserts: panel curves match the module closed forms
assert np.isclose(stationary_mean(0.2, Delta=3.8, Vc0=VC0), np.tanh(3.8 * 0.2 / VC0))
assert np.isclose(relaxation_time(0.0, tau0=1.0, Delta=3.8, Vc0=VC0), np.exp(3.8) / 2.0)
print("wrote 3 panels ->", PANELS)
