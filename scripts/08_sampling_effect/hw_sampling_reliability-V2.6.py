"""
Hardware Bernoulli sampling reliability for sMTJ-based stochastic units.
=======================================================================

Compares two paths to the smallest sampling count K_req(p, ε, α):

  (1) Binomial-CDF exact solver
        — covers every K ≥ 1 with no sampling noise, uses persistent-
          monotone bisection to handle the discrete coverage staircase;
  (2) Monte Carlo empirical solver
        — M independent p̂_K replicates per K, coverage estimated from
          the fraction within the ε-band; models what a hardware-in-the-
          loop calibration actually observes.

Figure has 6 panels (2 × 3):
  (a) Binomial PMF of p̂_K at p=0.5 for K ∈ {5,20,100,500}
  (b) σ(p̂_K) exact vs. MC markers vs. K^{-1/2}
  (c) Coverage at p=0.5: exact staircase + MC band, three ε curves
  (d) MC coverage estimator RMSE vs. replicate count M
  (e) K_req: exact vs. MC vs. CLT at three (p, ε) points
  (f) K_req(ε) at p=0.5 — exact points, CLT line, ε^{-2} scaling
"""

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec

# ── Style ────────────────────────────────────────────────────────────────
THU_DEEP, THU_MID, THU_SOFT, THU_PALE, THU_TINT = \
    "#660874", "#8B3A9E", "#A966BE", "#C99FD4", "#EFE0F7"
THU_GRID = "#DDD0E8"
CRIMSON, NAVY, TEAL, AMBER = "#A82038", "#1F5FA8", "#1A6B5A", "#C47A00"
CHARCOAL, NEAR_WHITE = "#2B2B2B", "#FFFFFF"

_FAMILY = ["Arial", "Liberation Sans"]
plt.rcParams.update({
    "font.family"          : "sans-serif",
    "font.sans-serif"      : _FAMILY,
    "font.size"            : 13,
    "axes.labelsize"       : 14,
    "axes.titlesize"       : 15,
    "axes.titlepad"        : 9,
    "figure.titlesize"     : 18,
    "figure.titleweight"   : "bold",
    "legend.fontsize"      : 12,
    "xtick.labelsize"      : 13,
    "ytick.labelsize"      : 13,
    "mathtext.fontset"     : "stix",
    "axes.linewidth"       : 0.9,
    "axes.edgecolor"       : CHARCOAL,
    "axes.facecolor"       : NEAR_WHITE,
    "axes.labelcolor"      : CHARCOAL,
    "axes.spines.top"      : False,
    "axes.spines.right"    : False,
    "axes.grid"            : True,
    "grid.color"           : THU_GRID,
    "grid.linewidth"       : 0.55,
    "grid.linestyle"       : "--",
    "grid.alpha"           : 0.7,
    "xtick.direction"      : "in",
    "ytick.direction"      : "in",
    "xtick.major.size"     : 4,
    "ytick.major.size"     : 4,
    "xtick.minor.visible"  : True,
    "ytick.minor.visible"  : True,
    "lines.linewidth"      : 1.5,
    "legend.frameon"       : True,
    "legend.framealpha"    : 0.92,
    "legend.edgecolor"     : THU_PALE,
    "legend.facecolor"     : NEAR_WHITE,
    "legend.handlelength"  : 1.8,
    "legend.handletextpad" : 0.5,
    "legend.labelspacing"  : 0.35,
    "figure.dpi"           : 150,
    "savefig.dpi"          : 300,
    "savefig.bbox"         : "tight",
    "figure.facecolor"     : NEAR_WHITE,
})

from pathlib import Path
OUTDIR = str(Path(__file__).resolve().parent) + "/"
RNG = np.random.default_rng(seed=20260420)

# Binomial coverage / K_req — exact + MC + CLT — from vgsot-sim analysis.
from vgsot_sim.analysis.sampling import (
    exact_coverage, K_required_persistent as find_K_req,
    mc_coverage, K_required_persistent_mc as find_K_req_mc,
)
def coverage_array(K_arr, p, eps):
    return np.array([exact_coverage(int(K), p, eps) for K in K_arr])

# ═════════════════════════════════════════════════════════════════════════
# 4. Working-point voltages
# ═════════════════════════════════════════════════════════════════════════
V_TH, BETA = 0.894, 44.6
def V_for_p(p):
    return V_TH + np.log(p / (1 - p)) / BETA

print("=" * 78)
print("  Operating-point voltages on the primary Sigmoid")
print("  (Device A, P→AP, 0.75 ns)")
print("=" * 78)
for p in [0.1, 0.5, 0.9]:
    print(f"  p = {p:.1f}   →   V = {V_for_p(p)*1e3:7.2f} mV")

# ═════════════════════════════════════════════════════════════════════════
# 5. Exact vs. MC vs. CLT  K_req  reconciliation
# ═════════════════════════════════════════════════════════════════════════
PS  = [0.1, 0.5, 0.9]
EPS = [0.02, 0.05, 0.10]
z95 = 1.959964
M_MC_BIG = 20_000

print()
print("=" * 95)
print("  K_req(p, ε)  —  exact (Binomial) vs. MC (persistence) vs. CLT")
print("=" * 95)
print(f"  {'p':>4} {'ε':>6}  {'K_req(ex)':>11}  {'K_req(MC)':>11}  "
      f"{'K_CLT':>8}  {'ΔMC':>8}  {'ΔCLT':>8}")
print("  " + "-" * 68)
K_REQ_EXACT, K_REQ_MC = {}, {}
for p in PS:
    for eps in EPS:
        K_ex = find_K_req(p, eps)
        K_mc = find_K_req_mc(p, eps, M=M_MC_BIG, rng=RNG)
        K_cl = z95**2 * p * (1 - p) / eps**2
        K_REQ_EXACT[(p, eps)] = K_ex
        K_REQ_MC[(p, eps)]    = K_mc
        dmc = (K_mc - K_ex) if K_mc else 0
        print(f"  {p:>4.1f} {eps:>6.2f}  {K_ex:>11d}  "
              f"{(str(K_mc) if K_mc else 'fail'):>11s}  "
              f"{K_cl:>8.1f}  {dmc:>+8d}  {K_cl - K_ex:>+8.1f}")

# ═════════════════════════════════════════════════════════════════════════
# 6. Figure  (2 × 3)
# ═════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(14.0, 9.0))
gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.30,
              left=0.06, right=0.98, top=0.90, bottom=0.08)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[0, 2])
ax_d = fig.add_subplot(gs[1, 0])
ax_e = fig.add_subplot(gs[1, 1])
ax_f = fig.add_subplot(gs[1, 2])

# ── (a) PMF of p̂_K at p=0.5 ─────────────────────────────────────────────
K_HIST = [5, 20, 100, 500]
COL_K  = [CRIMSON, AMBER, THU_DEEP, NAVY]
p_demo, eps_demo = 0.5, 0.05

ax_a.axvspan(p_demo - eps_demo, p_demo + eps_demo,
             color=AMBER, alpha=0.14, zorder=0,
             label=rf"$|{{\hat p}}-p|<\varepsilon={eps_demo:.2f}$")
ax_a.axvline(p_demo, color=CHARCOAL, lw=0.8, ls=":", alpha=0.55)

for K, col in zip(K_HIST, COL_K):
    k = np.arange(K + 1)
    pmf = stats.binom.pmf(k, K, p_demo)
    phat = k / K
    cov = exact_coverage(K, p_demo, eps_demo)
    mask = pmf >= 1e-4
    ax_a.plot(phat[mask], pmf[mask], "o-", color=col, lw=1.0,
              markersize=4.5, markerfacecolor="white", markeredgewidth=1.2,
              label=rf"$K={K}$,  cov.$={cov*100:.1f}\%$")

ax_a.set_xlabel(r"Empirical frequency  $\hat p_K = k/K$")
ax_a.set_ylabel(r"Probability mass  $P(\hat p_K = k/K)$")
ax_a.set_title(r"Binomial PMF of $\hat p_K$ at $p = 0.5$")
ax_a.set_xlim(-0.02, 1.02)
ax_a.set_yscale("log")
ax_a.set_ylim(8e-5, 0.5)
ax_a.legend(loc="lower center", ncol=2)

# ── (b) σ(p̂_K) — exact line + MC markers ────────────────────────────────
K_sweep = np.logspace(np.log10(2), np.log10(5000), 500)
COL_P = {0.1: NAVY, 0.5: CRIMSON, 0.9: TEAL}

for p in PS:
    sig = np.sqrt(p * (1 - p) / K_sweep)
    ax_b.loglog(K_sweep, sig, color=COL_P[p], lw=1.6,
                label=rf"$p = {p:.1f}$ (exact)")

K_mc_std = np.array([3, 8, 20, 50, 150, 500, 1500, 4500])
M_std = 3000
for p in PS:
    sig_mc = []
    for K in K_mc_std:
        phat = RNG.binomial(K, p, size=M_std) / K
        sig_mc.append(phat.std(ddof=1))
    ax_b.plot(K_mc_std, sig_mc, "s", color=COL_P[p],
              markersize=6.0, markerfacecolor="white", markeredgewidth=1.3,
              zorder=6)

ax_b.plot([], [], "s", color=CHARCOAL, markersize=6.0,
          markerfacecolor="white", markeredgewidth=1.3,
          label=rf"MC ($M={M_std}$)")

ax_b.axvline(100, color=CHARCOAL, lw=0.8, ls=":", alpha=0.55)

ax_b.set_xlabel(r"Bernoulli samples  $K$")
ax_b.set_ylabel(r"Std. error  $\sigma(\hat p_K)$")
ax_b.set_title(r"Std. error: exact vs. MC vs. $K^{-1/2}$ scaling")
ax_b.set_xlim(2, 5000)
ax_b.set_ylim(0.003, 0.4)
ax_b.legend(loc="upper right")
ax_b.xaxis.set_major_formatter(mticker.ScalarFormatter())

# ── (c) Coverage at p=0.5 — exact staircase + MC band ────────────────────
Ks_c = np.unique(np.concatenate([
    np.arange(2, 120),
    np.round(np.logspace(2.1, np.log10(3000), 200)).astype(int)
]))
COL_E = [NAVY, THU_DEEP, CRIMSON]
M_HITL = 500

for eps, col in zip(EPS, COL_E):
    cov_exact = coverage_array(Ks_c, 0.5, eps)
    ax_c.plot(Ks_c, cov_exact, color=col, lw=1.5, alpha=0.95,
              label=rf"$\varepsilon = {eps:.2f}$ (exact)")
    cov_mc = np.array([mc_coverage(K, 0.5, eps, M_HITL, RNG) for K in Ks_c])
    band = 1.96 * np.sqrt(cov_mc * (1 - cov_mc) / M_HITL)
    ax_c.fill_between(Ks_c, cov_mc - band, cov_mc + band,
                       color=col, alpha=0.15, lw=0)

    K_req = K_REQ_EXACT[(0.5, eps)]
    cov_req = exact_coverage(K_req, 0.5, eps)
    ax_c.scatter([K_req], [cov_req], s=140, marker="*",
                 facecolor=col, edgecolor="white", linewidths=1.4, zorder=9)
    ax_c.annotate(rf"$K_{{\rm req}}={K_req}$",
                  xy=(K_req, cov_req),
                  xytext=(K_req * 1.45, cov_req - 0.12),
                  color=col, fontweight="bold",
                  arrowprops=dict(arrowstyle="-", color=col, lw=0.7))

ax_c.axhline(0.95, color=CHARCOAL, lw=1.0, ls=":", alpha=0.8)
ax_c.text(2.3, 0.952, r"95% target", color=CHARCOAL, ha="left", va="bottom")
ax_c.plot([], [], color=CHARCOAL, lw=5, alpha=0.15,
          label=rf"MC 95% band ($M={M_HITL}$)")

ax_c.set_xscale("log")
ax_c.set_xlabel(r"Bernoulli samples  $K$")
ax_c.set_ylabel(r"Coverage  $P(|\hat p_K - p| < \varepsilon)$")
ax_c.set_title(r"Coverage at $p=0.5$: exact staircase vs. MC band")
ax_c.set_xlim(2, 3000)
ax_c.set_ylim(0, 1.03)
ax_c.legend(loc="lower right")
ax_c.xaxis.set_major_formatter(mticker.ScalarFormatter())

# ── (d) MC coverage RMSE vs M ────────────────────────────────────────────
K_fix, p_fix, eps_fix = 100, 0.5, 0.10
cov_true = exact_coverage(K_fix, p_fix, eps_fix)

M_sweep = np.logspace(np.log10(10), np.log10(30000), 18).astype(int)
R_rep   = 200
rms_err, mean_err = [], []
for M in M_sweep:
    errs = np.array([mc_coverage(K_fix, p_fix, eps_fix, M, RNG) - cov_true
                     for _ in range(R_rep)])
    rms_err.append(np.sqrt(np.mean(errs**2)))
    mean_err.append(np.mean(errs))
rms_err = np.array(rms_err)

rms_theory = np.sqrt(cov_true * (1 - cov_true) / M_sweep)

ax_d.loglog(M_sweep, rms_err, "o", color=THU_DEEP, markersize=6.5,
            markerfacecolor="white", markeredgewidth=1.6,
            label=rf"MC RMSE ({R_rep} seeds)", zorder=5)
ax_d.loglog(M_sweep, rms_theory, color=CHARCOAL, lw=1.3, ls=(0, (5, 2)),
            label=r"$\sqrt{\mathrm{cov}(1-\mathrm{cov})/M}$")
ax_d.axhline(0.01, color=AMBER, lw=1.0, ls=":", alpha=0.8)
ax_d.text(12, 0.0108, r"1% RMSE", color=AMBER, ha="left", va="bottom")

ax_d.set_xlabel(r"MC replicate count  $M$")
ax_d.set_ylabel(r"RMSE of $\hat{\mathrm{cov}}_{\rm MC}$")
ax_d.set_title(r"MC estimator noise: RMSE vs. $M$")
ax_d.set_xlim(8, 40000)
ax_d.set_ylim(1.5e-3, 0.3)
ax_d.legend(loc="upper right")
ax_d.xaxis.set_major_formatter(mticker.ScalarFormatter())

# ── (e) K_req bars — exact vs. MC vs. CLT ────────────────────────────────
x_pos = np.arange(len(EPS))
bar_w = 0.23
offsets = {0.1: -bar_w, 0.5: 0, 0.9: bar_w}

for p in PS:
    Ks_ex = [K_REQ_EXACT[(p, e)] for e in EPS]
    Ks_mc = [K_REQ_MC[(p, e)]    for e in EPS]
    col = COL_P[p]
    bars = ax_e.bar(x_pos + offsets[p], Ks_ex, bar_w,
                    color=col, alpha=0.85, edgecolor="white", lw=1.0,
                    label=rf"$p = {p:.1f}$ (exact)")
    for rect, v in zip(bars, Ks_ex):
        ax_e.text(rect.get_x() + rect.get_width() / 2, v * 1.12,
                  f"{v}", ha="center", va="bottom",
                  color=col, fontweight="bold")
    ax_e.plot(x_pos + offsets[p], Ks_mc, "s", color=col,
              markersize=5.2, markerfacecolor="white", markeredgewidth=1.3,
              zorder=7)
    Kcl = [z95**2 * p * (1 - p) / e**2 for e in EPS]
    ax_e.plot(x_pos + offsets[p], Kcl, "D", color=col,
              markersize=4.8, markerfacecolor=NEAR_WHITE,
              markeredgewidth=1.2, alpha=0.9, zorder=6)

ax_e.plot([], [], "s", color=CHARCOAL, markersize=5.2,
          markerfacecolor="white", markeredgewidth=1.3,
          label=rf"MC ($M={M_MC_BIG//1000}$k)")
ax_e.plot([], [], "D", color=CHARCOAL, markersize=4.8,
          markerfacecolor=NEAR_WHITE, markeredgewidth=1.2,
          label=r"CLT  $z^2 p(1-p)/\varepsilon^2$")

ax_e.set_yscale("log")
ax_e.set_xticks(x_pos)
ax_e.set_xticklabels([rf"$\varepsilon = {e:.2f}$" for e in EPS])
ax_e.set_ylabel(r"$K_{\rm req}$ at 95% confidence")
ax_e.set_title(r"$K_{\rm req}$: exact vs. MC vs. CLT")
ax_e.set_ylim(8, 15000)
ax_e.legend(loc="upper right", ncol=1)
ax_e.grid(axis="x", visible=False)

# ── (f) K_req(ε) at p=0.5 ────────────────────────────────────────────────
EPS_FINE = np.array([0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05,
                     0.07, 0.10, 0.15, 0.20])
Ks_fine = [find_K_req(0.5, e) for e in EPS_FINE]
K_cl_fine = z95**2 * 0.25 / EPS_FINE**2

ax_f.loglog(EPS_FINE * 100, Ks_fine, "o", color=CRIMSON,
            markersize=6.5, markerfacecolor="white", markeredgewidth=1.6,
            label=r"Exact $K_{\rm req}$ (Binomial)", zorder=5)
ax_f.loglog(EPS_FINE * 100, K_cl_fine, color=CHARCOAL, lw=1.3,
            ls=(0, (5, 2)),
            label=r"CLT  $\approx z^2/(4\varepsilon^2)$")

e0 = 0.05
K0 = find_K_req(0.5, e0)
eps_guide = np.logspace(np.log10(0.01), np.log10(0.20), 80)
K_guide = K0 * (e0 / eps_guide)**2
ax_f.loglog(eps_guide * 100, K_guide, color=THU_PALE, lw=1.0, ls=":",
            label=r"$\propto \varepsilon^{-2}$ guide", zorder=3)

ax_f.axhline(100, color=CHARCOAL, lw=0.8, ls=":", alpha=0.55)
ax_f.text(1.0, 110, r"$K = 100$",
          color=CHARCOAL, ha="left", va="bottom")

ax_f.set_xlabel(r"Target precision  $\varepsilon$ (%)")
ax_f.set_ylabel(r"$K_{\rm req}$ at $p = 0.5$")
ax_f.set_title(r"$K_{\rm req}(\varepsilon)$ at worst-case $p=0.5$")
ax_f.set_xlim(0.8, 25)
ax_f.set_ylim(20, 40000)
ax_f.legend(loc="upper right")
ax_f.xaxis.set_major_formatter(mticker.ScalarFormatter())

fig.suptitle(
    r"Hardware Bernoulli sampling reliability — Binomial-exact vs. Monte Carlo "
    r"(Device A, P$\to$AP, $t_w = 0.75$ ns)",
    y=0.965,
)

outpath = OUTDIR + "fig_hw_sampling_reliability.png"
fig.savefig(outpath, dpi=300, bbox_inches="tight")
print(f"\nSaved  {outpath}")
plt.close(fig)

# ═════════════════════════════════════════════════════════════════════════
# 7. §2.3.3 K = 100 context
# ═════════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("  Context: the §2.3.3 K = 100 characterisation protocol")
print("=" * 78)
K_REF = 100
for p in PS:
    lo, hi = 1e-4, 0.5
    for _ in range(60):
        mid = (lo + hi) / 2
        if exact_coverage(K_REF, p, mid) >= 0.95:
            hi = mid
        else:
            lo = mid
    eps_95 = hi
    print(f"  p = {p:.1f}:  at K = 100, 95% CI half-width  ε ≈ "
          f"{eps_95:.4f}  ({eps_95*100:.2f}%)")
