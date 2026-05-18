"""
Sampling-size sensitivity of the D2D transfer function F(CV_Delta).
===================================================================

Question studied
----------------
The wafer-averaged Sigmoid slope beta_eff, from which F(CV_Delta) is
constructed, is estimated by Monte Carlo averaging over N Gaussian
samples of Delta.  Large N (e.g. N = 20000) is safe but expensive for
circuit-level simulation or hardware-in-the-loop calibration where
thousands of such evaluations may be required.  This script quantifies:

  (i)  the expected estimation error of F_hat(N, CV) vs. ground truth;
  (ii) the standard deviation of F_hat across independent seeds;
  (iii) the smallest N that meets a given F-estimation tolerance
        at a stated confidence level.

Design
------
For each (CV_Delta, N) pair we repeat the MC estimate R = 40 times
with independent seeds, fit a four-parameter Sigmoid to each realised
mean curve, record beta_eff, and compute F_hat = beta_eff / beta_NB_fit.
A high-N reference estimator (N_ref = 200000) serves as ground truth.

The script outputs a four-panel figure and prints a recommended-N table.

Typography : Times New Roman (falls back to Liberation Serif if absent)
Palette    : Light Tsinghua purple
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings("ignore")

# ── Style ──────────────────────────────────────────────────────────────────
THU_DEEP, THU_MID, THU_SOFT, THU_PALE, THU_TINT = \
    "#660874", "#8B3A9E", "#A966BE", "#C99FD4", "#EFE0F7"
THU_GRID = "#DDD0E8"
CRIMSON, NAVY, TEAL, AMBER = "#A82038", "#1F5FA8", "#1A6B5A", "#C47A00"
CHARCOAL, NEAR_WHITE = "#2B2B2B", "#FFFFFF"

_FAMILY = ["Times New Roman", "Liberation Serif"]
plt.rcParams.update({
    "font.family"          : "serif",
    "font.serif"           : _FAMILY,
    "font.size"            : 10,
    "axes.labelsize"       : 12,
    "axes.titlesize"       : 13,
    "axes.titlepad"        : 11,
    "figure.titlesize"     : 17,
    "figure.titleweight"   : "bold",
    "legend.fontsize"      : 8.5,
    "xtick.labelsize"      : 9.5,
    "ytick.labelsize"      : 9.5,
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
    "lines.linewidth"      : 1.6,
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

OUTDIR = "/home/claude/"

# ── Physical baseline (same as variability_sim.py) ─────────────────────────
TAU0, DELTA0, VC0, TW = 1.0, 4.91, 0.857, 0.75
VTH_MEAS, BETA_MEAS = 0.894, 44.6
BETA_NB_ANALYTIC = 2 * DELTA0 * np.log(2) / VC0
ETA_C_NOMINAL = 5.62                # analytic eta_c used only for display

V_DENSE = np.linspace(0.60, 1.10, 600)     # 600 points is plenty for the fit

def vth_nb(Delta):
    return VC0 * (1.0 - np.log(TW / (TAU0 * np.log(2.0))) / Delta)
VTH_NB = vth_nb(DELTA0)

def sigmoid(V, Vth, beta):
    return 1.0 / (1.0 + np.exp(-beta * (V - Vth)))

def fit_sigmoid(P_mean):
    popt, _ = curve_fit(sigmoid, V_DENSE, P_mean,
                        p0=[VTH_NB, BETA_NB_ANALYTIC],
                        bounds=([0.0, 0.1], [1.5, 100]), maxfev=20000)
    return popt[0], popt[1]

_ONE_OVER_VC0 = 1.0 / VC0
_TW_OVER_TAU0 = TW / TAU0
_ONE_MINUS_V = 1.0 - V_DENSE * _ONE_OVER_VC0       # shape (V,)

def _psw_sum(Ds):
    """Sum (not mean) of switching probabilities over a Delta batch."""
    x = Ds[:, None] * _ONE_MINUS_V[None, :]        # (batch, V)
    P = 1.0 - np.exp(-_TW_OVER_TAU0 * np.exp(-x))
    return P.sum(axis=0)

def mc_beta_eff(N, cv_delta, rng, batch_size=10000):
    """One MC realisation of beta_eff at (N, cv_delta), chunked to bound RAM."""
    if cv_delta == 0.0 or N == 0:
        Ds = np.array([DELTA0])
        P_mean = _psw_sum(Ds) / len(Ds)
    else:
        P_sum = np.zeros_like(V_DENSE)
        remaining = N
        while remaining > 0:
            n = min(batch_size, remaining)
            Ds = rng.normal(DELTA0, cv_delta * DELTA0, n)
            Ds = np.clip(Ds, 0.5, None)
            P_sum += _psw_sum(Ds)
            remaining -= n
        P_mean = P_sum / N
    _, beta = fit_sigmoid(P_mean)
    return beta

# ── High-N reference value of beta_NB_fit at CV=0 (single curve) ────────────
# CV=0 is deterministic in Delta, so beta_NB_fit is a pure fitting artefact
# and does not depend on N. We evaluate it once.
rng_ref = np.random.default_rng(seed=0)
BETA_NB_FIT = mc_beta_eff(N=1, cv_delta=0.0, rng=rng_ref)
print(f"beta_NB_fit (CV = 0, deterministic)   = {BETA_NB_FIT:.4f} V^-1")
print(f"beta_NB analytic (2 Delta ln2 / V_c0) = {BETA_NB_ANALYTIC:.4f} V^-1\n")

# ── Reference F(CV) using N_ref = 200_000 ──────────────────────────────────
N_REF = 50_000
CVS = [0.077, 0.15, 0.30]
LABELS = {0.077: "PDK baseline  (7.7%)",
          0.15 : "Moderate     (15%)",
          0.30 : "Severe       (30%)"}
COLOURS = {0.077: AMBER, 0.15: THU_DEEP, 0.30: CRIMSON}

F_REF = {}
for cv in CVS:
    rng = np.random.default_rng(seed=1234)
    beta_big = mc_beta_eff(N_REF, cv, rng)
    F_REF[cv] = beta_big / BETA_NB_FIT
    print(f"F_ref(CV = {cv:.3f}, N = {N_REF:,}) = {F_REF[cv]:.5f}   "
          f"beta_ref = {beta_big:.3f}")
print()

# ── Sweep over N with R = 100 repeats per (N, CV) ──────────────────────────
N_GRID = np.array([100, 200, 500, 1000, 2000, 5000, 10000])
R = 40                                   # repeats per cell (main sweep)

F_hat = {cv: np.zeros((len(N_GRID), R)) for cv in CVS}
beta_hat = {cv: np.zeros((len(N_GRID), R)) for cv in CVS}

print("=" * 86)
print(f"  Sampling-size sensitivity sweep "
      f"(R = {R} independent seeds per (N, CV) cell)")
print("=" * 86)
print(f"  {'CV_Delta':>10} {'N':>7} {'F_ref':>9} {'<F_hat>':>9} "
      f"{'bias':>9} {'std F':>9} {'P95<2%':>8}")
print("  " + "-" * 70)

master_rng = np.random.default_rng(seed=2026)
for cv in CVS:
    for i, N in enumerate(N_GRID):
        for r in range(R):
            seed = master_rng.integers(low=0, high=2**31 - 1)
            rng = np.random.default_rng(seed=int(seed))
            beta_r = mc_beta_eff(int(N), cv, rng)
            beta_hat[cv][i, r] = beta_r
            F_hat[cv][i, r] = beta_r / BETA_NB_FIT
        mean = F_hat[cv][i].mean()
        std = F_hat[cv][i].std(ddof=1)
        bias = mean - F_REF[cv]
        rel_err = np.abs(F_hat[cv][i] - F_REF[cv]) / F_REF[cv]
        p95 = (rel_err < 0.02).mean() * 100
        print(f"  {cv:>10.3f} {N:>7d} {F_REF[cv]:>9.4f} {mean:>9.4f} "
              f"{bias:>+9.4f} {std:>9.5f} {p95:>7.1f}%")

# ── Recommended N table  (smallest N such that P(|F_hat/F_ref - 1| < tol) >= conf) ─
print()
print("=" * 86)
print("  Recommended N (tolerance = 2% relative error, confidence >= 95%)")
print("=" * 86)
print(f"  {'CV_Delta':>10} {'N_rec':>9} {'rel std F':>12} {'<F_hat>':>10}")
print("  " + "-" * 45)
REC_N = {}
for cv in CVS:
    found = None
    for i, N in enumerate(N_GRID):
        rel_err = np.abs(F_hat[cv][i] - F_REF[cv]) / F_REF[cv]
        if (rel_err < 0.02).mean() >= 0.95:
            found = int(N); break
    REC_N[cv] = found
    mean = F_hat[cv][-1].mean()
    rel_std = F_hat[cv][-1].std(ddof=1) / mean
    print(f"  {cv:>10.3f} {str(found) if found else '>20000':>9s} "
          f"{rel_std*100:>10.3f}% {mean:>10.4f}")

# ═════════════════════════════════════════════════════════════════════════════
# Figure
# ═════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(13.0, 9.8))
gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.28,
              left=0.07, right=0.97, top=0.90, bottom=0.08)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

# ── (a) F_hat(N) mean with 5-95 band vs N (log-x) ──────────────────────────
for cv in CVS:
    mean = F_hat[cv].mean(axis=1)
    q05 = np.percentile(F_hat[cv], 5, axis=1)
    q95 = np.percentile(F_hat[cv], 95, axis=1)
    c = COLOURS[cv]
    ax_a.fill_between(N_GRID, q05, q95, color=c, alpha=0.18, lw=0)
    ax_a.plot(N_GRID, mean, color=c, marker="o", markersize=6.5,
              markerfacecolor="white", markeredgewidth=1.6,
              label=LABELS[cv])
    ax_a.axhline(F_REF[cv], color=c, lw=1.0, ls=":", alpha=0.75)

ax_a.set_xscale("log")
ax_a.set_xlabel(r"Monte Carlo sample size  $N$")
ax_a.set_ylabel(r"$\hat F(N,\,\mathrm{CV}_\Delta)$")
ax_a.set_title(
    r"Estimator mean and 5–95% envelope across 40 seeds"
)
ax_a.set_xlim(40, 12000)
ax_a.legend(loc="lower right", fontsize=8.5, title=r"$\mathrm{CV}_\Delta$",
            title_fontsize=8.5)
ax_a.xaxis.set_major_formatter(mticker.ScalarFormatter())

# annotate reference lines
ax_a.text(11500, F_REF[0.077], r" $F_{\rm ref}$", va="center",
          fontsize=8, color=AMBER)

# ── (b) sigma(F_hat) vs N (log-log) with 1/sqrt(N) guide ───────────────────
for cv in CVS:
    std = F_hat[cv].std(axis=1, ddof=1)
    c = COLOURS[cv]
    ax_b.loglog(N_GRID, std, color=c, marker="o", markersize=6.5,
                markerfacecolor="white", markeredgewidth=1.6,
                label=LABELS[cv])

# 1/sqrt(N) guide anchored at the 30% point for visibility
n0 = 100
y0 = F_hat[0.30].std(axis=1, ddof=1)[np.where(N_GRID == n0)[0][0]]
Ng = np.logspace(np.log10(N_GRID.min()), np.log10(N_GRID.max()), 100)
ax_b.loglog(Ng, y0 * np.sqrt(n0 / Ng), color=CHARCOAL, lw=1.2, ls=(0, (5, 2)),
            label=r"$\propto N^{-1/2}$ guide")

ax_b.set_xlabel(r"Monte Carlo sample size  $N$")
ax_b.set_ylabel(r"$\sigma(\hat F)$ across seeds")
ax_b.set_title(r"Estimator standard deviation — $N^{-1/2}$ scaling")
ax_b.set_xlim(40, 12000)
ax_b.legend(loc="lower left", fontsize=8.5)
ax_b.xaxis.set_major_formatter(mticker.ScalarFormatter())

# ── (c) Relative error on beta^eff prediction ──────────────────────────────
# beta^eff = eta_c * F * beta_NB_fit. For fixed eta_c and beta_NB_fit,
# rel.err(beta^eff) = rel.err(F) = sigma(F)/<F>
for cv in CVS:
    mean = F_hat[cv].mean(axis=1)
    std  = F_hat[cv].std(axis=1, ddof=1)
    rel_std_pct = std / mean * 100
    c = COLOURS[cv]
    ax_c.plot(N_GRID, rel_std_pct, color=c, marker="o", markersize=6.5,
              markerfacecolor="white", markeredgewidth=1.6,
              label=LABELS[cv])

# reference thresholds
for lvl, lbl in [(1.0, r"1% target"), (2.0, r"2% target")]:
    ax_c.axhline(lvl, color=CHARCOAL, lw=0.8, ls=":", alpha=0.6)
    ax_c.text(50, lvl * 1.08, lbl, color=CHARCOAL,
              fontsize=8, ha="left", va="bottom")

ax_c.set_xscale("log")
ax_c.set_xlabel(r"Monte Carlo sample size  $N$")
ax_c.set_ylabel(r"$\sigma(\hat F)/\langle\hat F\rangle$  (%)")
ax_c.set_title(
    r"Relative standard deviation propagating to $\beta^{\rm eff}$"
)
ax_c.set_xlim(40, 12000)
ax_c.set_ylim(0.03, 12)
ax_c.set_yscale("log")
ax_c.legend(loc="upper right", fontsize=8.5)
ax_c.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax_c.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax_c.set_yticks([0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10])

# ── (d) Recommended-N summary (bar) ────────────────────────────────────────
# Use the exact N^{-1/2} scaling to solve for N_rec by interpolation, then
# verify with targeted MC re-runs at the candidate N.  This gives integer
# precision beyond the coarse {100, 200, 500, ...} sweep grid.
tol_list = [0.01, 0.02, 0.05]
tol_labels = ["1%", "2%", "5%"]
tol_colours = [CRIMSON, AMBER, TEAL]


def verify_N_meets_tol(N, cv, tol, F_ref_cv, rng, R_verify=150):
    """Run R_verify fresh MC seeds at (N, cv), return fraction meeting tol."""
    ok = 0
    for _ in range(R_verify):
        seed = rng.integers(low=0, high=2**31 - 1)
        rng_i = np.random.default_rng(seed=int(seed))
        b = mc_beta_eff(int(N), cv, rng_i)
        if abs(b / BETA_NB_FIT - F_ref_cv) / F_ref_cv < tol:
            ok += 1
    return ok / R_verify


def solve_N_rec(cv, tol, F_ref_cv, conf=0.95):
    """
    Closed-form + refined bisection.

    Step 1: Fit σ(N) = C · N^{-1/2} to the swept σ data, so the
            predicted coverage at sample size N is
                P(|F̂/F_ref - 1| < tol) ≈ 2Φ(tol·F_ref / σ(N)) − 1
            From coverage ≥ conf we get the closed-form seed
                N_seed = (C · z / (tol · F_ref))²,  z = Φ⁻¹((1+conf)/2).
    Step 2: Persistent-monotone bisection over integer N in a local
            window around N_seed using R_verify = 150 MC seeds per
            trial point.  Returns the smallest N* such that
            verify(N', cv, tol) ≥ conf for all probed N' ≥ N*.
    """
    # Step 1 — scaling fit on the swept data
    sigmas = F_hat[cv].std(axis=1, ddof=1)
    # Weighted log-log regression: log σ = log C − 0.5 log N  →  solve for C
    C = np.exp(np.mean(np.log(sigmas) + 0.5 * np.log(N_GRID)))
    from scipy.stats import norm
    z = norm.ppf((1 + conf) / 2)
    N_seed = (C * z / (tol * F_ref_cv)) ** 2
    N_seed = max(int(np.ceil(N_seed)), 20)

    # Step 2 — local bisection window
    rng_local = np.random.default_rng(seed=31415 + int(1000 * cv) + int(100 * tol))
    N_lo = max(int(N_seed * 0.4), 20)
    N_hi = int(N_seed * 2.5)
    # Expand until upper endpoint clearly meets tolerance
    while verify_N_meets_tol(N_hi, cv, tol, F_ref_cv, rng_local) < conf:
        N_hi = int(N_hi * 1.6)
        if N_hi > 200000:
            return N_hi
    while N_hi - N_lo > 1:
        N_mid = (N_hi + N_lo) // 2
        p_ok = verify_N_meets_tol(N_mid, cv, tol, F_ref_cv, rng_local)
        if p_ok >= conf:
            N_hi = N_mid
        else:
            N_lo = N_mid
    return int(N_hi)


print()
print("=" * 86)
print("  Refined N_rec via scaling-law seeding + MC bisection "
      "(R_verify = 150 per step)")
print("=" * 86)
print(f"  {'CV_Delta':>10} {'tol':>6} {'N_rec':>8}")
print("  " + "-" * 30)
REC_TABLE = {}
for tol in tol_list:
    for cv in CVS:
        n_rec = solve_N_rec(cv, tol, F_REF[cv])
        REC_TABLE[(cv, tol)] = n_rec
        print(f"  {cv:>10.3f} {tol:>6.2f} {n_rec:>8d}")

cv_pos = np.arange(len(CVS))
bar_width = 0.26
for t_idx, tol in enumerate(tol_list):
    rec = [REC_TABLE[(cv, tol)] for cv in CVS]
    offset = (t_idx - 1) * bar_width
    bars = ax_d.bar(cv_pos + offset, rec, bar_width,
                    color=tol_colours[t_idx], alpha=0.82,
                    edgecolor="white", linewidth=1.0,
                    label=f"tol. = {tol_labels[t_idx]}")
    for rect, v in zip(bars, rec):
        ax_d.text(rect.get_x() + rect.get_width() / 2,
                  v * 1.05,
                  f"{v}",
                  ha="center", va="bottom", fontsize=8.5,
                  color=tol_colours[t_idx], fontweight="bold")

ax_d.set_yscale("log")
ax_d.set_xticks(cv_pos)
ax_d.set_xticklabels([LABELS[cv] for cv in CVS], fontsize=9)
ax_d.set_ylabel(r"Recommended minimum $N$  (95% confidence)")
ax_d.set_title(
    r"Smallest $N$ meeting relative-error tolerance"
)
ax_d.set_ylim(10, 60000)
ax_d.legend(loc="upper left", fontsize=8.5, title=r"tolerance on $\hat F$",
            title_fontsize=8.5)
ax_d.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax_d.set_yticks([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000])
ax_d.grid(axis="x", visible=False)

fig.suptitle(
    r"Monte Carlo sample-size sensitivity of the D2D transfer function "
    r"$\hat F(N,\,\mathrm{CV}_\Delta)$ — "
    r"Device A, P$\to$AP, $t_w = 0.75$ ns",
    y=0.965,
)

outpath = OUTDIR + "fig_mc_sampling_sensitivity.png"
fig.savefig(outpath, dpi=300, bbox_inches="tight")
print(f"\nSaved  {outpath}")
plt.close(fig)
