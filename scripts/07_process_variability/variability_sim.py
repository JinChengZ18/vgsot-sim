"""
sMTJ variability Monte Carlo simulation  (v5)
=============================================
Quantifies the effect of device-to-device (D2D) variability in the thermal
stability factor Delta on the wafer-averaged switching probability curve,
with the variability budget now grounded in the PDK local-mismatch
parameters (Rp_mis, Rsot_mis, TMR_mis) via a Brinkman-model decomposition.

Changelog (v5)
--------------
1. Panel top-right now shows only C2C-calibrated curves (eta_c applied).
   Raw NB curves are no longer overlaid with the measured Sigmoid because
   doing so invited the wrong visual comparison — NB over-broadens C2C by
   construction, and the resolution is the C2C empirical factor, not a
   direct match between raw NB and measurement.
2. Title-hierarchy strengthened: suptitle (17 pt bold) > panel title
   (13 pt) > axis label (12 pt), with increased pad between levels.
3. Panel top-left replaced with the PDK-derived CV(Delta) budget stack:
   shows how CV(R_P) through Brinkman splits into CV(RA) and CV(D), then
   propagates through (H_k, M_s, V_mag) to yield the baseline CV(Delta).
4. PDK baseline (CV_Delta^PDK ~ 7.7 %) is now marked consistently on
   panels (c) and (d) as an amber star, replacing the previous arbitrary
   10 % reference point.

Typography : TeX Gyre Termes, math unified
Palette    : Light Tsinghua purple
"""

# ─────────────────────────────────────────────────────────────────────────────
# 0. Imports
# ─────────────────────────────────────────────────────────────────────────────
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings("ignore")

# Brinkman variance budget + D2D wafer-averaged P_sw + sigmoid slope fit
# all delegate to vgsot-sim's analysis sub-package.
from vgsot_sim.analysis import nb_fit
from vgsot_sim.analysis.variability import (
    PDKBudgetInputs, cv_delta_budget, wafer_average_psw, fit_sigmoid_to_average,
)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Style
# ─────────────────────────────────────────────────────────────────────────────
THU_DEEP   = "#660874"
THU_MID    = "#8B3A9E"
THU_SOFT   = "#A966BE"
THU_PALE   = "#C99FD4"
THU_TINT   = "#EFE0F7"
THU_GRID   = "#DDD0E8"
CRIMSON    = "#A82038"
NAVY       = "#1F5FA8"
TEAL       = "#1A6B5A"
AMBER      = "#C47A00"
CHARCOAL   = "#2B2B2B"
NEAR_WHITE = "#FFFFFF"

THU_LIGHT_CMAP = LinearSegmentedColormap.from_list("thu_light", [
    (0.00, "#FFFFFF"),
    (0.20, "#F5E8FA"),
    (0.40, "#DCBEEA"),
    (0.60, "#B47ACA"),
    (0.80, "#8A3EA8"),
    (1.00, "#660874"),
])

_FAMILY = ["Arial", "Liberation Sans"]
plt.rcParams.update({
    "font.family"          : "sans-serif",
    "font.sans-serif"      : _FAMILY,
    "font.size"            : 13,
    "axes.labelsize"       : 14,
    "axes.titlesize"       : 15,
    "axes.titlepad"        : 9,
    "figure.titlesize"     : 20,
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

# ─────────────────────────────────────────────────────────────────────────────
# 2. Baseline device parameters  (SAME-BATCH: Device A, AP→P, 0.75 ns)
# ─────────────────────────────────────────────────────────────────────────────
# At t_w = 0.75 ns the Psw(V) measurement was performed on the same device
# (Device A) and in the same batch as the hysteresis measurements that
# yielded Delta and V_c0.  The AP→P direction is chosen as primary reference
# because the P→AP data exhibits a back-hopping plateau near V ~ 940-1000 mV
# that corrupts a single-Sigmoid fit; the AP→P transition is clean
# (R^2 = 0.993).  A separate 5 ns Sigmoid batch (beta_s = 56.9 V^-1) is kept
# as an *external* reference, NOT as the baseline — see §2.3.5 of the
# main document.
TAU0      = 1.0           # ns (literature prior)
DELTA0    = 4.91          # AP→P, from same-batch hysteresis
VC0       = 0.857         # V
TW        = 0.75          # ns (same-batch operating point)
VTH_MEAS  = 0.894         # V, same-batch Sigmoid fit
BETA_MEAS = 44.6          # V^-1 (= 1/k with k = 22.43 mV), same-batch
BATCH_TAG = "Device A, AP→P, 0.75 ns (same batch)"
BETA_NB_ANALYTIC = nb_fit.beta_nb_analytic(DELTA0, VC0)  # ~7.94 V^-1

# ─────────────────────────────────────────────────────────────────────────────
# 3. PDK-based CV(Delta) budget  (Brinkman decomposition)
# ─────────────────────────────────────────────────────────────────────────────
# Inputs match the prior inline values (Hikstor PDK Rp_mis / Rsot_mis / TMR_mis
# + 300 mm process control on t_ox / phi / t_f / M_s). The Brinkman
# decomposition itself now lives in vgsot_sim.analysis.variability.
PDK_INPUTS = PDKBudgetInputs(
    CV_RP=0.07, CV_RSOT=0.07, CV_TMR=0.04,
    CV_TOX=0.003, CV_PHI=0.005, CV_TF=0.003, CV_MS=0.02,
    T_OX=1.0, PHI_BAR=0.6, KAPPA=10.25,
)
_budget = cv_delta_budget(PDK_INPUTS)
CV_RP   = PDK_INPUTS.CV_RP
CV_RSOT = PDK_INPUTS.CV_RSOT
CV_TMR  = PDK_INPUTS.CV_TMR
KAPPA   = PDK_INPUTS.KAPPA
T_OX    = PDK_INPUTS.T_OX
PHI_BAR = PDK_INPUTS.PHI_BAR
CV_TOX  = PDK_INPUTS.CV_TOX
CV_PHI  = PDK_INPUTS.CV_PHI
CV_TF   = PDK_INPUTS.CV_TF
CV_MS   = PDK_INPUTS.CV_MS
SENS_TOX = KAPPA * np.sqrt(PHI_BAR)
SENS_PHI = KAPPA * T_OX / (2 * np.sqrt(PHI_BAR))
CV_RA   = _budget.CV_RA
CV_D    = _budget.CV_D
CV_V    = _budget.CV_V
CV_HK   = _budget.CV_HK
CV_DELTA_PDK    = _budget.CV_Delta
CV_DELTA_PDK_SQ = CV_DELTA_PDK ** 2

print("=" * 72)
print("  PDK-based CV(Delta) budget  (Brinkman decomposition of CV(R_P))")
print("=" * 72)
print(f"  Brinkman sensitivities:")
print(f"    d ln(R A)/d t_ox = {SENS_TOX:.3f} /nm")
print(f"    d ln(R A)/d phi  = {SENS_PHI:.3f} /eV")
print(f"  CV(R A) from t_ox, phi fluctuations      = {CV_RA*100:5.2f} %")
print(f"  CV(D) from sqrt(CV(R_P)^2 - CV(RA)^2)/2   = {CV_D*100:5.2f} %  "
      f"  (sigma(D) = {CV_D*80:.2f} nm at D=80 nm)")
print(f"  CV(V_mag) = sqrt((2 CV_D)^2 + CV_tf^2)    = {CV_V*100:5.2f} %")
print(f"  CV(H_k)   (from TMR proxy)                = {CV_HK*100:5.2f} %")
print(f"  CV(M_s)   (literature)                    = {CV_MS*100:5.2f} %")
print(f"  CV(Delta)                                  = "
      f"{CV_DELTA_PDK*100:5.2f} %  <<<  PDK baseline")
print("=" * 72)

# ─────────────────────────────────────────────────────────────────────────────
# 4. Model definitions (thin script-local aliases over the analysis library)
# ─────────────────────────────────────────────────────────────────────────────
def psw_nb(V, Delta, Vc0=VC0, tau0=TAU0, tw=TW):
    """NB P_sw without the (1-V/Vc0) clamp — see thesis note on V_th>V_c0."""
    x = Delta * (1.0 - V / Vc0)
    return 1.0 - np.exp(-(tw / tau0) * np.exp(-x))


def vth_nb(Delta, Vc0=VC0, tau0=TAU0, tw=TW):
    return float(nb_fit.vth_nb(tw, Delta, Vc0, tau0=tau0))


def sigmoid(V, Vth, beta):
    return 1.0 / (1.0 + np.exp(-beta * (V - Vth)))

# ─────────────────────────────────────────────────────────────────────────────
# 5. Monte Carlo D2D sweep
# ─────────────────────────────────────────────────────────────────────────────
# V_DENSE must span the experimental operating voltage range.
# At t_w = 0.75 ns the Sigmoid mid-point is ~894 mV and the transition
# extends into the ~780-1000 mV range, so V is extrapolated beyond V_c0.
V_DENSE = np.linspace(0.60, 1.10, 1500)
VTH_NB  = vth_nb(DELTA0)

# Include the PDK baseline value explicitly in the sweep.
# Extended range up to 60% so that the 95% and 90% margin thresholds
# for the sharper 0.75-ns operating point remain within the scan.
CV_SWEEP = np.array([0.00, 0.03, 0.05, 0.07, CV_DELTA_PDK, 0.10, 0.13, 0.15,
                     0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60])
CV_SWEEP = np.sort(np.unique(np.round(CV_SWEEP, 4)))
N_SAMPLES = 20000
RNG = np.random.default_rng(seed=42)

mean_curves = {}
beta_eff    = np.zeros_like(CV_SWEEP)
vth_eff     = np.zeros_like(CV_SWEEP)

print(f"\n{'='*72}")
print(f"  Monte Carlo D2D sweep  (N = {N_SAMPLES}, t_w = {TW} ns)")
print(f"{'='*72}")
print(f"  {'CV_Delta':>9} {'sigma_D':>9} {'V_th_eff':>11} {'beta_eff':>11}")
print("  " + "-" * 45)

for i, cv in enumerate(CV_SWEEP):
    P_mean = wafer_average_psw(
        V_DENSE, cv, DELTA0, VC0,
        tw_ns=TW, tau0_ns=TAU0,
        N_samples=N_SAMPLES, rng=RNG,
    )
    mean_curves[cv] = P_mean
    try:
        Vth_fit, beta_fit = fit_sigmoid_to_average(V_DENSE, P_mean)
        vth_eff[i], beta_eff[i] = Vth_fit, beta_fit
    except RuntimeError:
        vth_eff[i], beta_eff[i] = np.nan, np.nan
    print(f"  {cv:>9.4f} {cv*DELTA0:>9.3f} {vth_eff[i]*1e3:>10.1f}  "
          f"{beta_eff[i]:>10.3f}")
print("=" * 72)

# ─────────────────────────────────────────────────────────────────────────────
# 6. C2C empirical factor and D2D transfer function
# ─────────────────────────────────────────────────────────────────────────────
BETA_SINGLE_FIT = beta_eff[0]
F_FUN           = beta_eff / BETA_SINGLE_FIT
ETA_C           = BETA_MEAS / BETA_SINGLE_FIT

idx_pdk = int(np.argmin(np.abs(CV_SWEEP - CV_DELTA_PDK)))

print(f"\n  NB single-device fit beta                 = {BETA_SINGLE_FIT:.3f} V^-1")
print(f"  NB single-device analytic beta            = {BETA_NB_ANALYTIC:.3f} V^-1")
print(f"  Measured Sigmoid beta                     = {BETA_MEAS:.3f} V^-1")
print(f"  Empirical C2C factor  eta_c               = {ETA_C:.3f}")
print()
print(f"  PDK baseline CV_Delta                     = {CV_DELTA_PDK*100:.2f} %")
print(f"  F(CV_PDK)                                 = {F_FUN[idx_pdk]:.3f}")
beta_pdk = ETA_C * F_FUN[idx_pdk] * BETA_SINGLE_FIT
print(f"  Predicted beta at CV_PDK                  = {beta_pdk:.3f} V^-1")
print(f"  Prediction / measurement ratio            = {beta_pdk/BETA_MEAS:.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# 7. Figure  (4 panels)
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13.0, 9.8))
gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.28,
              left=0.07, right=0.96, top=0.90, bottom=0.08)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

# ════════════════════════════════════════════════════════════════════════════
# Top-left — CV(Delta) budget from PDK-Brinkman decomposition
# ════════════════════════════════════════════════════════════════════════════
contrib_labels = [
    r"Area part of $V_{\rm mag}$" "\n" r"(via $R_P$ + Brinkman)",
    r"Anisotropy $H_k$" "\n" r"(via TMR proxy)",
    r"Magnetisation $M_s$" "\n" r"(literature)",
    r"Thickness $t_f$" "\n" r"(MBE precision)",
]
contrib_cv   = np.array([2 * CV_D, CV_HK, CV_MS, CV_TF])
contrib_var  = contrib_cv**2
contrib_pct  = contrib_var / contrib_var.sum() * 100
contrib_cols = [THU_DEEP, THU_MID, THU_SOFT, THU_PALE]

y_pos = np.arange(len(contrib_labels))
bars = ax_a.barh(y_pos, contrib_cv * 100, color=contrib_cols,
                 edgecolor="white", linewidth=1.2, height=0.60)

for bar, cv, pct in zip(bars, contrib_cv, contrib_pct):
    width = bar.get_width()
    ax_a.text(width + 0.15, bar.get_y() + bar.get_height() / 2,
              rf"{cv*100:.1f}%  ({pct:.0f}% of $\mathrm{{CV}}^2$)",
              va="center", ha="left", color=CHARCOAL)

ax_a.axvline(CV_DELTA_PDK * 100, color=CRIMSON, lw=1.8, ls="--", zorder=4)
ax_a.text(CV_DELTA_PDK * 100 + 0.15, 0.45,
          rf"Total" "\n" rf"$\mathrm{{CV}}(\Delta) = {CV_DELTA_PDK*100:.1f}$%",
          color=CRIMSON, va="center", ha="left", fontweight="bold",
          linespacing=1.15)

ax_a.set_yticks(y_pos)
ax_a.set_yticklabels(contrib_labels)
ax_a.invert_yaxis()
ax_a.set_xlabel(r"Relative standard deviation $\mathrm{CV}$ (%)")
ax_a.set_title(
    r"$\mathrm{CV}(\Delta)$ variance budget from PDK mismatch parameters"
)
ax_a.set_xlim(0, 9.0)
ax_a.grid(axis="y", visible=False)

# ════════════════════════════════════════════════════════════════════════════
# Top-right — C2C-calibrated wafer-averaged Sigmoid curves
# ════════════════════════════════════════════════════════════════════════════
CV_PLOT = [0.0, CV_DELTA_PDK, 0.20, 0.40, 0.60]
cv_range_max = 0.60
shade = lambda cv: THU_LIGHT_CMAP(0.20 + 0.75 * cv / cv_range_max)

# Storage for residual inset
calib_curves_for_inset = {}

for cv in CV_PLOT:
    idx = int(np.argmin(np.abs(CV_SWEEP - cv)))
    # Calibrated wafer model: at CV=0 the curve must equal the measured
    # Sigmoid (V_th_meas, beta_meas).  For CV>0, D2D broadens beta by
    # factor F(CV) and shifts V_th by the MC delta-shift relative to CV=0.
    vth_shift = vth_eff[idx] - vth_eff[0]
    vth_cal   = VTH_MEAS + vth_shift
    beta_cal  = BETA_MEAS * F_FUN[idx]
    P_cal     = sigmoid(V_DENSE, vth_cal, beta_cal)
    calib_curves_for_inset[cv] = P_cal

    if np.isclose(cv, CV_DELTA_PDK):
        ax_b.plot(V_DENSE * 1e3, P_cal, color=AMBER, lw=2.6, zorder=8,
                  label=rf"$\mathrm{{CV}}_\Delta = {cv*100:.1f}$% "
                        r"(PDK baseline)")
    elif cv == 0.0:
        ax_b.plot(V_DENSE * 1e3, P_cal, color=CRIMSON, lw=1.8, ls=":",
                  label=rf"$\mathrm{{CV}}_\Delta = 0$ ")
    else:
        ax_b.plot(V_DENSE * 1e3, P_cal, color=shade(cv), lw=1.8,
                  label=rf"$\mathrm{{CV}}_\Delta = {int(cv*100)}$%")

V_fit = np.linspace(780, 1000, 400) * 1e-3
ax_b.plot(V_fit * 1e3, sigmoid(V_fit, VTH_MEAS, BETA_MEAS),
          color=TEAL, lw=2.4, ls=(0, (5, 2)), zorder=10,
          label=r"Measured Sigmoid")

ax_b.axhline(0.5, color=CHARCOAL, lw=0.5, ls=":", alpha=0.5)
ax_b.set_xlabel(r"Write voltage $|V_{\rm SOT}|$ (mV)")
ax_b.set_ylabel(r"Switching probability $P_{\rm sw}(V)$")
ax_b.set_title(
    r"C2C-calibrated wafer-averaged Sigmoid curves"
)
ax_b.set_xlim(800, 990)
ax_b.set_ylim(-0.03, 1.05)
ax_b.legend(loc="lower right")

# ── Inset: residual P_sw(V; CV) - P_sw(V; 0) to make broadening visible ───
# At t_w = 0.75 ns the operating point is so close to the deterministic
# regime (V_th > V_c0) that the bare curves visually overlap; the residual
# panel makes the millivolt-level broadening readable.
ax_b_ins = ax_b.inset_axes([0.08, 0.55, 0.42, 0.40])
P_ref = calib_curves_for_inset[0.0]
for cv in CV_PLOT:
    if cv == 0.0: continue
    diff = (calib_curves_for_inset[cv] - P_ref) * 100   # in %
    if np.isclose(cv, CV_DELTA_PDK):
        ax_b_ins.plot(V_DENSE * 1e3, diff, color=AMBER, lw=2.0,
                       label=rf"{cv*100:.1f}%")
    else:
        ax_b_ins.plot(V_DENSE * 1e3, diff, color=shade(cv), lw=1.4,
                       label=rf"{int(cv*100)}%")
ax_b_ins.axhline(0, color=CHARCOAL, lw=0.5, ls=":", alpha=0.5)
ax_b_ins.set_xlim(820, 970)
ax_b_ins.set_xlabel(r"$|V_{\rm SOT}|$ (mV)", labelpad=1)
ax_b_ins.set_ylabel(r"$\Delta P_{\rm sw}$ (%)", labelpad=1)
ax_b_ins.tick_params(axis="both", labelsize=7, pad=1.5)
ax_b_ins.set_title(r"Deviation from $\mathrm{CV}_\Delta = 0$", pad=3, fontsize=12)
ax_b_ins.legend(loc="upper right", fontsize=6,
                 title=r"$\mathrm{CV}_\Delta$", title_fontsize=10.5,
                 framealpha=0.92, handlelength=1.2, labelspacing=0.2,
                 borderpad=0.3, handletextpad=0.3, ncol=2)
ax_b_ins.grid(True, color="#DDD0E8", lw=0.4, ls="--", alpha=0.6)
for spine in ax_b_ins.spines.values():
    spine.set_linewidth(0.7)

# ════════════════════════════════════════════════════════════════════════════
# Bottom-left — D2D transfer function F(CV_Delta)
# ════════════════════════════════════════════════════════════════════════════
ax_c.plot(CV_SWEEP * 100, F_FUN, color=THU_DEEP, lw=2.0, marker="o",
          markersize=7, markerfacecolor="white", markeredgewidth=1.8,
          label=r"MC simulation (N = 20000)")
ax_c.axhline(1.0, color=CHARCOAL, lw=0.8, ls=":", alpha=0.6)
ax_c.fill_between(CV_SWEEP * 100, F_FUN, 1.0,
                   color=THU_TINT, alpha=0.6, zorder=0,
                   label=r"D2D broadening region")

ax_c.axvline(CV_DELTA_PDK * 100, color=AMBER, lw=1.2, ls="--", alpha=0.9)
ax_c.scatter([CV_DELTA_PDK * 100], [F_FUN[idx_pdk]], s=180, marker="*",
             facecolor=AMBER, edgecolor="white", linewidths=1.4, zorder=9,
             label=rf"PDK baseline  $\mathrm{{CV}}_\Delta={CV_DELTA_PDK*100:.1f}$%"
                   rf", $\mathcal{{F}}={F_FUN[idx_pdk]:.3f}$")

ax_c.set_xlabel(r"$\mathrm{CV}_\Delta$ (%)")
ax_c.set_ylabel(r"$\mathcal{F}(\mathrm{CV}_\Delta) = \beta_{\rm eff}/\beta_{\rm NB}^{\rm fit}$")
ax_c.set_title(
    r"D2D transfer function — $\beta_{\rm eff}$ vs $\mathrm{CV}_\Delta$"
)
ax_c.set_xlim(-1, 62)
ax_c.set_ylim(0.78, 1.04)
ax_c.legend(loc="lower left")

# ════════════════════════════════════════════════════════════════════════════
# Bottom-right — combined D2D + C2C prediction vs measurement
# ════════════════════════════════════════════════════════════════════════════
ax_d.axhline(BETA_SINGLE_FIT, color=NAVY, lw=1.8, ls=":",
             label=rf"NB single device (fit): $\beta_{{\rm NB}}^{{\rm fit}} = "
                   rf"{BETA_SINGLE_FIT:.2f}$ V$^{{-1}}$")

ax_d.plot(CV_SWEEP * 100, F_FUN * BETA_SINGLE_FIT,
          color=NAVY, lw=2.0, marker="s", markersize=6,
          markerfacecolor="white", markeredgewidth=1.5,
          label=r"Wafer NB: $\mathcal{F}(\mathrm{CV}_\Delta)\,\beta_{\rm NB}^{\rm fit}$")

ax_d.axhline(BETA_MEAS, color=TEAL, lw=1.8, ls=(0, (5, 2)),
             label=rf"Measured: $\beta_s = {BETA_MEAS:.2f}$ V$^{{-1}}$ "
                   rf"($\eta_c = {ETA_C:.2f}$)")

ax_d.plot(CV_SWEEP * 100, ETA_C * F_FUN * BETA_SINGLE_FIT,
          color=CRIMSON, lw=2.2, marker="o", markersize=7,
          markerfacecolor="white", markeredgewidth=1.8,
          label=r"Combined: $\eta_c\,\mathcal{F}(\mathrm{CV}_\Delta)\,\beta_{\rm NB}^{\rm fit}$")

ax_d.axvline(CV_DELTA_PDK * 100, color=AMBER, lw=1.2, ls="--", alpha=0.9)
ax_d.scatter([CV_DELTA_PDK * 100], [beta_pdk], s=180, marker="*",
             facecolor=AMBER, edgecolor="white", linewidths=1.4, zorder=9)
ax_d.annotate(rf"$\beta^{{\rm eff}} = {beta_pdk:.1f}$ V$^{{-1}}$" "\n"
              rf"({beta_pdk/BETA_MEAS*100:.1f}% of measured)",
              xy=(CV_DELTA_PDK*100, beta_pdk),
              xytext=(22, 55),
              color=AMBER, fontweight="bold",
              arrowprops=dict(arrowstyle="->", color=AMBER, lw=1.0))

ax_d.set_xlabel(r"$\mathrm{CV}_\Delta$ (%)")
ax_d.set_ylabel(r"Sigmoid slope $\beta$ (V$^{-1}$)")
ax_d.set_title(r"Combined prediction vs measurement")
ax_d.set_xlim(-1, 62)
# Log scale so the NB single-device band (~8 V^-1) and the calibrated
# wafer band (~45 V^-1) are both legible.
ax_d.set_yscale("log")
ax_d.set_ylim(6, 70)
ax_d.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax_d.set_yticks([8, 10, 20, 30, 40, 60])
ax_d.legend(loc="center")

fig.suptitle(
    r"Process variability impact on sMTJ probabilistic response",
    y=0.97
)

fig.savefig(OUTDIR + "Chapter02_local_16.png", dpi=300, bbox_inches="tight")
print(f"\nSaved  Chapter02_local_16.png  to {OUTDIR}")
plt.close(fig)

# ─────────────────────────────────────────────────────────────────────────────
# 8. Summary
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print("  Key numerical results")
print("=" * 72)
print(f"  PDK-Brinkman derived CV(Delta)           = {CV_DELTA_PDK*100:.2f} %")
print(f"    V_mag area term (geometry, from R_P)    = "
      f"{(2*CV_D)**2/CV_DELTA_PDK_SQ*100:5.1f} %")
print(f"    H_k     (interface, from TMR)           = "
      f"{CV_HK**2/CV_DELTA_PDK_SQ*100:5.1f} %")
print(f"    M_s     (bulk magnetisation)            = "
      f"{CV_MS**2/CV_DELTA_PDK_SQ*100:5.1f} %")
print(f"    t_f     (thickness part of V_mag)       = "
      f"{CV_TF**2/CV_DELTA_PDK_SQ*100:5.1f} %")
print()
print(f"  NB single-device fit beta                 = {BETA_SINGLE_FIT:.3f} V^-1")
print(f"  Measured beta                             = {BETA_MEAS:.3f} V^-1")
print(f"  Empirical C2C factor  eta_c               = {ETA_C:.3f}")
print()
print(f"  Monte Carlo at CV_Delta = {CV_DELTA_PDK*100:.1f} % (PDK baseline):")
print(f"    F(CV_PDK)                                = {F_FUN[idx_pdk]:.4f}")
print(f"    Predicted beta                           = {beta_pdk:.3f} V^-1")
print(f"    Prediction / measurement ratio           = {beta_pdk/BETA_MEAS:.4f}")
print()
print("  CV thresholds for beta_eff degradation:")
for target in [0.99, 0.95, 0.90]:
    if F_FUN.min() < target:
        cv_at = np.interp(-target, -F_FUN, CV_SWEEP)
        print(f"    beta_eff / beta_NB_fit >= {target:.2f}  ->  "
              f"CV_Delta <= {cv_at*100:.1f} %")
    else:
        print(f"    beta_eff / beta_NB_fit >= {target:.2f}  ->  "
              f"CV_Delta > {CV_SWEEP[-1]*100:.0f} % (outside scan)")
print("=" * 72)
