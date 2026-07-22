"""
sMTJ Device Characterisation from Experimental Hysteresis Data
==============================================================
Direct fit of the Néel-Brown model to real measured hysteresis loops
from the 300 mm wafer platform, producing a three-panel summary:

  top-left   Resistance-voltage hysteresis loops at four pulse widths
  top-right  Critical voltage V_c vs ln(t_w) with log-linear fit
  bottom     Extracted NB model joint probability map  P_sw(V, t_w)

Primary device : Device A (four pulse widths: 0.75, 1, 2, 5 ns)
Cross-check    : Device B (three pulse widths: 1, 2, 5 ns)
External field : H_x = 200 Oe in-plane

Changelog (v4)
--------------
1. (a)/(b)/(c) labels removed from titles; titles describe content directly.
2. Title font enlarged (titlesize 13, figure titlesize 15).
3. Legends repositioned to avoid overlap with data/curves.
4. Heatmap colormap replaced with a light white-to-Tsinghua-purple scale;
   contour lines darkened and annotation boxes lightened for consistency.
5. Internal device codes renamed: device 4 -> Device A, device 2 -> Device B.
"""

# ─────────────────────────────────────────────────────────────────────────────
# 0. Imports
# ─────────────────────────────────────────────────────────────────────────────
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings("ignore")

# Physics + NB-fit math come from the vgsot-sim analysis sub-package:
#   - nb_fit.extract_Vc / loglin_fit / nb_params / fit_direction
#   - nb_fit.psw_nb / vth_nb / tau_ret / beta_nb_analytic
from vgsot_sim.analysis import nb_fit

# ─────────────────────────────────────────────────────────────────────────────
# 1. Global style
# ─────────────────────────────────────────────────────────────────────────────
THU_DEEP   = "#660874"
THU_MID    = "#8B3A9E"
THU_SOFT   = "#A966BE"
THU_PALE   = "#C99FD4"
THU_TINT   = "#EFE0F7"
THU_GRID   = "#DDD0E8"
AMBER      = "#C47A00"
TEAL       = "#1A6B5A"
CRIMSON    = "#A82038"
NAVY       = "#1F5FA8"
CHARCOAL   = "#2B2B2B"
NEAR_WHITE = "#FFFFFF"

# ── Light white-to-THU-purple colormap (replaces the earlier dark diverging) ─
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
    "axes.titleweight"     : "normal",
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
    "xtick.color"          : CHARCOAL,
    "ytick.color"          : CHARCOAL,
    "xtick.direction"      : "in",
    "ytick.direction"      : "in",
    "xtick.major.size"     : 4,
    "ytick.major.size"     : 4,
    "xtick.minor.size"     : 2.5,
    "ytick.minor.size"     : 2.5,
    "xtick.minor.visible"  : True,
    "ytick.minor.visible"  : True,
    "lines.linewidth"      : 1.6,
    "lines.markersize"     : 6,
    "legend.frameon"       : True,
    "legend.framealpha"    : 0.92,
    "legend.edgecolor"     : THU_PALE,
    "legend.facecolor"     : NEAR_WHITE,
    "legend.handlelength"  : 1.8,
    "legend.handletextpad" : 0.5,
    "legend.labelspacing"  : 0.35,
    "legend.borderpad"     : 0.45,
    "figure.dpi"           : 150,
    "savefig.dpi"          : 300,
    "savefig.bbox"         : "tight",
    "figure.facecolor"     : NEAR_WHITE,
})

from pathlib import Path
# Hysteresis-loop data files live alongside the 07 process-broadening
# scripts (`device{2,4}pulse_width_*.txt`); see 05_experimental_raw_data/raw/MAPPING.md
# for the device-ID ↔ paper-label convention.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJ_ROOT  = SCRIPT_DIR.parent
DATADIR    = str(PROJ_ROOT / "07_process_variability") + "/"
OUTDIR     = str(SCRIPT_DIR) + "/"

# ─────────────────────────────────────────────────────────────────────────────
# 2. Data loading — internal codes renamed to public labels
# ─────────────────────────────────────────────────────────────────────────────
DEVICE_FILES = {
    "A": [
        (0.75, "device4pulse_width_0_750_ns_200_Oe.txt"),
        (1.00, "device4pulse_width_1_000_ns_200_Oe.txt"),
        (2.00, "device4pulse_width_2_000_ns_200_Oe.txt"),
        (5.00, "device4pulse_width_5_000_ns_200_Oe.txt"),
    ],
    "B": [
        (1.00, "device2pulse_width_1_000_ns_200_Oe.txt"),
        (2.00, "device2pulse_width_2_000_ns_200_Oe.txt"),
        (5.00, "device2pulse_width_5_000_ns_200_Oe.txt"),
    ],
}


def load_loop(path):
    V, R = np.loadtxt(path, unpack=True)
    return V, R


DATA, TABLE = {}, []
for dev, files in DEVICE_FILES.items():
    DATA[dev] = []
    for tw, fname in files:
        V, R = load_loop(os.path.join(DATADIR, fname))
        Vm, Vp = nb_fit.extract_Vc(V, R)
        DATA[dev].append((tw, V, R, Vm, Vp))
        TABLE.append((dev, tw, Vm, Vp))

print("=" * 72)
print("  Critical-voltage extraction from measured hysteresis loops")
print("=" * 72)
print(f"  {'Device':<10} {'t_w (ns)':>10} {'V_th- (V)':>12} {'V_th+ (V)':>12}")
print("  " + "-" * 50)
for dev, tw, vm, vp in TABLE:
    print(f"  Device {dev:<5} {tw:>10.3f} {vm:>12.3f} {vp:>12.3f}")
print("=" * 72)

# ─────────────────────────────────────────────────────────────────────────────
# 3. Log-linear fits and NB model parameters
# ─────────────────────────────────────────────────────────────────────────────
KBT  = 0.02585
TAU0 = 1.0


FITS = {}
print()
print("=" * 72)
print(f"  Log-linear and NB fits (tau_0 = {TAU0:.1f} ns)")
print("=" * 72)
for dev in DATA:
    tws = np.array([r[0] for r in DATA[dev]])
    Vps = np.array([r[4] for r in DATA[dev]])
    Vms = np.array([np.abs(r[3]) for r in DATA[dev]])
    fit_p = nb_fit.fit_direction(tws, Vps, tau0=TAU0)
    fit_m = nb_fit.fit_direction(tws, Vms, tau0=TAU0)
    FITS[dev] = {
        "tws": tws,
        "pos": dict(a=fit_p.a, b=fit_p.b, D=fit_p.Delta, Vc=fit_p.Vc0, V=Vps,
                     tret=fit_p.tau_ret_ns),
        "neg": dict(a=fit_m.a, b=fit_m.b, D=fit_m.Delta, Vc=fit_m.Vc0, V=Vms,
                     tret=fit_m.tau_ret_ns),
    }
    print(f"\n  Device {dev}")
    for lbl, d in [("P→AP", FITS[dev]["pos"]), ("AP→P", FITS[dev]["neg"])]:
        print(f"    {lbl}:  V = {d['a']:.3f} - {d['b']:.3f} ln(t/ns)")
        print(f"           Δ = {d['D']:.2f}, V_c0 = {d['Vc']*1e3:.0f} mV, "
              f"E_b = {d['D']*KBT*1e3:.0f} meV, τ_ret = {d['tret']:.0f} ns")
print("=" * 72)

# Re-export common NB-fit functions under the script's historic names so the
# plotting code below remains unchanged.
psw = nb_fit.psw_nb
vth_nb = nb_fit.vth_nb

# ─────────────────────────────────────────────────────────────────────────────
# 5. Main three-panel figure
# ─────────────────────────────────────────────────────────────────────────────
DEV_MAIN = "A"
MAIN = FITS[DEV_MAIN]
loops = DATA[DEV_MAIN]

fig = plt.figure(figsize=(12.5, 9.6))
gs  = GridSpec(2, 2, figure=fig,
               width_ratios=[1.0, 1.0], height_ratios=[1.0, 1.15],
               hspace=0.38, wspace=0.24,
               left=0.07, right=0.95, top=0.93, bottom=0.08)

ax_top_l = fig.add_subplot(gs[0, 0])
ax_top_r = fig.add_subplot(gs[0, 1])
ax_bot   = fig.add_subplot(gs[1, :])

# ════════════════════════════════════════════════════════════════════════════
# Top-left — hysteresis loops
# ════════════════════════════════════════════════════════════════════════════
LOOP_COLORS  = ["#D9B8E6", THU_SOFT, THU_MID, THU_DEEP]
LOOP_MARKERS = ["o", "s", "^", "D"]

for (tw, V, R, Vm, Vp), col, mk in zip(loops, LOOP_COLORS, LOOP_MARKERS):
    ax_top_l.plot(V, R * 1e-3,
                  color=col, lw=1.3, marker=mk, markersize=3.6,
                  markerfacecolor=col, markeredgecolor="white",
                  markeredgewidth=0.4,
                  label=rf"$t_w = {tw:g}$ ns")

ax_top_l.set_xlabel(r"Pulse voltage (V)")
ax_top_l.set_ylabel(r"Device resistance (k$\Omega$)")
ax_top_l.set_title(r"Resistance–voltage hysteresis at four pulse widths")
ax_top_l.set_xlim(-1.15, 1.15)
ax_top_l.set_ylim(4.2, 10.6)
ax_top_l.axvline(0, color=CHARCOAL, lw=0.5, linestyle=":", alpha=0.4)
# Legend placed in the hysteresis neutral zone (|V| < 0.35, 7–8 kΩ), a region
# never traversed by data for this device.  Compact 2-column layout.
ax_top_l.legend(loc="center", bbox_to_anchor=(0.5, 0.50),
                fontsize=11, ncol=2,
                title=rf"Device {DEV_MAIN}  ($H_x = 200$ Oe)",
                title_fontsize=11,
                columnspacing=0.9, handletextpad=0.4)

# ════════════════════════════════════════════════════════════════════════════
# Top-right — V_c vs ln(t_w)
# ════════════════════════════════════════════════════════════════════════════
T_CURVE = np.logspace(np.log10(0.5), np.log10(20), 400)

ax_top_r.scatter(MAIN["tws"], MAIN["pos"]["V"],
                 s=75, marker="o", facecolor="white",
                 edgecolor=CRIMSON, linewidths=1.8, zorder=5,
                 label=r"Exp. $V_{\rm th+}$ (P$\to$AP)")

a_p, b_p = MAIN["pos"]["a"], MAIN["pos"]["b"]
ax_top_r.plot(T_CURVE, a_p - b_p * np.log(T_CURVE),
              color=CRIMSON, lw=2.0,
              label=rf"Fit: $V(t)={a_p:.2f}-{b_p:.2f}\ln(t)$")

ax_top_r.scatter(MAIN["tws"], -MAIN["neg"]["V"],
                 s=75, marker="s", facecolor="white",
                 edgecolor=NAVY, linewidths=1.8, zorder=5,
                 label=r"Exp. $V_{\rm th-}$ (AP$\to$P)")

a_m, b_m = MAIN["neg"]["a"], MAIN["neg"]["b"]
ax_top_r.plot(T_CURVE, -(a_m - b_m * np.log(T_CURVE)),
              color=NAVY, lw=2.0,
              label=rf"Fit: $V(t)=-{a_m:.2f}+{b_m:.2f}\ln(t)$")

ax_top_r.axhline(0, color=CHARCOAL, lw=0.5, linestyle=":", alpha=0.5)
ax_top_r.set_xscale("log")
ax_top_r.set_xlabel(r"Pulse width $t_w$ (ns)")
ax_top_r.set_ylabel(r"$V_{\rm c}$ (V)")
ax_top_r.set_title(r"Critical voltage vs pulse width — log-linear fit")
ax_top_r.set_xlim(0.5, 20)
ax_top_r.set_ylim(-1.05, 1.05)
ax_top_r.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax_top_r.set_xticks([0.5, 1, 2, 5, 10, 20])
# Legend anchored at center-left (the two fit lines pass through ±0.8 V
# near the origin; the y≈0 center-left region is empty between them).
ax_top_r.legend(loc="center left", fontsize=10.5, framealpha=0.96)

D_p, Vc_p      = MAIN["pos"]["D"], MAIN["pos"]["Vc"]
D_m, Vc_m      = MAIN["neg"]["D"], MAIN["neg"]["Vc"]
tret_p, tret_m = MAIN["pos"]["tret"], MAIN["neg"]["tret"]

annot = (rf"$\tau_0 = {TAU0:.0f}$ ns assumed" "\n"
         rf"P$\to$AP: $\Delta={D_p:.2f}$, $V_{{c0}}={Vc_p*1e3:.0f}$ mV, "
         rf"$\tau_{{\rm ret}}={tret_p:.0f}$ ns" "\n"
         rf"AP$\to$P: $\Delta={D_m:.2f}$, $V_{{c0}}={Vc_m*1e3:.0f}$ mV, "
         rf"$\tau_{{\rm ret}}={tret_m:.0f}$ ns")
# Parameter box: bottom-right quadrant, below the AP→P line at long t_w.
ax_top_r.text(0.98, 0.02, annot, transform=ax_top_r.transAxes,
              fontsize=10, va="bottom", ha="right",
              bbox=dict(boxstyle="round,pad=0.35",
                        facecolor=THU_TINT, edgecolor=THU_PALE, alpha=0.96))

# ════════════════════════════════════════════════════════════════════════════
# Bottom — joint P_sw(V, t_w) heatmap  (light colormap, dark contours)
# ════════════════════════════════════════════════════════════════════════════
V_grid  = np.linspace(0.0, 1.10, 400)
tw_grid = np.logspace(np.log10(0.4), np.log10(50), 400)
VV, TT  = np.meshgrid(V_grid, tw_grid)
# Use the unclipped NB form for contour-smoothness in the heatmap; the
# canonical `psw_nb` clips x = max(Δ(1−V/V_c0), 0) at V_c0 to enforce the
# deterministic-switching regime, which creates a visible knee in iso-Psw
# contours at V = V_c0. Removing the clip gives smooth contours
# (asymptotically same value because exp(-(t/τ)·e^{−x}) → 1 fast for x<0).
PM = 1.0 - np.exp(-(TT / TAU0) * np.exp(-D_p * (1.0 - VV / Vc_p)))

im = ax_bot.pcolormesh(VV * 1e3, TT, PM, cmap=THU_LIGHT_CMAP,
                       vmin=0, vmax=1, shading="auto", rasterized=True)

LEVELS   = [0.05, 0.10, 0.30, 0.50, 0.70, 0.90, 0.95]
CLBL_FMT = {lv: f"{int(lv*100)}%" for lv in LEVELS}
# Dark purple contour lines — high contrast against the light backdrop,
# replacing the earlier white lines that relied on a dark heatmap.
# Contour lines: dark purple in the light (low-P) region automatically stays
# readable; in the dark (high-P) region the same colour would vanish, so we
# overlay white contour lines there with a cutoff at the 50% level.
CS_low = ax_bot.contour(VV * 1e3, TT, PM,
                        levels=[0.05, 0.10, 0.30, 0.50],
                        colors=[THU_DEEP], linewidths=0.95)
CS_high = ax_bot.contour(VV * 1e3, TT, PM,
                         levels=[0.70, 0.90, 0.95],
                         colors=["white"], linewidths=0.95)
ax_bot.clabel(CS_low, fmt={lv: f"{int(lv*100)}%" for lv in [0.05, 0.10, 0.30, 0.50]},
              fontsize=10.5, inline=True, inline_spacing=3)
ax_bot.clabel(CS_high, fmt={lv: f"{int(lv*100)}%" for lv in [0.70, 0.90, 0.95]},
              fontsize=10.5, inline=True, inline_spacing=3)

# 50%-probability locus — charcoal dashed line, visible on all backgrounds
t_loc = np.logspace(np.log10(0.4), np.log10(50), 600)
v_loc = vth_nb(t_loc, D_p, Vc_p)
mask  = v_loc > 0
ax_bot.plot(v_loc[mask] * 1e3, t_loc[mask],
            color=CHARCOAL, lw=2.2, ls="--", zorder=6,
            label=r"$P_{\rm sw}=50\%$ locus (NB model)")

for dev, mk, lbl in [("A", "o", r"Device A (P$\to$AP)"),
                     ("B", "^", r"Device B (P$\to$AP)")]:
    tws_d = FITS[dev]["tws"]
    Vs_d  = FITS[dev]["pos"]["V"]
    ax_bot.scatter(Vs_d * 1e3, tws_d,
                   s=70, c="white", marker=mk,
                   edgecolors=CRIMSON, linewidths=1.8, zorder=7,
                   label=lbl)

ax_bot.set_yscale("log")
ax_bot.set_xlabel(r"Write voltage $V_{\rm SOT}$ (mV)")
ax_bot.set_ylabel(r"Pulse width $t_w$ (ns)")
ax_bot.set_title(
    r"Joint switching probability $P_{\rm sw}(V,\,t_w)$ — "
    r"Néel-Brown model calibrated on Device A (P$\to$AP)"
)
ax_bot.set_xlim(0, 1080)
ax_bot.set_ylim(0.4, 50)
ax_bot.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax_bot.set_yticks([0.5, 1, 2, 5, 10, 20, 50])

# Legend placed in lower-left (empty low-probability region on the new cmap).
ax_bot.legend(loc="lower left", fontsize=11.5,
              framealpha=0.95, edgecolor=THU_PALE, ncol=1)

cbar = fig.colorbar(im, ax=ax_bot, fraction=0.035, pad=0.015)
cbar.set_label(r"$P_{\rm sw}$", fontsize=12.5)
cbar.ax.tick_params(labelsize=11)
cbar.set_ticks([0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0])

ann = (rf"$\Delta = {D_p:.2f}$, $V_{{c0}} = {Vc_p*1e3:.0f}$ mV" "\n"
       rf"$\tau_0 = {TAU0:.0f}$ ns,  $\tau_{{\rm ret}} = {tret_p:.0f}$ ns" "\n"
       rf"$\beta_s^{{\rm NB}} = {2*D_p*np.log(2)/Vc_p:.2f}$ V$^{{-1}}$")
# Parameter box anchored in the upper-left quadrant — the low-V, long-t_w
# region where no iso-probability contours run, and away from the legend.
ax_bot.text(0.02, 0.97, ann, transform=ax_bot.transAxes,
            fontsize=10.5, ha="left", va="top", color=CHARCOAL,
            bbox=dict(boxstyle="round,pad=0.4",
                      facecolor=NEAR_WHITE, edgecolor=THU_PALE, alpha=0.95))

fig.savefig(OUTDIR + "Chapter02_local_13.png", dpi=300, bbox_inches="tight")
print(f"\nSaved  Chapter02_local_13.png  to {OUTDIR}")
plt.close(fig)

# ─────────────────────────────────────────────────────────────────────────────
# 6. Device-to-device consistency figure
# ─────────────────────────────────────────────────────────────────────────────
fig2, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

ax = axes[0]
for dev, col, mk in [("A", CRIMSON, "o"),
                     ("B", THU_DEEP, "^")]:
    tws_d = FITS[dev]["tws"]
    Vs_d  = FITS[dev]["pos"]["V"]
    a, b  = FITS[dev]["pos"]["a"], FITS[dev]["pos"]["b"]
    D, Vc = FITS[dev]["pos"]["D"], FITS[dev]["pos"]["Vc"]
    ax.scatter(tws_d, Vs_d, s=70, marker=mk,
               facecolor="white", edgecolor=col, linewidths=1.9,
               label=rf"Device {dev}: $\Delta$={D:.2f}, $V_{{c0}}$={Vc*1e3:.0f} mV",
               zorder=5)
    ax.plot(T_CURVE, a - b * np.log(T_CURVE), color=col, lw=1.9, alpha=0.85)

ax.set_xscale("log")
ax.set_xlabel(r"Pulse width $t_w$ (ns)")
ax.set_ylabel(r"$V_{\rm th+}$ (V)  [P$\to$AP]")
ax.set_title(r"Device-to-device consistency (P$\to$AP)")
ax.set_xlim(0.5, 20)
ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.set_xticks([0.5, 1, 2, 5, 10, 20])
ax.legend(loc="upper right", fontsize=10.5)

ax = axes[1]
dev_labels = list(DATA.keys())
x_pos = np.arange(len(dev_labels))
widths = 0.28

params = {
    r"$\Delta$ (P$\to$AP)": [FITS[d]["pos"]["D"] for d in dev_labels],
    r"$\Delta$ (AP$\to$P)": [FITS[d]["neg"]["D"] for d in dev_labels],
}
colours = [CRIMSON, NAVY]

for i, (lbl, vals) in enumerate(params.items()):
    ax.bar(x_pos + (i - 0.5) * widths, vals, widths,
           color=colours[i], alpha=0.82, edgecolor="white",
           label=lbl)

for i, d in enumerate(dev_labels):
    ax.text(i, max(FITS[d]["pos"]["D"], FITS[d]["neg"]["D"]) + 0.25,
            rf"$V_{{c0}}\!=$ {FITS[d]['pos']['Vc']*1e3:.0f} mV",
            ha="center", fontsize=10, color=CHARCOAL)

ax.set_xticks(x_pos)
ax.set_xticklabels([f"Device {d}" for d in dev_labels])
ax.set_ylabel(r"Thermal stability factor $\Delta$")
ax.set_title(r"Extracted $\Delta$ per device  ($\tau_0 = 1$ ns)")
ax.set_ylim(0, 7.5)
ax.legend(fontsize=10.5, loc="upper right")
ax.grid(axis="x", visible=False)

fig2.savefig(OUTDIR + "Chapter02_local_15.png", dpi=300, bbox_inches="tight")
print(f"Saved  Chapter02_local_15.png  to {OUTDIR}")
plt.close(fig2)

# ─────────────────────────────────────────────────────────────────────────────
# 7. Final parameter summary
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print("  Final Parameter Summary (direct fit to measured hysteresis data)")
print("=" * 72)
print(f"  τ_0 = {TAU0:.1f} ns (literature prior)")
for dev in DATA:
    print(f"\n  Device {dev}:")
    for lbl, key in [("P→AP", "pos"), ("AP→P", "neg")]:
        d = FITS[dev][key]
        print(f"    {lbl}: a = {d['a']:.3f} V, b = {d['b']:.3f} V, "
              f"Δ = {d['D']:.2f}, V_c0 = {d['Vc']*1e3:.0f} mV, "
              f"E_b = {d['D']*KBT*1e3:.0f} meV, τ_ret = {d['tret']:.0f} ns")
print("=" * 72)
