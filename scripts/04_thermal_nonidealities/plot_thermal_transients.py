"""
Transient thermal evolution of an 80 nm VGSOT-MTJ device.

Refactored: physics now comes from
`vgsot_sim.temperature_cases.thermal_transient`, which evaluates the
lumped one-dimensional RC thermal network derived in §2.2.2.1:
  C_v * t_MTJ * dT/dt = Q_MTJ + Q_SOT - (lambda_MgO / t_MgO) (T - T_0)
The closed-form T(t) = T_0 + ΔT_eq (1 − exp(−t/τ_th)) and the
operating-mode comparison (STT-only / SOT-only / STT+SOT) are both
returned ready to plot.

Generates fig_04_thermal_transients.png with two stacked panels:
  (top)    three operating modes at V = 0.8 V
  (bottom) pure SOT mode, V_SOT swept from 0.4 V to 1.2 V
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.temperature_cases import thermal_transient
from vgsot_sim.thermal import thermal_time_constant

# ── Style ────────────────────────────────────────────────────────────────
THU_DEEP, THU_MID, THU_SOFT, THU_PALE, THU_TINT = \
    "#660874", "#8B3A9E", "#A966BE", "#C99FD4", "#EFE0F7"
THU_GRID = "#DDD0E8"
CRIMSON, NAVY, TEAL, AMBER = "#A82038", "#1F5FA8", "#1A6B5A", "#C47A00"
CHARCOAL, NEAR_WHITE = "#2B2B2B", "#FFFFFF"

_FAMILY = ["Arial", "Liberation Serif"]
plt.rcParams.update({
    "font.family"          : "serif",
    "font.serif"           : _FAMILY,
    "font.size"            : 12,
    "axes.labelsize"       : 13,
    "axes.titlesize"       : 14,
    "axes.titlepad"        : 7,
    "figure.titlesize"     : 16,
    "figure.titleweight"   : "bold",
    "legend.fontsize"      : 9.5,
    "xtick.labelsize"      : 12,
    "ytick.labelsize"      : 12,
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

# ── Physical evaluation through vgsot-sim ────────────────────────────────
cc = PhysicalConstantsConfig()
T0 = 300.0
t_ps = np.linspace(0.0, 100.0, 600)
t_seconds = t_ps * 1e-12

tau_th_ps = thermal_time_constant(cc) * 1e12

# Three operating modes at 0.8 V
modes = []
for label_tex, vmtj, vsot, colour, dash in [
    (r'STT only ($V_\mathrm{MTJ}{=}0.8$ V)', 0.8, 0.0, '#888888', (4, 2)),
    (r'SOT only ($V_\mathrm{SOT}{=}0.8$ V)', 0.0, 0.8, THU_MID,   (5, 2, 1, 2)),
    (r'STT+SOT (0.8 V / 0.8 V)',             0.8, 0.8, THU_DEEP,  None),
]:
    r = thermal_transient(cc, V_MTJ=vmtj, V_SOT=vsot, T_0=T0, t_seconds=t_seconds, label=label_tex)
    modes.append((r, colour, dash))

# SOT-only voltage scan
V_list = [0.4, 0.6, 0.8, 1.0, 1.2]
sot_curves = [
    thermal_transient(cc, V_MTJ=0.0, V_SOT=v, T_0=T0, t_seconds=t_seconds,
                       label=fr'$V_\mathrm{{SOT}}{{=}}{v:.1f}$ V')
    for v in V_list
]

# ── Figure (1×2 side-by-side, compact, no overlaps) ──────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4),
                          gridspec_kw={"wspace": 0.28})

# (a) three operating modes -------------------------------------------
ax = axes[0]
for r, colour, dash in modes:
    label = rf'{r.label}, $T_\mathrm{{eq}}{{=}}{r.T_eq_K:.1f}$ K'
    if dash is None:
        ax.plot(t_ps, r.T_K, '-', color=colour, linewidth=1.6, label=label)
    else:
        ax.plot(t_ps, r.T_K, color=colour, linewidth=1.6, dashes=dash, label=label)

ax.axvline(tau_th_ps, color=THU_DEEP, linestyle=':', linewidth=1.2)
# Bold τ_th annotation in the upper-left clear zone, boxed for prominence
ax.annotate(rf'$\tau_\mathrm{{th}}\approx {tau_th_ps:.1f}$ ps',
            xy=(tau_th_ps, 326), xytext=(tau_th_ps + 5, 340),
            fontsize=12.5, fontweight='bold', color=THU_DEEP,
            arrowprops=dict(arrowstyle='->', color=THU_DEEP, lw=1.1),
            bbox=dict(boxstyle="round,pad=0.25", fc="#FFFFFF",
                      ec=THU_DEEP, lw=0.9, alpha=0.95))
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Device temperature (K)')
ax.set_title('Three operating modes (0.8 V drive)')
# Legend in lower right — single column to avoid overlap with curves
ax.legend(frameon=True, loc='lower right', fontsize=9.5,
          framealpha=0.95, edgecolor=THU_PALE)
ax.grid(True, linestyle='--', linewidth=0.4, color=THU_GRID, alpha=0.8)
ax.set_xlim(0, 100)
# leave headroom for the τ_th label above the highest curve (~341 K)
ax.set_ylim(298, 348)

# (b) SOT-only swept ------------------------------------------------
colours = plt.cm.Purples(np.linspace(0.35, 0.92, len(V_list)))
ax = axes[1]
for r, colour in zip(sot_curves, colours):
    ax.plot(t_ps, r.T_K, color=colour, linewidth=1.6,
            label=rf'{r.label} ($T_\mathrm{{eq}}{{=}}{r.T_eq_K:.1f}$ K)')

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Device temperature (K)')
ax.set_title(r'Pure-SOT mode: $V_\mathrm{SOT}$ scan')
# Legend in upper-left where the rising curves do not reach (T < ~330 K)
ax.legend(frameon=True, loc='upper left', fontsize=9.5,
          framealpha=0.95, edgecolor=THU_PALE)
ax.grid(True, linestyle='--', linewidth=0.4, color=THU_GRID, alpha=0.8)
ax.set_xlim(0, 100)
# top curve hits ~363 K, give a little headroom so legend can sit inside
ax.set_ylim(298, 375)

plt.tight_layout()
out_path = str(Path(__file__).resolve().parent / 'fig_04_thermal_transients.png')
plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out_path}')
print(f'tau_th = {tau_th_ps:.2f} ps')
for r, *_ in modes:
    print(f'  {r.label}: T_eq = {r.T_eq_K:.2f} K (ΔT = {r.T_eq_K - T0:.2f} K)')
for r in sot_curves:
    print(f'  V_SOT = {r.label[-7:-5]} V: T_eq = {r.T_eq_K:.2f} K (ΔT = {r.T_eq_K - T0:.2f} K)')
