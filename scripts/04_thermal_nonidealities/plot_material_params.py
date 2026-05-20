"""
Temperature-dependent material parameters and TMR bias response for an 80 nm
VGSOT-MTJ device.

Refactored: physics now comes from vgsot-sim
  - `vgsot_sim.temperature_cases.sweep_material_temperature` — M_s(T), K_i(T), η(T)
  - `vgsot_sim.temperature_cases.tmr_voltage_sweep` — TMR(V) under PDK / Lorentzian

This script only carries layout, labels and style; the math lives in the
library so it can be regenerated from the same code path used in the chapter
simulations.

Generates fig_05_material_params.png as a 2x2 panel:
  (a) Saturation magnetization  M_s(T)  — Bloch T^{3/2} law
  (b) Interfacial anisotropy     K_i(T)  — modified Callen-Callen (exponent 2.18)
  (c) Spin polarization ratio    eta(T)  — follows M_s(T)
  (d) TMR ratio vs bias voltage — two compact-model forms compared
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.temperature_cases import (
    sweep_material_temperature,
    tmr_voltage_sweep,
)

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
    "axes.titlepad"        : 9,
    "figure.titlesize"     : 18,
    "figure.titleweight"   : "bold",
    "legend.fontsize"      : 10,
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

# ── Run the physics through vgsot-sim ────────────────────────────────────
cc = PhysicalConstantsConfig()             # 80 nm SOT-MTJ baseline parameters
mat = sweep_material_temperature(cc, T_min=250.0, T_max=500.0, n_points=400)
tmr_curve = tmr_voltage_sweep(cc, V_min=0.0, V_max=1.5, n_points=400)

T      = mat.T_K
T_RT   = cc.T_RT
Ms_RT  = cc.Ms
Ki_RT  = cc.Ki
eta_RT = cc.P

V       = tmr_curve.V
TMR_pdk = np.maximum(tmr_curve.TMR_pdk, 0.0) * 100.0   # percent
TMR_lor = tmr_curve.TMR_lor * 100.0                      # percent
TMR_0   = cc.TMR
V_h     = cc.Vh

# ------------------------- Figure -------------------------------
fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.2))

def style(ax, title):
    ax.set_title(title)
    ax.grid(True, linestyle='--', linewidth=0.4, color=THU_GRID, alpha=0.8)

# --- (a) M_s(T) ---
ax = axes[0, 0]
ax.plot(T, mat.M_s / 1e3, color=THU_DEEP, linewidth=1.6)
ax.plot([T_RT], [Ms_RT / 1e3], 'o', color=THU_DEEP, markersize=5, zorder=3)
ax.annotate('RT (300 K)', xy=(T_RT, Ms_RT / 1e3),
            xytext=(T_RT + 12, Ms_RT / 1e3 + 5))
ax.set_xlabel('Temperature (K)')
ax.set_ylabel(r'$M_s$ (kA/m)')
style(ax, r'Saturation magnetization $M_s(T)$')

# --- (b) K_i(T) ---
ax = axes[0, 1]
ax.plot(T, mat.K_i * 1e3, color=THU_DEEP, linewidth=1.6)
ax.plot([T_RT], [Ki_RT * 1e3], 'o', color=THU_DEEP, markersize=5, zorder=3)
ax.annotate('RT (300 K)', xy=(T_RT, Ki_RT * 1e3),
            xytext=(T_RT + 12, Ki_RT * 1e3 + 0.006))
ax.set_xlabel('Temperature (K)')
ax.set_ylabel(r'$K_i$ (mJ/m$^2$)')
style(ax, r'Interfacial anisotropy $K_i(T)$')

# --- (c) eta(T) ---
ax = axes[1, 0]
ax.plot(T, mat.eta, color=THU_DEEP, linewidth=1.6)
ax.plot([T_RT], [eta_RT], 'o', color=THU_DEEP, markersize=5, zorder=3)
ax.annotate('RT (300 K)', xy=(T_RT, eta_RT),
            xytext=(T_RT + 12, eta_RT + 0.005))
ax.set_xlabel('Temperature (K)')
ax.set_ylabel(r'$\eta$')
style(ax, r'Spin polarization ratio $\eta(T)$')

# --- (d) TMR(V) — two compact-model forms compared ---
ax = axes[1, 1]
ax.plot(V, TMR_pdk, color=THU_DEEP, linewidth=1.7,
        label='Quadratic-rational (PDK)')
ax.plot(V, TMR_lor, color=THU_MID, linewidth=1.7, linestyle='--',
        label=r'Lorentzian, $V_h=0.5$ V')

# Mark V_h and TMR_0/2 on the Lorentzian curve
ax.axhline(TMR_0 * 100 / 2, color=THU_GRID, linewidth=0.6, linestyle=':')
ax.axvline(V_h,              color=THU_GRID, linewidth=0.6, linestyle=':')
ax.annotate(r'$V_h$', xy=(V_h, 5), color='#555555')
ax.annotate(r'$\mathrm{TMR}_0/2$',
            xy=(1.35, TMR_0 * 100 / 2 + 3),
            color='#555555', ha='right')

ax.set_xlabel(r'MTJ bias voltage $V_\mathrm{MTJ}$ (V)')
ax.set_ylabel('TMR (%)')
ax.set_xlim(0, 1.5)
# y-limit auto-scales to the calibrated TMR_0 (now 100% to match experiment)
ax.set_ylim(-5, TMR_0 * 100.0 * 1.15)
ax.legend(frameon=False, loc='upper right')
style(ax, 'TMR ratio vs bias voltage')

plt.tight_layout()
out_path = str(Path(__file__).resolve().parent / 'fig_05_material_params.png')
plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out_path}')

# ------------------------- Numeric comparison -------------------
from vgsot_sim.tmr import tmr_eff
print()
print(f"{'V (V)':>7} | {'TMR_PDK (%)':>12} | {'TMR_Lor (%)':>12} | diff (pp)")
print('-' * 50)
for v in [0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50]:
    cc.tmr_model = "pdk"
    tp = max(tmr_eff(v, cc), 0.0) * 100.0
    cc.tmr_model = "lorentzian"
    tl = tmr_eff(v, cc) * 100.0
    print(f"{v:7.2f} | {tp:12.2f} | {tl:12.2f} | {tp - tl:+8.2f}")
