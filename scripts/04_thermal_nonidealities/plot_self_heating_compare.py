"""
图 2.14 — 自热反馈下的磁化轨迹对比（self-heating ON vs OFF），紧凑 2×2 布局。

驱动 vgsot-sim 的时间步进求解器两次（相同 I_SOT、V_MTJ、脉宽、环境温度、
RNG 种子），一次 `enable_self_heating=False`，一次 `enable_self_heating=True`。
ON 支每步推进 RC 热扩散方程更新 T(t)，并以 T-修正后的 M_s(T)、K_i(T)
喂入各向异性场。两条轨迹共享同一确定性噪声实现，二者之差完全由温度反馈引入。

输出：fig_self_heating_trajectory.png — 紧凑 2×2 面板：
  (a) m_z(t) 对比 + (b) T(t)
  (c) R_MTJ(t) 对比 + (d) 材料漂移 ΔM_s/M_s、ΔK_i/K_i (%)

工作点：
  pap=1（初始平行态 m_z ≈ -1），V_MTJ=0（纯 SOT），
  I_SOT=-1500 µA  → V_SOT≈1.16 V → ΔT_eq≈59 K  → ΔK_i/K_i≈-11%
  t_pulse=3 ns，t_relax=5 ns（pulse 关断后让 T 衰减回环境温度）
  vnv=0, non=0（确定性比较，无热噪声）
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation
from vgsot_sim.thermal import (
    thermal_time_constant, steady_state_temperature, q_mtj, q_sot,
)
from vgsot_sim.initialize import compute_Rp

# ── Style (matches fig_material_params.py reference) ─────────────────────
THU_DEEP, THU_MID, THU_SOFT, THU_PALE = "#660874", "#8B3A9E", "#A966BE", "#C99FD4"
THU_GRID = "#DDD0E8"
CRIMSON, NAVY, TEAL, AMBER = "#A82038", "#1F5FA8", "#1A6B5A", "#C47A00"
CHARCOAL, NEAR_WHITE = "#2B2B2B", "#FFFFFF"

plt.rcParams.update({
    "font.family"          : "serif",
    "font.serif"           : ["Arial", "Liberation Serif"],
    "font.size"            : 12,
    "axes.labelsize"       : 13,
    "axes.titlesize"       : 14,
    "axes.titlepad"        : 7,
    "legend.fontsize"      : 10,
    "xtick.labelsize"      : 12,
    "ytick.labelsize"      : 12,
    "mathtext.fontset"     : "stix",
    "axes.linewidth"       : 0.9,
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
    "lines.linewidth"      : 1.7,
    "legend.frameon"       : True,
    "legend.framealpha"    : 0.92,
    "legend.edgecolor"     : THU_PALE,
    "figure.dpi"           : 150,
    "savefig.dpi"          : 300,
})

# ── Sim setup ────────────────────────────────────────────────────────────
cc = PhysicalConstantsConfig()
PULSE_NS, RELAX_NS = 3.0, 5.0
I_SOT = -1500e-6
V_MTJ = 0.0
T0 = 300.0

sim_end = int(round((PULSE_NS + RELAX_NS) * 1e-9 / cc.t_step))
mid1 = int(round(PULSE_NS * 1e-9 / cc.t_step))

common = dict(
    sim_start_step=1, sim_mid1_step=mid1, sim_mid2_step=sim_end, sim_end_step=sim_end,
    pap=1,
    v_mtj_stage1=V_MTJ, v_mtj_stage2=V_MTJ, v_mtj_stage3=V_MTJ,
    i_sot_stage1=I_SOT, i_sot_stage2=0.0, i_sot_stage3=0.0,
    estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1, estt_stage3=0, esot_stage3=1,
    vnv=0, non=0,
    r_sot_fl_dl=0.83,
    constants=cc,
    show_progress=False,
)
np.random.seed(2026); res_off = run_piecewise_direct_excitation(**common, enable_self_heating=False)
np.random.seed(2026); res_on  = run_piecewise_direct_excitation(**common, enable_self_heating=True, T_ambient_K=T0)

Rp_RT = compute_Rp(cc)
Q_steady = q_mtj(V_MTJ, Rp_RT, cc) + q_sot(I_SOT * cc.R_W, cc, R_SOT=cc.R_W)
T_eq_pred = steady_state_temperature(Q_steady, T0, cc)
tau_th_ps = thermal_time_constant(cc) * 1e12

t_ns = res_off.time_s * 1e9

# ── Figure: 2×2 compact (enlarged hspace so row titles never overlap) ───
fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(
    2, 2, figsize=(11.0, 7.4),
    gridspec_kw={"wspace": 0.28, "hspace": 0.50},
)


# Panel-label helper: places a single "(x)" tag in the upper-left corner
# without baking the letter into the panel title. Lets us reorder / rename
# panels later by editing only the title string, never the (a)/(b)/(c)
# enumeration. Article caption references the same abstract enumeration,
# so labels stay portable across script iterations.
def panel_tag(ax, tag, *, color=CHARCOAL, fontsize=13):
    ax.text(-0.10, 1.06, f"({tag})", transform=ax.transAxes,
            ha="left", va="top", fontsize=fontsize, fontweight="bold",
            color=color)


# (a) m_z(t)
ax_a.plot(t_ns, res_off.mz, color=CHARCOAL, lw=1.7, label="self-heating OFF")
ax_a.plot(t_ns, res_on.mz,  color=CRIMSON,  lw=1.7, ls="--", label="self-heating ON")
ax_a.axvline(PULSE_NS, color="gray", lw=0.7, ls=":")
ax_a.set_xlabel("Time (ns)")
ax_a.set_ylabel(r"$m_z(t)$")
ax_a.set_ylim(-1.15, 1.15)
ax_a.set_title("Magnetisation trajectory")
ax_a.legend(loc="lower right", frameon=True)

# (b) T(t)
ax_b.plot(t_ns, res_on.T_K, color=THU_DEEP, lw=1.7, label=r"$T(t)$ (stepper)")
ax_b.axhline(T_eq_pred, color=AMBER, lw=1.2, ls=":",
             label=rf"$T_{{\mathrm{{eq}}}} = {T_eq_pred:.1f}\,\mathrm{{K}}$")
ax_b.axhline(T0, color="gray", lw=0.7, ls=":")
ax_b.axvline(PULSE_NS, color="gray", lw=0.7, ls=":")
ax_b.text(0.05, T_eq_pred - 1.5,
          rf"$\tau_{{\mathrm{{th}}}} \approx {tau_th_ps:.1f}\,\mathrm{{ps}}$",
          color=THU_DEEP, fontsize=11, va="top")
ax_b.set_xlabel("Time (ns)")
ax_b.set_ylabel(r"$T$ (K)")
ax_b.set_title("Device temperature")
ax_b.legend(loc="center right", frameon=True)

# (c) R_MTJ(t)
ax_c.plot(t_ns, res_off.r_mtj * 1e-3, color=CHARCOAL, lw=1.7, label="OFF")
ax_c.plot(t_ns, res_on.r_mtj  * 1e-3, color=CRIMSON,  lw=1.7, ls="--", label="ON")
ax_c.axvline(PULSE_NS, color="gray", lw=0.7, ls=":")
ax_c.set_xlabel("Time (ns)")
ax_c.set_ylabel(r"$R_{\mathrm{MTJ}}$ (k$\Omega$)")
ax_c.set_title("Resistance trajectory")
ax_c.legend(loc="upper right", frameon=True)

# (d) Material drift
Ms_drift = (res_on.Ms_T - cc.Ms) / cc.Ms * 100.0
Ki_drift = (res_on.Ki_T - cc.Ki) / cc.Ki * 100.0
ax_d.plot(t_ns, Ms_drift, color=NAVY, lw=1.7, label=r"$\Delta M_s/M_s$")
ax_d.plot(t_ns, Ki_drift, color=TEAL, lw=1.7, label=r"$\Delta K_i/K_i$")
ax_d.axhline(0.0, color="gray", lw=0.6)
ax_d.axvline(PULSE_NS, color="gray", lw=0.7, ls=":")
ax_d.set_xlabel("Time (ns)")
ax_d.set_ylabel("Drift (%)")
ax_d.set_title("Material-parameter drift")
ax_d.legend(loc="lower right", frameon=True)

for ax in (ax_a, ax_b, ax_c, ax_d):
    ax.grid(True, linestyle="--", linewidth=0.4, color=THU_GRID, alpha=0.7)

plt.suptitle(rf"Pure-SOT switching, $I_{{\mathrm{{SOT}}}} = {I_SOT*1e6:.0f}\,\mu\mathrm{{A}}$, "
             rf"$V_{{\mathrm{{MTJ}}}} = {V_MTJ:.1f}\,\mathrm{{V}}$, "
             rf"$t_p = {PULSE_NS:.0f}\,\mathrm{{ns}} + t_{{\mathrm{{relax}}}} = {RELAX_NS:.0f}\,\mathrm{{ns}}$",
             fontsize=14, y=0.995)
plt.subplots_adjust(top=0.92)

out_path = Path(__file__).resolve().parent / "fig_self_heating_trajectory.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved {out_path}")

print()
print("=" * 70)
print("  Self-heating impact summary (pure-SOT operating point)")
print("=" * 70)
print(f"  Predicted τ_th       = {tau_th_ps:.2f} ps")
print(f"  Predicted T_eq       = {T_eq_pred:.2f} K   (ΔT = {T_eq_pred-T0:.2f} K)")
idx_peak = int(np.argmax(res_on.T_K))
print(f"  Stepper peak T       = {res_on.T_K[idx_peak]:.2f} K @ t = {t_ns[idx_peak]:.2f} ns")
print(f"  ΔM_s/M_s peak        = {(res_on.Ms_T[idx_peak]-cc.Ms)/cc.Ms*100:+.2f}%")
print(f"  ΔK_i/K_i peak        = {(res_on.Ki_T[idx_peak]-cc.Ki)/cc.Ki*100:+.2f}%")
print(f"  Net ΔH_PMA/H_PMA     ≈ "
      f"{((res_on.Ki_T[idx_peak]-cc.Ki)/cc.Ki - (res_on.Ms_T[idx_peak]-cc.Ms)/cc.Ms) * 100:+.2f}%")
print(f"  m_z @ end (OFF / ON) = {res_off.mz[-1]:+.4f} / {res_on.mz[-1]:+.4f}")
print("=" * 70)
