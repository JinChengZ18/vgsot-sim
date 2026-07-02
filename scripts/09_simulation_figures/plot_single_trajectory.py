"""
图 2.10 — vgsot-sim 在 t_w = 0.75 ns 写入脉冲下的单次 (m_z, R_MTJ, I_SOT)
演化事件可视化（与 §2.3.3 详细 P_sw 测试协议一致）。

三面板：
  (a) m_z(t)
  (b) R_MTJ(t) — 经 PDK TMR 模型从 m_z 直接换算（R_P ≈ 5 kΩ, R_AP ≈ 10 kΩ）
  (c) I_SOT(t) — 驱动脉冲

风格紧跟 04 章节 fig_material_params 范式（紧凑布局、字号合理、无重叠）。

输出路径：本目录 + 同步到 ../../article/00_chapter_drafts/figs/Chapter02_local_08.png
"""
from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation

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
    "axes.titlepad"        : 8,
    "legend.fontsize"      : 10.5,
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
    "lines.linewidth"      : 1.6,
    "legend.frameon"       : True,
    "legend.framealpha"    : 0.92,
    "legend.edgecolor"     : THU_PALE,
    "legend.facecolor"     : NEAR_WHITE,
    "figure.dpi"           : 150,
    "savefig.dpi"          : 300,
})

# ── Sim setup (matches §2.3.3 experimental protocol) ─────────────────────
cc = PhysicalConstantsConfig()
PULSE_NS = 0.75           # write pulse — matches detailed-P_sw measurement
RELAX_NS = 3.25           # relaxation tail (3 ns is enough for m to settle)
T_AMBIENT = 300.0

sim_end = int(round((PULSE_NS + RELAX_NS) * 1e-9 / cc.t_step))
mid1 = int(round(PULSE_NS * 1e-9 / cc.t_step))

# Four I_SOT spanning sub-/marginal/just-above/super-threshold. At the runtime
# θ_SH=0.066, integrated with the Cayley scheme (run_piecewise_direct_excitation's
# default — this script does not override integrator=), the pilot ensemble gives
# I_50 ≈ 1281 µA, consistent with the experimental V_th = 894 mV / R_W ≈ 1152 µA.
I_LIST_UA = [-600, -1100, -1300, -2000]
COLOURS   = [THU_PALE, THU_SOFT, THU_MID, THU_DEEP]
LABELS    = [rf"$I_{{\mathrm{{SOT}}}} = {i}\,\mu\mathrm{{A}}$" for i in I_LIST_UA]

SEED = 1     # representative seed (each I_SOT then exhibits its modal outcome
              # in the SER MC distribution — sub-threshold no-switch through
              # super-threshold switch — instead of a single rare draw)

results = []
for i_sot_uA in I_LIST_UA:
    np.random.seed(SEED)
    res = run_piecewise_direct_excitation(
        sim_start_step=1, sim_mid1_step=mid1, sim_mid2_step=sim_end, sim_end_step=sim_end,
        pap=1,
        v_mtj_stage1=0.0, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
        i_sot_stage1=i_sot_uA*1e-6, i_sot_stage2=0.0, i_sot_stage3=0.0,
        estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1, estt_stage3=0, esot_stage3=1,
        vnv=0, non=1, r_sot_fl_dl=0.83,
        constants=cc, show_progress=False,
        enable_self_heating=True, T_ambient_K=T_AMBIENT,
    )
    results.append(res)

t_ns = results[0].time_s * 1e9

# ── Figure (3 stacked panels, compact like fig_material_params) ──────────
fig, axes = plt.subplots(3, 1, figsize=(8.4, 7.2), sharex=True,
                          gridspec_kw={"hspace": 0.22})

# (a) m_z(t)
ax = axes[0]
for res, c, lab in zip(results, COLOURS, LABELS):
    ax.plot(t_ns, res.mz, color=c, lw=1.7, label=lab)
ax.axhline(0, color="gray", lw=0.5, ls=":")
ax.axvline(PULSE_NS, color=CRIMSON, lw=1.2, ls=(0, (4, 2)))
ax.text(PULSE_NS + 0.04, 1.05, "pulse off",
        color=CRIMSON, fontsize=11, fontweight="bold", ha="left", va="top")
ax.set_ylabel(r"$m_z(t)$")
ax.set_ylim(-1.15, 1.15)
ax.set_title(rf"Single-trajectory output, pure SOT, $t_p = {PULSE_NS:g}\,\mathrm{{ns}}$"
              rf" (+ {RELAX_NS:g} ns relax)")
ax.legend(loc="center right", frameon=True, ncol=1, fontsize=10)

# (b) R_MTJ(t)
ax = axes[1]
for res, c in zip(results, COLOURS):
    ax.plot(t_ns, res.r_mtj * 1e-3, color=c, lw=1.7)
ax.axvline(PULSE_NS, color=CRIMSON, lw=1.2, ls=(0, (4, 2)))
ax.set_ylabel(r"$R_{\mathrm{MTJ}}$ (k$\Omega$)")
ax.set_ylim(4.5, 10.5)

# (c) I_SOT(t)
ax = axes[2]
for res, c in zip(results, COLOURS):
    ax.plot(t_ns, res.i_sot * 1e6, color=c, lw=1.7)
ax.axvline(PULSE_NS, color=CRIMSON, lw=1.2, ls=(0, (4, 2)))
ax.set_ylabel(r"$I_{\mathrm{SOT}}$ ($\mu$A)")
ax.set_xlabel("Time (ns)")

for ax in axes:
    ax.grid(True, linestyle="--", linewidth=0.4, color=THU_GRID, alpha=0.7)

plt.tight_layout()

out_path = Path(__file__).resolve().parent / "Chapter02_local_08.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
chapter_fig = Path(__file__).resolve().parent.parent.parent / "article" / "figs" / "Chapter02_local_08.png"
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)
print(f"Saved {out_path}")
print(f"Synced {chapter_fig}")

# ── Summary ─────────────────────────────────────────────────────────────
print()
print("=" * 72)
print(f"  Single-trajectory summary (t_w = {PULSE_NS} ns, self-heating ON)")
print("=" * 72)
print(f"  {'I_SOT (µA)':>12} {'mz_init':>10} {'mz_end':>10} {'ΔT_peak':>10}  {'switched?':>11}")
print("  " + "-" * 60)
for i_uA, res in zip(I_LIST_UA, results):
    dT_peak = float(np.max(res.T_K) - T_AMBIENT)
    switched = (res.mz[0] * res.mz[-1] < 0) and abs(res.mz[-1]) > 0.5
    print(f"  {i_uA:>12d} {res.mz[0]:>+10.4f} {res.mz[-1]:>+10.4f} "
          f"{dT_peak:>9.1f}K  {'yes' if switched else 'no':>11}")
print("=" * 72)
