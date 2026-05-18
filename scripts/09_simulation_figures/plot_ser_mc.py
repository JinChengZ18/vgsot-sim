"""
图 2.11 — vgsot-sim 在 t_w = 0.75 ns 写入脉冲下 **P_sw** vs |I_SOT| 蒙特卡罗
扫描（P_sw = 1 − SER，与 §2.3.3 实测 Sigmoid P_sw 曲线在同一量纲上对比）。

单面板（inset 嵌套格式）：主图为宽范围 (300-3500 µA) 扫描覆盖三段；
inset 子图聚焦阈值附近 (800-1400 µA)，在实验 I_th=1152 µA 附近密集取点、
两侧稀疏，以还原 Sigmoid 形状。

可通过命令行参数 `--metric=ser` 切换为旧的 SER 表示。

风格沿用 fig_material_params 范式（紧凑尺寸 + 合理字号 + 无重叠）。
"""
import argparse
from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--metric", choices=("psw", "ser"), default="psw",
                    help="y-axis metric: 'psw' (default, switching probability) or 'ser' (legacy)")
parser.add_argument("--rng-mode", choices=("legacy", "generator"), default="legacy",
                    help="MC trial seeding: 'legacy' (np.random.seed, default — matches "
                         "chapter figures) or 'generator' (np.random.default_rng, cleaner, "
                         "different bit-stream but identical statistics in the limit)")
parser.add_argument("--integrator", choices=("euler_spherical", "cayley"),
                    default="euler_spherical",
                    help="LLG step type: 'euler_spherical' (default, matches the θ_SH=0.04 "
                         "experimental calibration) or 'cayley' (norm-preserving vector form, "
                         "for precision studies; threshold ~10%% higher than Euler)")
ARGS, _ = parser.parse_known_args()
METRIC    = ARGS.metric          # "psw" or "ser"
RNG_MODE  = ARGS.rng_mode        # "legacy" or "generator"
INTEGRATOR = ARGS.integrator     # "euler_spherical" or "cayley"

from vgsot_sim.configs import PhysicalConstantsConfig, SerSotNoVcmaThermalConfig
from vgsot_sim.ser_cases import ser_sot_no_vcma_thermal

# ── Style ────────────────────────────────────────────────────────────────
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
    "legend.fontsize"      : 10,
    "xtick.labelsize"      : 12,
    "ytick.labelsize"      : 12,
    "mathtext.fontset"     : "stix",
    "axes.linewidth"       : 0.9,
    "axes.edgecolor"       : CHARCOAL,
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
    "figure.dpi"           : 150,
    "savefig.dpi"          : 300,
})

cc = PhysicalConstantsConfig()
PULSE_NS, TOTAL_NS = 0.75, 4.0
sim_end = int(round(TOTAL_NS * 1e-9 / cc.t_step))
mid1    = int(round(PULSE_NS * 1e-9 / cc.t_step))

# Wide-range sweep (main panel)
I_SOT_LIST = np.array([-300, -600, -900, -1100, -1300, -1500, -1800,
                        -2200, -2800, -3500], dtype=float) * 1e-6

# Threshold-region fine sweep (inset). Dense around exp I_th = 1152 µA,
# sparser on either side, so the sigmoid drop is well sampled.
I_SOT_INSET = np.array([-800, -900, -1000, -1080, -1120, -1140, -1160, -1180,
                         -1220, -1280, -1350, -1400], dtype=float) * 1e-6

N_TRIALS = 80
N_TRIALS_FINE = 80          # bumped so inset CI is comparable to main

cfg = SerSotNoVcmaThermalConfig(
    i_sot_list=tuple(I_SOT_LIST), trials=N_TRIALS,
    sim_start_step=1, sim_mid1_step=mid1, sim_end_step=sim_end,
    pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
    target_mz=1.0, failure_tol=0.2, constants=cc,
)

print(f"Running SER MC at t_w = {PULSE_NS} ns "
      f"(rng_mode={RNG_MODE}, integrator={INTEGRATOR}) ...")
print("[wide / self-heating OFF]")
res_off = ser_sot_no_vcma_thermal(cfg, show_progress=False,
                                   enable_self_heating=False, seed=2026,
                                   rng_mode=RNG_MODE, integrator=INTEGRATOR)
print("[wide / self-heating ON]")
res_on  = ser_sot_no_vcma_thermal(cfg, show_progress=False,
                                   enable_self_heating=True, T_ambient_K=300.0,
                                   seed=2026,
                                   rng_mode=RNG_MODE, integrator=INTEGRATOR)

print(f"[fine sweep around exp I_th=1152 µA / {len(I_SOT_INSET)} × {N_TRIALS_FINE} runs]")
cfg_fine = SerSotNoVcmaThermalConfig(
    i_sot_list=tuple(I_SOT_INSET), trials=N_TRIALS_FINE,
    sim_start_step=1, sim_mid1_step=mid1, sim_end_step=sim_end,
    pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
    target_mz=1.0, failure_tol=0.2, constants=cc,
)
res_fine = ser_sot_no_vcma_thermal(cfg_fine, show_progress=False,
                                     enable_self_heating=True,
                                     T_ambient_K=300.0, seed=4096,
                                     rng_mode=RNG_MODE, integrator=INTEGRATOR)

# ── Wilson 95% CI ────────────────────────────────────────────────────────
def wilson(p, n, z=1.96):
    p = np.clip(p, 0, 1)
    denom = 1 + z**2/n
    centre = (p + z**2/(2*n))/denom
    half   = z*np.sqrt(p*(1-p)/n + z**2/(4*n**2))/denom
    return centre - half, centre + half

# Map SER → display metric. When METRIC=="psw", y is the *switching*
# probability P_sw = 1 − SER, which is what §2.3.3 reports directly. The
# Wilson interval is symmetric under the transform p ↔ 1 − p, so [lo, hi]
# of P_sw is just [1 − hi_SER, 1 − lo_SER].
def to_y(ser):                     # SER array → y-axis array
    return (1.0 - np.asarray(ser)) if METRIC == "psw" else np.asarray(ser)

def to_lohi(lo, hi):
    if METRIC == "psw":
        return 1.0 - np.asarray(hi), 1.0 - np.asarray(lo)
    return lo, hi

lo_off, hi_off = wilson(res_off.ser, N_TRIALS)
lo_on,  hi_on  = wilson(res_on.ser,  N_TRIALS)
lo_fi,  hi_fi  = wilson(res_fine.ser, N_TRIALS_FINE)

y_off,  y_on,  y_fi  = to_y(res_off.ser), to_y(res_on.ser), to_y(res_fine.ser)
ylo_off, yhi_off = to_lohi(lo_off, hi_off)
ylo_on,  yhi_on  = to_lohi(lo_on,  hi_on)
ylo_fi,  yhi_fi  = to_lohi(lo_fi,  hi_fi)

Y_LABEL_MAIN  = (r"$P_{\mathrm{sw}}$  $\equiv P(\mathrm{switch})$"
                 if METRIC == "psw" else r"SER  $\equiv P(\mathrm{fail})$")
Y_LABEL_INSET = (r"$P_{\mathrm{sw}}$" if METRIC == "psw" else r"SER")
TITLE_METRIC = "P_{\\mathrm{sw}}" if METRIC == "psw" else "SER"

# ── Single panel + inset (nested layout) ────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(8.4, 5.8))

I_uA = -I_SOT_LIST * 1e6
ax.fill_between(I_uA, ylo_off, yhi_off, color=NAVY,    alpha=0.15, lw=0)
ax.plot(I_uA, y_off, "o-", color=NAVY, lw=1.8, markersize=6,
        markerfacecolor="white", markeredgewidth=1.5, label="self-heating OFF")
ax.fill_between(I_uA, ylo_on, yhi_on,  color=CRIMSON,  alpha=0.15, lw=0)
ax.plot(I_uA, y_on,  "s--", color=CRIMSON, lw=1.8, markersize=6,
        markerfacecolor="white", markeredgewidth=1.5, label="self-heating ON")
ax.axhline(0.5, color="gray", lw=0.5, ls=":")

# Experimental reference
I_th_exp = 894e-3 / cc.R_W * 1e6
ax.axvline(I_th_exp, color=TEAL, lw=1.4, ls=(0, (4, 2)),
           label=rf"exp $I_{{\mathrm{{th}}}} = {I_th_exp:.0f}\,\mu\mathrm{{A}}$")

ax.set_xlabel(r"$|I_{\mathrm{SOT}}|$ ($\mu$A)")
ax.set_ylabel(Y_LABEL_MAIN)
ax.set_title(rf"vgsot-sim Monte-Carlo ${TITLE_METRIC}$ sweep, $t_p = {PULSE_NS:g}\,\mathrm{{ns}}$, $V_{{\mathrm{{MTJ}}}}=0$")
ax.set_ylim(-0.05, 1.08)
ax.set_xlim(0, max(I_uA) + 50)
ax.legend(loc="upper left" if METRIC == "psw" else "lower left",
          fontsize=10, framealpha=0.95)
ax.grid(True, linestyle="--", linewidth=0.4, color=THU_GRID, alpha=0.7)

# ── Inset: threshold-region sigmoid (zoom 800-1400 µA) ──────────────────
inset = ax.inset_axes([0.43, 0.10 if METRIC == "psw" else 0.42, 0.55, 0.50])
I_fi_uA = -I_SOT_INSET * 1e6
inset.fill_between(I_fi_uA, ylo_fi, yhi_fi, color=AMBER, alpha=0.20, lw=0)
inset.plot(I_fi_uA, y_fi, "o-", color=AMBER, lw=1.8, markersize=5.5,
           markerfacecolor="white", markeredgewidth=1.4)
inset.axvline(I_th_exp, color=TEAL, lw=1.2, ls=(0, (4, 2)))
inset.axhline(0.5, color="gray", lw=0.5, ls=":")
inset.set_xlabel(r"$|I_{\mathrm{SOT}}|$ ($\mu$A)", fontsize=10, labelpad=2)
inset.set_ylabel(Y_LABEL_INSET, fontsize=10, labelpad=2)
inset.tick_params(axis="both", labelsize=9)
inset.set_title("Threshold-region sigmoid (zoom)", fontsize=10.5, pad=4)
inset.set_xlim(min(I_fi_uA)-20, max(I_fi_uA)+20)
inset.set_ylim(-0.05, 1.08)
inset.grid(True, linestyle="--", linewidth=0.4, color=THU_GRID, alpha=0.6)

plt.tight_layout()
out_path = Path(__file__).resolve().parent / "fig_2_11_ser_mc.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
# Sync to vgsot-sim/article/00_chapter_drafts/figures/ (the canonical chapter copy)
chapter_fig = Path(__file__).resolve().parent.parent.parent / "article" / "00_chapter_drafts" / "figures" / "fig_2_11_ser_mc.png"
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)
print(f"Saved {out_path}")
print(f"Synced {chapter_fig}  (metric = {METRIC})")

print()
print("=" * 78)
metric_lbl = "P_sw" if METRIC == "psw" else "SER"
print(f"  Wide-range {metric_lbl}  (N = {N_TRIALS}/point)")
print("=" * 78)
for I, soff, son in zip(I_uA, res_off.ser, res_on.ser):
    if METRIC == "psw":
        soff_d, son_d = 1.0 - soff, 1.0 - son
    else:
        soff_d, son_d = soff, son
    print(f"  {I:>6.0f} µA   OFF={soff_d:.3f}   ON={son_d:.3f}   "
          f"Δ{metric_lbl}={son_d-soff_d:+.3f}")
print()
print(f"  Inset fine sweep around exp I_th={I_th_exp:.0f} µA  (N = {N_TRIALS_FINE}/point)")
for I, s in zip(I_fi_uA, res_fine.ser):
    s_d = 1.0 - s if METRIC == "psw" else s
    marker = " ← exp I_th" if abs(I - I_th_exp) < 25 else ""
    print(f"    |I_SOT| = {I:>5.0f} µA  →  {metric_lbl} = {s_d:.3f}{marker}")
print("=" * 78)
