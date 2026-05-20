"""
图 2.7 — vgsot-sim 三层架构示意图（清理版）。

清理改动相对前一版：
  * Pip Install 节点改为顶层内部子节点，避免与三层"飞线"
  * Configuration Objects → kernels 的箭头改为正交折线（先垂直再水平）
  * Updates per time step 改为单一逆向反馈箭头（无横向跨越）
  * Input data / Trace 输出仍位于左右外侧但通过顶部边线进入
  * 所有框尺寸严格对齐网格，避免越界

输出：本目录 + 同步到 ../article/00_chapter_drafts/figures/fig_07_vgsot_sim_architecture.png
"""
from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ── Palette ──────────────────────────────────────────────────────────────
THU_DEEP, THU_MID, THU_SOFT, THU_PALE, THU_TINT = "#660874", "#8B3A9E", "#A966BE", "#C99FD4", "#EFE0F7"
LAYER_TOP = "#F7ECFB"           # User-interface band (lightest)
LAYER_MID = "#E5C9EF"           # Configuration band
LAYER_LOW = "#C99FD4"           # Kernel band (deepest)
CHARCOAL, NEAR_WHITE = "#2B2B2B", "#FFFFFF"

plt.rcParams.update({
    "font.family"          : "sans-serif",
    "font.sans-serif"      : ["Arial", "Liberation Sans"],
    "font.size"            : 11,
    "mathtext.fontset"     : "stix",
    "savefig.dpi"          : 300,
    "figure.dpi"           : 150,
})

fig, ax = plt.subplots(figsize=(13.5, 8.0))
ax.set_xlim(0, 16.5)
ax.set_ylim(0, 11.0)
ax.axis("off")


# ─── Helpers ────────────────────────────────────────────────────────────
def box(x, y, w, h, label, *, fc=NEAR_WHITE, ec=CHARCOAL, fontsize=11,
        bold=True, lw=1.0, sub=None, sub_fontsize=9.5):
    bbox = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.16",
        facecolor=fc, edgecolor=ec, linewidth=lw,
    )
    ax.add_patch(bbox)
    weight = "bold" if bold else "normal"
    if sub is None:
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=fontsize, fontweight=weight, color=CHARCOAL)
    else:
        ax.text(x + w / 2, y + h * 0.66, label, ha="center", va="center",
                fontsize=fontsize, fontweight=weight, color=CHARCOAL)
        ax.text(x + w / 2, y + h * 0.30, sub, ha="center", va="center",
                fontsize=sub_fontsize, color=CHARCOAL, style="italic")


def varrow(x, y0, y1, color=CHARCOAL, lw=1.2, mut=12):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>",
                                  color=color, lw=lw, mutation_scale=mut,
                                  shrinkA=2, shrinkB=2))


def harrow(x0, x1, y, color=CHARCOAL, lw=1.2, mut=12, style="-|>"):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle=style,
                                  color=color, lw=lw, mutation_scale=mut,
                                  shrinkA=2, shrinkB=2))


def line(x0, y0, x1, y1, color=CHARCOAL, lw=1.0):
    ax.plot([x0, x1], [y0, y1], color=color, lw=lw, solid_capstyle="round")


# ─── Three coloured layer bands ─────────────────────────────────────────
LAYER_X, LAYER_W = 1.05, 14.50
TOP_Y,   TOP_H   = 8.65, 1.85
MID_Y,   MID_H   = 5.55, 2.85
LOW_Y,   LOW_H   = 1.50, 3.80

for y, h, c in [(TOP_Y, TOP_H, LAYER_TOP),
                (MID_Y, MID_H, LAYER_MID),
                (LOW_Y, LOW_H, LAYER_LOW)]:
    ax.add_patch(patches.FancyBboxPatch(
        (LAYER_X, y), LAYER_W, h,
        boxstyle="round,pad=0.05,rounding_size=0.22",
        facecolor=c, edgecolor="none"))

# Layer titles (top-left of each band, italic accent on the right)
ax.text(LAYER_X + 0.30, TOP_Y + TOP_H - 0.18,
        "User Interfaces & IO",
        ha="left", va="top", fontsize=12.5, fontweight="bold", color=CHARCOAL)
ax.text(LAYER_X + LAYER_W - 0.30, TOP_Y + TOP_H - 0.18,
        "CLI · Python API",
        ha="right", va="top", fontsize=10, style="italic", color=CHARCOAL)

ax.text(LAYER_X + 0.30, MID_Y + MID_H - 0.18,
        "Experiment Scenarios & Configurations",
        ha="left", va="top", fontsize=12.5, fontweight="bold", color=CHARCOAL)
ax.text(LAYER_X + LAYER_W - 0.30, MID_Y + MID_H - 0.18,
        "cases · configs",
        ha="right", va="top", fontsize=10, style="italic", color=CHARCOAL)

ax.text(LAYER_X + 0.30, LOW_Y + LOW_H - 0.18,
        "Physical Kernels & Core Modules",
        ha="left", va="top", fontsize=12.5, fontweight="bold", color=CHARCOAL)
ax.text(LAYER_X + LAYER_W - 0.30, LOW_Y + LOW_H - 0.18,
        "single-step update loop",
        ha="right", va="top", fontsize=10, style="italic", color=CHARCOAL)


# ════════════════════════════════════════════════════════════════════════
# TOP LAYER — CLI, Python API, Pip Install (all inside the band)
# ════════════════════════════════════════════════════════════════════════
TOP_BOX_Y, TOP_BOX_H = TOP_Y + 0.20, 0.95
box(1.40, TOP_BOX_Y, 3.30, TOP_BOX_H, "Pip Install",
    sub="pip install vgsot-sim",
    fc=NEAR_WHITE, ec=CHARCOAL, fontsize=11.5, sub_fontsize=9.5)
box(5.10, TOP_BOX_Y, 4.80, TOP_BOX_H, "Command-Line Interface",
    fc=NEAR_WHITE, ec=CHARCOAL, fontsize=11.5)
box(10.30, TOP_BOX_Y, 4.85, TOP_BOX_H, "Python API",
    fc=NEAR_WHITE, ec=CHARCOAL, fontsize=11.5)


# ════════════════════════════════════════════════════════════════════════
# MIDDLE LAYER — Standardized cases (left) + Configuration Objects (right)
# ════════════════════════════════════════════════════════════════════════
# Standardized cases outer frame holding 4 sub-cases (2x2 grid)
CASES_X, CASES_Y, CASES_W, CASES_H = 1.40, 5.75, 8.30, 1.90
ax.add_patch(FancyBboxPatch(
    (CASES_X, CASES_Y), CASES_W, CASES_H,
    boxstyle="round,pad=0.04,rounding_size=0.16",
    facecolor=NEAR_WHITE, edgecolor=CHARCOAL, lw=1.0))
# Sub-cases 2x2 — fill the cases container without a redundant inner title
sub_w, sub_h = 3.85, 0.70
sub_y_top    = CASES_Y + CASES_H - sub_h - 0.20
sub_y_bot    = CASES_Y + 0.20
sub_left_x   = CASES_X + 0.25
sub_right_x  = CASES_X + 0.25 + sub_w + 0.20
for (sx, sy, label) in [
    (sub_left_x,  sub_y_top, "SOT baseline case"),
    (sub_right_x, sub_y_top, "VCMA-assisted case"),
    (sub_left_x,  sub_y_bot, "Optimized VGSOT case"),
    (sub_right_x, sub_y_bot, "SER Monte Carlo case"),
]:
    box(sx, sy, sub_w, sub_h, label,
        fc=NEAR_WHITE, ec=CHARCOAL, fontsize=10.5, bold=False, lw=0.8)

# Configuration Objects (right)
CFG_X, CFG_Y, CFG_W, CFG_H = 10.30, 6.05, 4.85, 1.30
box(CFG_X, CFG_Y, CFG_W, CFG_H,
    "Configuration Objects",
    sub=r"$(J,\ V,\ t_w,\ \mathrm{etc.})$",
    fc=NEAR_WHITE, ec=CHARCOAL, fontsize=11.5, sub_fontsize=10.5)

# cases → configuration arrow (horizontal, inside middle band)
harrow(CASES_X + CASES_W + 0.05, CFG_X - 0.05,
       (CFG_Y + CFG_H / 2),
       color=CHARCOAL, lw=1.2, mut=12)


# ════════════════════════════════════════════════════════════════════════
# BOTTOM LAYER — 5 kernel modules in a row, all inside band
# ════════════════════════════════════════════════════════════════════════
K_Y, K_H = 2.00, 1.85

# dynamic_switching (largest, on the left with LLG formula)
DSx, DSw = 1.40, 4.30
ax.add_patch(FancyBboxPatch(
    (DSx, K_Y), DSw, K_H,
    boxstyle="round,pad=0.04,rounding_size=0.16",
    facecolor=NEAR_WHITE, edgecolor=CHARCOAL, lw=1.0))
ax.text(DSx + DSw / 2, K_Y + K_H * 0.78, "dynamic_switching",
        ha="center", va="center", fontsize=12, fontweight="bold", color=CHARCOAL)
ax.text(DSx + DSw / 2, K_Y + K_H * 0.52, "(LLG Integration)",
        ha="center", va="center", fontsize=10, color=CHARCOAL, style="italic")
ax.text(DSx + DSw / 2, K_Y + K_H * 0.22,
        r"$\mathbf{m}_{n+1} = \mathcal{F}(\mathbf{m}_n,\,\mathbf{H}_{\rm eff},\,\mathbf{H}_{\rm th})$",
        ha="center", va="center", fontsize=11, color=CHARCOAL)

# Four narrower modules (anisotropy / stochastic / electronic / tmr)
MOD_W = 2.25
MOD_GAP = 0.10
MOD_X0 = DSx + DSw + 0.20
MOD_Y = K_Y + 0.20
MOD_H = K_H - 0.40

def mod_box(idx, title, sub):
    x = MOD_X0 + idx * (MOD_W + MOD_GAP)
    ax.add_patch(FancyBboxPatch(
        (x, MOD_Y), MOD_W, MOD_H,
        boxstyle="round,pad=0.04,rounding_size=0.14",
        facecolor=NEAR_WHITE, edgecolor=CHARCOAL, lw=1.0))
    ax.text(x + MOD_W / 2, MOD_Y + MOD_H * 0.66, title,
            ha="center", va="center", fontsize=11, fontweight="bold", color=CHARCOAL)
    ax.text(x + MOD_W / 2, MOD_Y + MOD_H * 0.30, sub,
            ha="center", va="center", fontsize=9.5, color=CHARCOAL, style="italic")
    return x

ani_x  = mod_box(0, "anisotropy",  "(Effective Field)")
sto_x  = mod_box(1, "stochastic",  "(Thermal Noise)")
ele_x  = mod_box(2, "electronic",  "(Transport)")
tmr_x  = mod_box(3, "tmr",         "(Resistance)")


# ════════════════════════════════════════════════════════════════════════
# Inter-layer arrows (clean, top-to-bottom, no diagonal flying)
# ════════════════════════════════════════════════════════════════════════
# Top → Middle: CLI/API both feed the cases (one arrow per top box, vertical)
varrow(5.10 + 4.80 / 2, TOP_BOX_Y, MID_Y + MID_H,  # CLI → middle
       color=THU_DEEP, lw=1.4, mut=12)
varrow(10.30 + 4.85 / 2, TOP_BOX_Y, CFG_Y + CFG_H,  # Python API → cfg
       color=THU_DEEP, lw=1.4, mut=12)
varrow(1.40 + 3.30 / 2, TOP_BOX_Y, MID_Y + MID_H,  # Pip → middle
       color=CHARCOAL, lw=1.0, mut=10)

# Middle → Bottom: Configuration Objects fans down into each kernel via
# an orthogonal route (drop straight down to a "bus" line, then each
# kernel takes a short branch up).  No diagonals.
BUS_Y = LOW_Y + LOW_H + 0.20   # just above bottom band
CFG_DROP_X = CFG_X + CFG_W / 2
# Vertical drop from cfg to bus
line(CFG_DROP_X, CFG_Y, CFG_DROP_X, BUS_Y,
     color=CHARCOAL, lw=1.0)
# Horizontal bus
bus_left  = DSx + DSw / 2
bus_right = CFG_DROP_X
line(bus_left, BUS_Y, bus_right, BUS_Y, color=CHARCOAL, lw=1.0)
# Vertical branches from bus into each kernel top edge
for cx in [DSx + DSw / 2,
           ani_x + MOD_W / 2,
           sto_x + MOD_W / 2,
           ele_x + MOD_W / 2,
           tmr_x + MOD_W / 2]:
    varrow(cx, BUS_Y, K_Y + K_H + 0.04,
           color=CHARCOAL, lw=0.9, mut=10)

# Update loop: single curved feedback from tmr back to dynamic_switching
# (drawn as two right-angle segments along the bottom band's lower edge)
FEEDBACK_Y = LOW_Y + 0.20
line(tmr_x + MOD_W / 2, K_Y - 0.05, tmr_x + MOD_W / 2, FEEDBACK_Y,
     color=THU_DEEP, lw=1.2)
line(tmr_x + MOD_W / 2, FEEDBACK_Y, DSx + DSw / 2, FEEDBACK_Y,
     color=THU_DEEP, lw=1.2)
varrow(DSx + DSw / 2, FEEDBACK_Y, K_Y - 0.05,
       color=THU_DEEP, lw=1.2, mut=12)
ax.text((DSx + DSw / 2 + tmr_x + MOD_W / 2) / 2,
        FEEDBACK_Y - 0.20, "update per time step",
        ha="center", va="top", fontsize=10, style="italic", color=THU_DEEP)


# ════════════════════════════════════════════════════════════════════════
# External IO arrows — labelled at top corners with short stems
# ════════════════════════════════════════════════════════════════════════
# Input data — far-left, side-entry into dynamic_switching (bottom layer)
ax.text(0.20, K_Y + K_H * 0.50, "Input\ndata", ha="left", va="center",
        fontsize=11, fontweight="bold", color=CHARCOAL,
        multialignment="left")
harrow(0.95, DSx - 0.05, K_Y + K_H * 0.50,
       color=CHARCOAL, lw=1.2, mut=12)

# Output — far-right, leaves tmr → out to label
ax.text(16.45, K_Y + K_H * 0.50,
        "Trace\nProbability\nSER stats",
        ha="right", va="center", fontsize=10, fontweight="bold", color=CHARCOAL,
        multialignment="left", linespacing=1.15)
harrow(tmr_x + MOD_W + 0.04, 15.60, K_Y + K_H * 0.50,
       color=CHARCOAL, lw=1.2, mut=12)


# ════════════════════════════════════════════════════════════════════════
# Save + sync
# ════════════════════════════════════════════════════════════════════════
out_path = Path(__file__).resolve().parent / "fig_07_vgsot_sim_architecture.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight",
            facecolor="white", pad_inches=0.20)
plt.close()
chapter_fig = (Path(__file__).resolve().parent.parent /
               "article" / "00_chapter_drafts" / "figures" / "fig_07_vgsot_sim_architecture.png")
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)
print(f"Saved  {out_path}")
print(f"Synced {chapter_fig}")
