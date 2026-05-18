"""
图 2.2 — 面向概率计算的sMTJ行为级模型分层架构。

自下而上四个层次：
  Device layer：PMA sMTJ + 三类驱动（V_MTJ, I_SOT, thermal fluctuations）
  Physical-model layer：LLG 方程、调制能垒 ΔE_b、Néel-Brown 模型 → P_sw
  Behavioral layer：Sigmoid 近似 P_sw(u) ≈ σ(β_s (u - u_th))
  Computation layer：(t_w, I_SOT, V_MTJ) → u → Bernoulli 采样 → 随机比特流

清华紫配色，紧凑布局，向上的抽象提升箭头连接四层。

输出：fig_2_1_2_behavioral_layers.png（本目录 + 同步到 ../article/00_chapter_drafts/figures/）
"""
from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

# ── Palette (THU purple family) ─────────────────────────────────────────
THU_DEEP, THU_MID, THU_SOFT, THU_PALE = "#660874", "#8B3A9E", "#A966BE", "#C99FD4"
THU_TINT, THU_LIGHT = "#EFE0F7", "#F7ECFB"
THU_GRID = "#DDD0E8"
CRIMSON, NAVY, TEAL, AMBER = "#A82038", "#1F5FA8", "#1A6B5A", "#C47A00"
CHARCOAL, NEAR_WHITE = "#2B2B2B", "#FFFFFF"

# Layer-band colours (light to deeper bottom-up — gradient stays in THU family)
BAND_TOP    = "#F7ECFB"   # Computation
BAND_BEHAV  = "#E5C9EF"   # Behavioral
BAND_PHYS   = "#C99FD4"   # Physical-model
BAND_DEVICE = "#A966BE"   # Device

plt.rcParams.update({
    "font.family"      : "sans-serif",
    "font.sans-serif"  : ["Arial", "Liberation Sans"],
    "font.size"        : 11,
    "mathtext.fontset" : "stix",
    "savefig.dpi"      : 300,
    "figure.dpi"       : 150,
})

# Compact landscape layout
fig, ax = plt.subplots(figsize=(12.5, 8.0))
ax.set_xlim(0, 16.0)
ax.set_ylim(0, 10.5)
ax.axis("off")


# ─── Helpers ────────────────────────────────────────────────────────────
def band(y, h, color, title, *, title_color=CHARCOAL, title_fontsize=12.5):
    ax.add_patch(FancyBboxPatch(
        (LAYER_X, y), LAYER_W, h,
        boxstyle="round,pad=0.05,rounding_size=0.20",
        facecolor=color, edgecolor="none"))
    ax.text(LAYER_X + 0.30, y + h - 0.20, title,
            ha="left", va="top", fontsize=title_fontsize,
            fontweight="bold", color=title_color)


def box(x, y, w, h, label, *, fc=NEAR_WHITE, ec=CHARCOAL, fontsize=10.5,
        bold=True, lw=1.0, sub=None, sub_fontsize=9.0):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.14",
        facecolor=fc, edgecolor=ec, linewidth=lw))
    weight = "bold" if bold else "normal"
    if sub is None:
        ax.text(x + w / 2, y + h / 2, label,
                ha="center", va="center", fontsize=fontsize, fontweight=weight,
                color=CHARCOAL)
    else:
        ax.text(x + w / 2, y + h * 0.72, label,
                ha="center", va="center", fontsize=fontsize, fontweight=weight,
                color=CHARCOAL)
        ax.text(x + w / 2, y + h * 0.28, sub,
                ha="center", va="center", fontsize=sub_fontsize,
                color=CHARCOAL, style="italic")


def arrow(x0, y0, x1, y1, *, color=CHARCOAL, lw=1.4, style="-|>", mut=14):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                        color=color, lw=lw, mutation_scale=mut,
                        shrinkA=2, shrinkB=2)
    ax.add_patch(a)


# ─── Layer geometry — tighter rows, wider boxes ─────────────────────────
LAYER_X, LAYER_W = 0.4, 15.2
DEVICE_Y, DEVICE_H = 0.5, 2.10
PHYS_Y,   PHYS_H   = 3.0, 2.10
BEHAV_Y,  BEHAV_H  = 5.5, 2.10
COMP_Y,   COMP_H   = 8.0, 2.20

# Bands (bottom-up, header text in top-left of each band)
band(DEVICE_Y, DEVICE_H, BAND_DEVICE, "Device layer",
     title_color=NEAR_WHITE, title_fontsize=12.5)
band(PHYS_Y, PHYS_H, BAND_PHYS, "Physical-model layer",
     title_color=NEAR_WHITE, title_fontsize=12.5)
band(BEHAV_Y, BEHAV_H, BAND_BEHAV, "Behavioral layer",
     title_color=CHARCOAL, title_fontsize=12.5)
band(COMP_Y, COMP_H, BAND_TOP, "Computation layer",
     title_color=CHARCOAL, title_fontsize=12.5)

# Layer-tagline at top-right of each band (deeper info, italic)
def tagline(y, h, text, color):
    ax.text(LAYER_X + LAYER_W - 0.30, y + h - 0.20, text,
            ha="right", va="top", fontsize=10.5,
            style="italic", color=color)
tagline(DEVICE_Y, DEVICE_H, "PMA sMTJ + drives", NEAR_WHITE)
tagline(PHYS_Y,   PHYS_H,   "LLG + Néel–Brown", NEAR_WHITE)
tagline(BEHAV_Y,  BEHAV_H,  "Sigmoid approximation", CHARCOAL)
tagline(COMP_Y,   COMP_H,   "programmable Bernoulli source", CHARCOAL)


# ════════════════════════════════════════════════════════════════════════
# Device layer — sMTJ pillar (left) + three drive boxes (right)
# ════════════════════════════════════════════════════════════════════════
# sMTJ pillar schematic (left side)
SCH_X0, SCH_X1 = 0.9, 3.7
HM_Y0, HM_Y1   = DEVICE_Y + 0.20, DEVICE_Y + 0.45
ax.add_patch(Rectangle((SCH_X0, HM_Y0), SCH_X1 - SCH_X0, HM_Y1 - HM_Y0,
                       facecolor="#6E6E6E", edgecolor=CHARCOAL, lw=0.9))
ax.text((SCH_X0 + SCH_X1) / 2.0, (HM_Y0 + HM_Y1) / 2.0, "HM channel",
        ha="center", va="center", fontsize=8.5, color=NEAR_WHITE)
MTJ_CX = (SCH_X0 + SCH_X1) / 2.0
MTJ_HW = 0.48
RL_Y0, RL_Y1 = HM_Y1, HM_Y1 + 0.28
ax.add_patch(Rectangle((MTJ_CX - MTJ_HW, RL_Y0), 2 * MTJ_HW, RL_Y1 - RL_Y0,
                       facecolor="#3F5F8C", edgecolor=CHARCOAL, lw=0.9))
ax.annotate("", xy=(MTJ_CX - 0.15, RL_Y1 - 0.04),
            xytext=(MTJ_CX - 0.15, RL_Y0 + 0.04),
            arrowprops=dict(arrowstyle="-|>", color=NEAR_WHITE, lw=1.1))
MGO_Y0, MGO_Y1 = RL_Y1, RL_Y1 + 0.16
ax.add_patch(Rectangle((MTJ_CX - MTJ_HW, MGO_Y0), 2 * MTJ_HW, MGO_Y1 - MGO_Y0,
                       facecolor="#D9C28A", edgecolor=CHARCOAL, lw=0.9))
FL_Y0, FL_Y1 = MGO_Y1, MGO_Y1 + 0.28
ax.add_patch(Rectangle((MTJ_CX - MTJ_HW, FL_Y0), 2 * MTJ_HW, FL_Y1 - FL_Y0,
                       facecolor="#3F5F8C", edgecolor=CHARCOAL, lw=0.9))
ax.annotate("", xy=(MTJ_CX - 0.15, FL_Y1 - 0.04),
            xytext=(MTJ_CX - 0.15, FL_Y0 + 0.04),
            arrowprops=dict(arrowstyle="<|-|>", color=NEAR_WHITE, lw=1.1))
TOP_Y0, TOP_Y1 = FL_Y1, FL_Y1 + 0.12
ax.add_patch(Rectangle((MTJ_CX - MTJ_HW, TOP_Y0), 2 * MTJ_HW, TOP_Y1 - TOP_Y0,
                       facecolor=CHARCOAL))

# Three drive boxes (right of pillar)
DRV_Y = DEVICE_Y + 0.40
DRV_H = 1.20
DRV_W = 3.55
box(4.4,  DRV_Y, DRV_W, DRV_H, "VCMA bias",
    sub=r"$V_{\rm MTJ}$ on MgO barrier",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)
box(8.2,  DRV_Y, DRV_W, DRV_H, "SOT current",
    sub=r"$I_{\rm SOT}$ through HM channel",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)
box(12.0, DRV_Y, DRV_W, DRV_H, "Thermal noise",
    sub=r"$\mathbf{H}_{\rm th}(t)$ (Brown 1963)",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)


# ════════════════════════════════════════════════════════════════════════
# Physical-model layer
# ════════════════════════════════════════════════════════════════════════
PHY_Y0 = PHYS_Y + 0.40
PHY_H = 1.20

box(0.9, PHY_Y0, 4.8, PHY_H, "Stochastic LLG",
    sub=r"$d\mathbf{m}/dt$ with SOT, Gilbert, $\mathbf{H}_{\rm th}$",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)
box(6.0, PHY_Y0, 4.5, PHY_H, "Modulated barrier",
    sub=r"$\Delta E_b(V_{\rm MTJ}, I_{\rm SOT})$  via VCMA + SHE",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)
box(10.8, PHY_Y0, 4.4, PHY_H, r"Néel–Brown $\to P_{\rm sw}$",
    sub=r"$P_{\rm sw}=1-e^{-t_w/\tau_0\cdot e^{-\Delta_{\rm eff}}}$",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)


# ════════════════════════════════════════════════════════════════════════
# Behavioral layer — sigmoid sketch (left) + parameter boxes (right)
# ════════════════════════════════════════════════════════════════════════
# Sigmoid sketch as inline plot
SIG_X0, SIG_Y0, SIG_W, SIG_H = 0.9, BEHAV_Y + 0.35, 3.8, 1.30
ax.add_patch(FancyBboxPatch(
    (SIG_X0, SIG_Y0), SIG_W, SIG_H,
    boxstyle="round,pad=0.04,rounding_size=0.14",
    facecolor=NEAR_WHITE, edgecolor=CHARCOAL, lw=1.0))
# Mini sigmoid curve
mini_x = np.linspace(-4, 4, 200)
mini_y = 1.0 / (1.0 + np.exp(-mini_x))
mx0, my0 = SIG_X0 + 0.45, SIG_Y0 + 0.18
mw, mh   = SIG_W - 0.60, SIG_H - 0.55
ax.plot(mx0 + (mini_x + 4) / 8.0 * mw, my0 + mini_y * mh,
        color=THU_DEEP, lw=2.0)
ax.plot([mx0, mx0 + mw], [my0, my0], color=CHARCOAL, lw=0.9)
ax.plot([mx0, mx0], [my0, my0 + mh], color=CHARCOAL, lw=0.9)
ax.text(mx0 + mw / 2.0, my0 - 0.13, r"$u - u_{\rm th}$",
        ha="center", va="top", fontsize=9.5)
ax.text(mx0 - 0.10, my0 + mh / 2.0, r"$P_{\rm sw}$",
        ha="right", va="center", fontsize=9.5)
ax.text(SIG_X0 + SIG_W / 2.0, SIG_Y0 + SIG_H - 0.18,
        r"$P_{\rm sw}(u)=\sigma(\beta_s(u-u_{\rm th}))$",
        ha="center", va="top", fontsize=11, fontweight="bold", color=CHARCOAL)

box(5.0, BEHAV_Y + 0.35, 4.8, 1.30, "Behavioral parameters",
    sub=r"$\beta_s\!=\!2\kappa\ln 2/(k_BT)$;  $u_{\rm th}\!\downarrow$ with $t_w,T$",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)

box(10.1, BEHAV_Y + 0.35, 5.1, 1.30, "Physics-to-parameter map",
    sub=r"$\beta_s = f(\Delta,\xi,\theta_{\rm SH})$;  D2D $\to\beta_s\!\downarrow$",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)


# ════════════════════════════════════════════════════════════════════════
# Computation layer
# ════════════════════════════════════════════════════════════════════════
CMP_Y0 = COMP_Y + 0.50
CMP_H = 1.30

box(0.9,  CMP_Y0, 3.5, CMP_H, "External drives",
    sub=r"$(t_w,\ I_{\rm SOT},\ V_{\rm MTJ})$",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10.5)
box(4.7,  CMP_Y0, 3.0, CMP_H, r"$u$-mapping",
    sub="effective drive",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)
box(8.0,  CMP_Y0, 3.5, CMP_H, "Bernoulli sampler",
    sub=r"$m\sim\mathrm{Bern}(P_{\rm sw}(u))$",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)
box(11.8, CMP_Y0, 3.4, CMP_H, "Random bit stream",
    sub="stochastic source",
    fc=NEAR_WHITE, fontsize=11.5, sub_fontsize=10)

# Horizontal flow arrows
for x0, x1 in [(4.4, 4.7), (7.7, 8.0), (11.5, 11.8)]:
    arrow(x0, CMP_Y0 + CMP_H / 2, x1, CMP_Y0 + CMP_H / 2,
          color=THU_DEEP, lw=1.4, mut=12)


# ════════════════════════════════════════════════════════════════════════
# Inter-layer upward arrows (compact — 3 columns to keep visual order)
# ════════════════════════════════════════════════════════════════════════
def vertical_arrow(x, y0, y1, color=THU_DEEP, lw=1.4):
    arrow(x, y0, x, y1, color=color, lw=lw, mut=12)

# Device -> Physics
for x in (3.0, 8.0, 13.0):
    vertical_arrow(x, DEVICE_Y + DEVICE_H, PHYS_Y)
# Physics -> Behavioral
for x in (3.0, 8.0, 13.0):
    vertical_arrow(x, PHYS_Y + PHYS_H, BEHAV_Y)
# Behavioral -> Computation
for x in (3.0, 8.0, 13.0):
    vertical_arrow(x, BEHAV_Y + BEHAV_H, COMP_Y)


# ════════════════════════════════════════════════════════════════════════
# Save + sync
# ════════════════════════════════════════════════════════════════════════
out_path = Path(__file__).resolve().parent / "fig_2_1_2_behavioral_layers.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white",
            pad_inches=0.15)
plt.close()
chapter_fig = (Path(__file__).resolve().parent.parent /
               "article" / "00_chapter_drafts" / "figures" / "fig_2_1_2_behavioral_layers.png")
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)
print(f"Saved  {out_path}")
print(f"Synced {chapter_fig}")
