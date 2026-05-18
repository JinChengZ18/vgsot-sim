"""
图 2.3 — sMTJ温度效应与非理想性建模总览。

四条非理想通道（左侧），共同汇聚到右侧 P_sw 响应：
  (i)   自热源 → 热扩散 → T(t)
  (ii)  T-依赖材料参数：M_s(T), K_i(T), η(T)
  (iii) TMR 与电输运非线性：R_MTJ(V, T, θ)
  (iv)  退磁与形状效应：N_x, N_y, N_z

中间节点表示 LLG 演化 (m(t))。所有通道汇至 P_sw，
P_sw 响应有两个观察量：阈值漂移 u_th(T) 与斜率展宽 β_s(T, D2D)。

输出：fig_2_2_1_thermal_nonidealities.png（本目录 + 同步到 ../article/00_chapter_drafts/figures/）
"""
from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# Palette
THU_DEEP, THU_MID, THU_SOFT, THU_PALE = "#660874", "#8B3A9E", "#A966BE", "#C99FD4"
THU_TINT, THU_GRID = "#EFE0F7", "#DDD0E8"
CRIMSON, NAVY, TEAL, AMBER = "#A82038", "#1F5FA8", "#1A6B5A", "#C47A00"
CHARCOAL, NEAR_WHITE = "#2B2B2B", "#FFFFFF"

# Per-channel accent colours
SH_COLOR = "#C0392B"   # self-heating
TM_COLOR = "#1F5FA8"   # T-dep materials
TR_COLOR = "#1A6B5A"   # TMR / transport
DM_COLOR = "#C47A00"   # demag / shape

plt.rcParams.update({
    "font.family"      : "sans-serif",
    "font.sans-serif"  : ["Arial", "Liberation Sans"],
    "font.size"        : 11,
    "mathtext.fontset" : "stix",
    "savefig.dpi"      : 300,
    "figure.dpi"       : 150,
})

# More compact aspect, less empty space
fig, ax = plt.subplots(figsize=(12.5, 6.4))
ax.set_xlim(0, 14.5)
ax.set_ylim(0, 8.5)
ax.axis("off")


def box(x, y, w, h, label, *, fc=NEAR_WHITE, ec=CHARCOAL, lw=1.2,
        fontsize=11.5, sub=None, sub_fontsize=10):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.16",
        facecolor=fc, edgecolor=ec, linewidth=lw))
    if sub is None:
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", color=CHARCOAL)
    else:
        ax.text(x + w / 2, y + h * 0.66, label, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", color=CHARCOAL)
        ax.text(x + w / 2, y + h * 0.30, sub, ha="center", va="center",
                fontsize=sub_fontsize, color=CHARCOAL, style="italic")


def arrow(x0, y0, x1, y1, *, color=CHARCOAL, lw=1.6, style="-|>", mut=16):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                        color=color, lw=lw, mutation_scale=mut,
                        shrinkA=3, shrinkB=3)
    ax.add_patch(a)


# ════════════════════════════════════════════════════════════════════════
# Title (top, full-width)
# ════════════════════════════════════════════════════════════════════════
ax.text(7.25, 8.20, "Thermal & non-ideal effects in sMTJ behavioural modelling",
        ha="center", va="top", fontsize=13.5, fontweight="bold", color=CHARCOAL)


# ════════════════════════════════════════════════════════════════════════
# Left column — four non-ideality channels (tighter rows)
# ════════════════════════════════════════════════════════════════════════
CH_X, CH_W = 0.5, 4.20
CH_H = 1.30
CH_Y_LIST = [6.30, 4.65, 3.00, 1.35]
CH_LABELS = [
    ("Self-heating",
     r"$C_v\,dT/dt=\lambda\nabla^2T+Q\to T(t)$",
     SH_COLOR),
    ("T-dependent materials",
     r"$M_s(T),\ K_i(T),\ \eta(T)$",
     TM_COLOR),
    ("TMR & transport",
     r"$R_{\rm MTJ}(V,T,\theta)$",
     TR_COLOR),
    ("Demag & shape",
     r"$\mathbf{H}_{\rm D}=-N_{\!xyz}\,\mathbf{M}$",
     DM_COLOR),
]
for (y, (title, sub, col)) in zip(CH_Y_LIST, CH_LABELS):
    box(CH_X, y, CH_W, CH_H, title, sub=sub,
        fc=NEAR_WHITE, ec=col, lw=1.7,
        fontsize=12, sub_fontsize=10.5)
    # accent stripe
    ax.add_patch(plt.Rectangle((CH_X, y), 0.10, CH_H, facecolor=col,
                               edgecolor="none"))

# Channel index labels (i)–(iv)
for i, y in enumerate(CH_Y_LIST):
    label = ["(i)", "(ii)", "(iii)", "(iv)"][i]
    ax.text(CH_X - 0.30, y + CH_H / 2, label,
            ha="right", va="center", fontsize=12,
            fontweight="bold", color=CHARCOAL)


# ════════════════════════════════════════════════════════════════════════
# Middle column — LLG dynamic core (taller, more central)
# ════════════════════════════════════════════════════════════════════════
MID_X, MID_Y, MID_W, MID_H = 5.55, 3.25, 3.80, 2.80
ax.add_patch(FancyBboxPatch(
    (MID_X, MID_Y), MID_W, MID_H,
    boxstyle="round,pad=0.06,rounding_size=0.22",
    facecolor=THU_TINT, edgecolor=THU_DEEP, linewidth=1.8))
ax.text(MID_X + MID_W / 2, MID_Y + MID_H - 0.35,
        "Magnetisation dynamics",
        ha="center", va="top", fontsize=13, fontweight="bold", color=THU_DEEP)
ax.text(MID_X + MID_W / 2, MID_Y + MID_H * 0.50,
        r"$\dfrac{d\mathbf{m}}{dt}=-\gamma\mathbf{m}\!\times\!\mathbf{H}_{\rm eff}"
        r"\,+\,\alpha\mathbf{m}\!\times\!\dfrac{d\mathbf{m}}{dt}\,+\,"
        r"\boldsymbol{\tau}_{\rm SOT}\,+\,\boldsymbol{\tau}_{\rm th}$",
        ha="center", va="center", fontsize=11, color=CHARCOAL)
ax.text(MID_X + MID_W / 2, MID_Y + MID_H * 0.13,
        r"updated each $\Delta t$ with $T$-corrected parameters",
        ha="center", va="center", fontsize=10, color=CHARCOAL, style="italic")


# ════════════════════════════════════════════════════════════════════════
# Right column — P_sw response with two observables
# ════════════════════════════════════════════════════════════════════════
PSW_X, PSW_Y, PSW_W, PSW_H = 9.80, 2.80, 4.40, 3.70
ax.add_patch(FancyBboxPatch(
    (PSW_X, PSW_Y), PSW_W, PSW_H,
    boxstyle="round,pad=0.06,rounding_size=0.22",
    facecolor="#FFF8DC", edgecolor=CRIMSON, linewidth=1.8))
ax.text(PSW_X + PSW_W / 2, PSW_Y + PSW_H - 0.30,
        r"$P_{\rm sw}(u)$ response",
        ha="center", va="top", fontsize=13, fontweight="bold", color=CRIMSON)

# Mini sigmoid plot inside (larger, more usable space)
sig_x = np.linspace(-4, 4, 200)
sig_y_nom = 1.0 / (1.0 + np.exp(-1.4 * sig_x))
sig_y_drift = 1.0 / (1.0 + np.exp(-0.9 * (sig_x + 0.8)))

mx0, my0 = PSW_X + 0.55, PSW_Y + 0.55
mw, mh   = PSW_W - 0.90, PSW_H - 1.30
ax.plot(mx0 + (sig_x + 4) / 8.0 * mw, my0 + sig_y_nom * mh,
        color=CHARCOAL, lw=1.8)
ax.plot(mx0 + (sig_x + 4) / 8.0 * mw, my0 + sig_y_drift * mh,
        color=CRIMSON, lw=1.8, ls="--")
ax.plot([mx0, mx0 + mw], [my0, my0], color=CHARCOAL, lw=0.9)
ax.plot([mx0, mx0], [my0, my0 + mh], color=CHARCOAL, lw=0.9)
ax.text(mx0 + mw / 2, my0 - 0.20, r"$u$",
        ha="center", va="top", fontsize=11)
ax.text(mx0 - 0.10, my0 + mh / 2, r"$P_{\rm sw}$",
        ha="right", va="center", fontsize=11)

# Threshold shift annotation
ax.annotate("",
            xy=(mx0 + 0.30 * mw, my0 + 0.50 * mh),
            xytext=(mx0 + 0.55 * mw, my0 + 0.50 * mh),
            arrowprops=dict(arrowstyle="<|-", color=CRIMSON, lw=1.3))
ax.text(mx0 + 0.07 * mw, my0 + 0.40 * mh, "threshold shift",
        ha="left", va="top", fontsize=9, color=CRIMSON, style="italic")

ax.annotate("",
            xy=(mx0 + 0.82 * mw, my0 + 0.30 * mh),
            xytext=(mx0 + 0.82 * mw, my0 + 0.85 * mh),
            arrowprops=dict(arrowstyle="<->", color=CRIMSON, lw=1.3))
ax.text(mx0 + 0.84 * mw, my0 + 0.55 * mh, "slope\nbroadening",
        ha="left", va="center", fontsize=9, color=CRIMSON, style="italic")


# ════════════════════════════════════════════════════════════════════════
# Arrows: each channel → middle LLG core (colour-matched)
# ════════════════════════════════════════════════════════════════════════
def channel_to_core(idx, col):
    y_c = CH_Y_LIST[idx] + CH_H / 2.0
    y_mid_attach = MID_Y + MID_H * (0.86 - 0.22 * idx)
    arrow(CH_X + CH_W + 0.05, y_c,
          MID_X - 0.05, y_mid_attach,
          color=col, lw=1.6, mut=16)

for i, (_, _, c) in enumerate(CH_LABELS):
    channel_to_core(i, c)

# LLG → P_sw arrow
arrow(MID_X + MID_W + 0.05, MID_Y + MID_H * 0.50,
      PSW_X - 0.05, PSW_Y + PSW_H * 0.55,
      color=THU_DEEP, lw=2.0, mut=18)
ax.text((MID_X + MID_W + PSW_X) / 2.0, MID_Y + MID_H * 0.50 + 0.35,
        r"sampling + $t_w$ window",
        ha="center", va="bottom", fontsize=10, color=THU_DEEP, style="italic")


# ════════════════════════════════════════════════════════════════════════
# Save + sync
# ════════════════════════════════════════════════════════════════════════
out_path = Path(__file__).resolve().parent / "fig_2_2_1_thermal_nonidealities.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white",
            pad_inches=0.15)
plt.close()
chapter_fig = (Path(__file__).resolve().parent.parent /
               "article" / "00_chapter_drafts" / "figures" / "fig_2_2_1_thermal_nonidealities.png")
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)
print(f"Saved  {out_path}")
print(f"Synced {chapter_fig}")
