"""
图 2.1.1 — 三端SOT-sMTJ器件结构与T型等效电路网络。

左侧子图：物理堆叠示意。HM沟道（重金属）为底层水平条带；其上为MTJ柱
（参考层 / MgO势垒 / 自由层），自由层与参考层均以PMA箭头标出。MTJ顶电极
为 T1，HM沟道两端为 T2、T3。写电流 I_SOT 经 T2→T3 沿沟道流过产生
横向自旋积累 (σ̂)；读路径仅在 T1 与沟道接地端之间施加小偏压。

右侧子图：T型等效电阻网络。R_MTJ 为隧穿电阻（用可变电阻符号表示），
HM沟道沿长度方向均分为两段 R_SOT/2，三段电阻在中间节点 N 处汇合。
读出操作在 T1—N 之间获取 V_MTJ；写入操作在 T2—T3 之间施加 V_2 − V_3
驱动 I_SOT = (V_2 − V_3) / R_SOT。

输出：fig_01_t_circuit.png（本目录 + 同步到 ../article/00_chapter_drafts/figures/）
"""
from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch

# ── Palette ─────────────────────────────────────────────────────────────
THU_DEEP, THU_MID, THU_SOFT, THU_PALE = "#660874", "#8B3A9E", "#A966BE", "#C99FD4"
THU_TINT, THU_GRID = "#EFE0F7", "#DDD0E8"
CRIMSON, NAVY, TEAL, AMBER = "#A82038", "#1F5FA8", "#1A6B5A", "#C47A00"
CHARCOAL, NEAR_WHITE = "#2B2B2B", "#FFFFFF"
HM_GREY, FM_BLUE, MGO_TAN = "#6E6E6E", "#3F5F8C", "#D9C28A"

plt.rcParams.update({
    "font.family"      : "serif",
    "font.serif"       : ["Arial", "Liberation Serif"],
    "font.size"        : 11,
    "mathtext.fontset" : "stix",
    "savefig.dpi"      : 300,
    "figure.dpi"       : 150,
})

fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(13.0, 5.6),
                                  gridspec_kw={"wspace": 0.15,
                                               "width_ratios": [1.0, 1.0]})

# ════════════════════════════════════════════════════════════════════════
# Left — Physical device structure
# ════════════════════════════════════════════════════════════════════════
ax_l.set_xlim(0, 10)
ax_l.set_ylim(0, 8.0)
ax_l.set_aspect("equal")
ax_l.axis("off")

# HM channel (bottom slab) — placed higher to leave room for terminal stems
HM_X0, HM_X1, HM_Y0, HM_Y1 = 1.0, 9.0, 2.4, 3.15
ax_l.add_patch(Rectangle((HM_X0, HM_Y0), HM_X1 - HM_X0, HM_Y1 - HM_Y0,
                          facecolor=HM_GREY, edgecolor=CHARCOAL, lw=1.2))
ax_l.text(5.0, (HM_Y0 + HM_Y1) / 2.0, "Heavy Metal (HM) channel",
          ha="center", va="center", fontsize=11.5, fontweight="bold",
          color=NEAR_WHITE)
ax_l.text(HM_X0 - 0.10, HM_Y0 - 0.30, r"$L_{\rm SOT}$, $W_{\rm SOT}$, $T_{\rm SOT}$",
          ha="left", va="top", fontsize=9.5, color=CHARCOAL, style="italic")

# MTJ pillar (centred on HM)
MTJ_CX, MTJ_HW = 5.0, 0.75
# Reference layer (bottom of pillar, on top of HM)
RL_Y0, RL_Y1 = HM_Y1, HM_Y1 + 0.55
ax_l.add_patch(Rectangle((MTJ_CX - MTJ_HW, RL_Y0), 2 * MTJ_HW, RL_Y1 - RL_Y0,
                          facecolor=FM_BLUE, edgecolor=CHARCOAL, lw=1.0))
ax_l.text(MTJ_CX + MTJ_HW + 0.15, (RL_Y0 + RL_Y1) / 2.0, "Reference  (PMA, fixed)",
          ha="left", va="center", fontsize=10, color=CHARCOAL)
# PMA arrow (reference layer, fixed up)
ax_l.annotate("", xy=(MTJ_CX - 0.30, RL_Y1 - 0.08),
              xytext=(MTJ_CX - 0.30, RL_Y0 + 0.08),
              arrowprops=dict(arrowstyle="-|>", color=NEAR_WHITE, lw=1.6))

# MgO barrier
MGO_Y0, MGO_Y1 = RL_Y1, RL_Y1 + 0.25
ax_l.add_patch(Rectangle((MTJ_CX - MTJ_HW, MGO_Y0), 2 * MTJ_HW, MGO_Y1 - MGO_Y0,
                          facecolor=MGO_TAN, edgecolor=CHARCOAL, lw=1.0))
ax_l.text(MTJ_CX + MTJ_HW + 0.15, (MGO_Y0 + MGO_Y1) / 2.0,
          r"MgO  ($t_{\rm ox}$)",
          ha="left", va="center", fontsize=10, color=CHARCOAL)

# Free layer (top of pillar)
FL_Y0, FL_Y1 = MGO_Y1, MGO_Y1 + 0.55
ax_l.add_patch(Rectangle((MTJ_CX - MTJ_HW, FL_Y0), 2 * MTJ_HW, FL_Y1 - FL_Y0,
                          facecolor=FM_BLUE, edgecolor=CHARCOAL, lw=1.0))
ax_l.text(MTJ_CX + MTJ_HW + 0.15, (FL_Y0 + FL_Y1) / 2.0,
          r"Free  (PMA, $t_f$)",
          ha="left", va="center", fontsize=10, color=CHARCOAL)
# PMA arrow (free layer, can flip — show two-headed double arrow)
ax_l.annotate("", xy=(MTJ_CX - 0.30, FL_Y1 - 0.08),
              xytext=(MTJ_CX - 0.30, FL_Y0 + 0.08),
              arrowprops=dict(arrowstyle="<|-|>", color=NEAR_WHITE, lw=1.6))
ax_l.text(MTJ_CX + 0.10, (FL_Y0 + FL_Y1) / 2.0, r"$\mathbf{m}$",
          ha="left", va="center", fontsize=11, color=NEAR_WHITE, fontweight="bold")

# Top electrode (T1)
TOP_Y0, TOP_Y1 = FL_Y1, FL_Y1 + 0.25
ax_l.add_patch(Rectangle((MTJ_CX - MTJ_HW, TOP_Y0), 2 * MTJ_HW, TOP_Y1 - TOP_Y0,
                          facecolor=CHARCOAL, edgecolor=CHARCOAL, lw=1.0))

# T1 wire going up
T1_X, T1_Y = MTJ_CX, TOP_Y1 + 1.85
ax_l.plot([MTJ_CX, MTJ_CX], [TOP_Y1, T1_Y], color=CHARCOAL, lw=1.6)
ax_l.plot(T1_X, T1_Y, "o", ms=9, mec=CHARCOAL, mfc=NEAR_WHITE, mew=1.5, zorder=3)
ax_l.text(T1_X, T1_Y + 0.30, r"T$_1$  ($V_1$)", ha="center", va="bottom",
          fontsize=12, fontweight="bold", color=CHARCOAL)

# T2 wire (left side of HM)
T2_X = HM_X0 - 0.55
ax_l.plot([HM_X0, T2_X], [HM_Y0 + 0.30, HM_Y0 + 0.30], color=CHARCOAL, lw=1.6)
ax_l.plot([T2_X, T2_X], [HM_Y0 + 0.30, HM_Y0 + 0.30 - 0.80], color=CHARCOAL, lw=1.6)
ax_l.plot(T2_X, HM_Y0 + 0.30 - 0.80, "o", ms=9, mec=CHARCOAL, mfc=NEAR_WHITE,
          mew=1.5, zorder=3)
ax_l.text(T2_X, HM_Y0 + 0.30 - 0.80 - 0.30, r"T$_2$  ($V_2$)",
          ha="center", va="top", fontsize=12, fontweight="bold", color=CHARCOAL)

# T3 wire (right side of HM)
T3_X = HM_X1 + 0.55
ax_l.plot([HM_X1, T3_X], [HM_Y0 + 0.30, HM_Y0 + 0.30], color=CHARCOAL, lw=1.6)
ax_l.plot([T3_X, T3_X], [HM_Y0 + 0.30, HM_Y0 + 0.30 - 0.80], color=CHARCOAL, lw=1.6)
ax_l.plot(T3_X, HM_Y0 + 0.30 - 0.80, "o", ms=9, mec=CHARCOAL, mfc=NEAR_WHITE,
          mew=1.5, zorder=3)
ax_l.text(T3_X, HM_Y0 + 0.30 - 0.80 - 0.30, r"T$_3$  ($V_3$)",
          ha="center", va="top", fontsize=12, fontweight="bold", color=CHARCOAL)

# I_SOT arrow inside the channel — labelled outside the HM block to avoid
# colliding with the "Heavy Metal (HM) channel" caption.
ax_l.annotate("", xy=(7.6, HM_Y1 - 0.20),
              xytext=(2.4, HM_Y1 - 0.20),
              arrowprops=dict(arrowstyle="-|>", color=NEAR_WHITE, lw=2.0))
ax_l.text(5.0, HM_Y0 - 0.30, r"$I_{\rm SOT}$ (write)",
          ha="center", va="top", fontsize=10, color=CHARCOAL, fontweight="bold")

# Spin-accumulation arrow under the free layer (perpendicular to current)
ax_l.annotate("", xy=(MTJ_CX + 0.05, HM_Y1 + 0.55),
              xytext=(MTJ_CX - 0.05, HM_Y1 + 0.55),
              arrowprops=dict(arrowstyle="-|>", color=THU_DEEP, lw=1.4))
ax_l.text(MTJ_CX - 1.50, HM_Y1 + 0.65, r"$\hat{\sigma}$",
          ha="right", va="bottom", fontsize=11, color=THU_DEEP, fontweight="bold")

# I_MTJ small arrow (read path, downward through MTJ)
ax_l.annotate("", xy=(MTJ_CX + 0.55, FL_Y0 - 0.10),
              xytext=(MTJ_CX + 0.55, T1_Y - 0.50),
              arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.2))
ax_l.text(MTJ_CX + 0.95, (FL_Y0 + T1_Y - 0.5) / 2, r"$I_{\rm MTJ}$ (read)",
          ha="left", va="center", fontsize=9.5, color=CRIMSON)

# Sub-title without panel tag
ax_l.text(0.5, 1.00, "Physical device structure",
          transform=ax_l.transAxes, ha="center", va="top",
          fontsize=13, fontweight="bold", color=CHARCOAL)


# ════════════════════════════════════════════════════════════════════════
# Panel (b) — T-equivalent circuit
# ════════════════════════════════════════════════════════════════════════
ax_r.set_xlim(0, 10)
ax_r.set_ylim(0, 8.0)
ax_r.set_aspect("equal")
ax_r.axis("off")


def zigzag(ax, x0, y0, x1, y1, n_zigs=8, amp=0.18, color=CHARCOAL, lw=1.6):
    pts = []
    for i in range(n_zigs + 1):
        t = i / n_zigs
        x = x0 + t * (x1 - x0)
        y = y0 + t * (y1 - y0)
        if 0 < i < n_zigs:
            dx, dy = (x1 - x0) / n_zigs, (y1 - y0) / n_zigs
            length = np.hypot(dx, dy)
            nx, ny = -dy / length, dx / length
            sign = +1 if (i % 2 == 1) else -1
            x += sign * amp * nx
            y += sign * amp * ny
        pts.append((x, y))
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=color, lw=lw, solid_joinstyle="miter")


def vresistor(ax, x_center, y_top, y_bot, n_zigs=10, amp=0.20, color=CHARCOAL,
              lw=1.7, variable=False):
    zigzag(ax, x_center, y_top, x_center, y_bot,
           n_zigs=n_zigs, amp=amp, color=color, lw=lw)
    if variable:
        ax.annotate("", xy=(x_center + 0.55, y_top - 0.25),
                    xytext=(x_center - 0.35, y_bot + 0.15),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.2))


def terminal_dot(ax, x, y, label, label_offset=(0.0, 0.30), fontsize=12):
    ax.plot(x, y, "o", ms=9, mec=CHARCOAL, mfc=NEAR_WHITE, mew=1.5, zorder=3)
    ax.text(x + label_offset[0], y + label_offset[1], label,
            ha="center", va="center", fontsize=fontsize, fontweight="bold",
            color=CHARCOAL)

# T1 terminal (top) — aligned vertically with left panel T1
T1_RX, T1_RY = 5.0, 7.0
terminal_dot(ax_r, T1_RX, T1_RY, r"T$_1$  ($V_1$)", label_offset=(0.0, 0.30))
ax_r.plot([T1_RX - 0.45, T1_RX + 0.45], [T1_RY - 0.50, T1_RY - 0.50],
          color=CHARCOAL, lw=2.6)
ax_r.plot([T1_RX, T1_RX], [T1_RY - 0.08, T1_RY - 0.50],
          color=CHARCOAL, lw=1.6)
# I_MTJ arrow
ax_r.annotate("", xy=(T1_RX - 0.55, T1_RY - 0.85),
              xytext=(T1_RX - 0.55, T1_RY - 0.18),
              arrowprops=dict(arrowstyle="-|>", color=CHARCOAL, lw=1.5))
ax_r.text(T1_RX - 1.10, T1_RY - 0.55, r"$I_{\rm MTJ}$",
          ha="center", va="center", fontsize=11)

# MTJ resistor (variable-R)
MTJ_Y_TOP, MTJ_Y_BOT = 6.25, 3.80
vresistor(ax_r, T1_RX, MTJ_Y_TOP, MTJ_Y_BOT, n_zigs=10, amp=0.22,
          color=CHARCOAL, lw=1.8, variable=True)
ax_r.text(T1_RX - 1.05, (MTJ_Y_TOP + MTJ_Y_BOT) / 2.0, r"$R_{\rm MTJ}$",
          ha="center", va="center", fontsize=12.5)

# V_MTJ measurement arrow
v_arrow_x = T1_RX + 0.85
ax_r.annotate("", xy=(v_arrow_x, T1_RY - 0.20),
              xytext=(v_arrow_x, MTJ_Y_BOT - 0.10),
              arrowprops=dict(arrowstyle="<-", color=CHARCOAL, lw=1.4))
ax_r.text(v_arrow_x + 0.45, (T1_RY + MTJ_Y_BOT) / 2.0 - 0.10,
          r"$V_{\rm MTJ}$", ha="center", va="center", fontsize=12)

# Mid node N
MID_RX, MID_RY = T1_RX, 3.45
ax_r.plot(MID_RX, MID_RY, "o", ms=8, mec=CHARCOAL, mfc=CHARCOAL, zorder=3)
ax_r.plot([T1_RX, T1_RX], [MTJ_Y_BOT - 0.05, MID_RY], color=CHARCOAL, lw=1.6)
ax_r.text(MID_RX + 0.20, MID_RY - 0.10, r"N", ha="left", va="top",
          fontsize=11, fontweight="bold", color=CHARCOAL)

# Horizontal H-resistor zigzag
def hzigzag(ax, x0, x1, y, n_zigs=6, amp=0.18, color=CHARCOAL, lw=1.6):
    pts = []
    for i in range(n_zigs + 1):
        t = i / n_zigs
        x = x0 + t * (x1 - x0)
        yy = y
        if 0 < i < n_zigs:
            yy = y + (amp if (i % 2 == 1) else -amp)
        pts.append((x, yy))
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=color, lw=lw)

LEFT_RX, RIGHT_RX = 1.4, 8.6
ax_r.plot([MID_RX, MID_RX - 1.30], [MID_RY, MID_RY], color=CHARCOAL, lw=1.6)
ax_r.plot([MID_RX, MID_RX + 1.30], [MID_RY, MID_RY], color=CHARCOAL, lw=1.6)
hzigzag(ax_r, MID_RX - 1.30, MID_RX - 2.90, MID_RY)
hzigzag(ax_r, MID_RX + 1.30, MID_RX + 2.90, MID_RY)
ax_r.text(MID_RX - 2.10, MID_RY + 0.45, r"$R_{\rm SOT}/2$",
          ha="center", va="center", fontsize=11.5)
ax_r.text(MID_RX + 2.10, MID_RY + 0.45, r"$R_{\rm SOT}/2$",
          ha="center", va="center", fontsize=11.5)
ax_r.plot([MID_RX - 2.90, LEFT_RX],  [MID_RY, MID_RY], color=CHARCOAL, lw=1.6)
ax_r.plot([MID_RX + 2.90, RIGHT_RX], [MID_RY, MID_RY], color=CHARCOAL, lw=1.6)
terminal_dot(ax_r, LEFT_RX,  MID_RY, r"T$_2$  ($V_2$)", label_offset=(-0.05, +0.45))
terminal_dot(ax_r, RIGHT_RX, MID_RY, r"T$_3$  ($V_3$)", label_offset=(+0.05, +0.45))

# I_2 / I_3 arrows
ax_r.annotate("", xy=(LEFT_RX + 0.55, MID_RY),
              xytext=(LEFT_RX + 0.10, MID_RY),
              arrowprops=dict(arrowstyle="-|>", color=CHARCOAL, lw=1.5))
ax_r.text(LEFT_RX + 0.65, MID_RY - 0.35, r"$I_2$",
          ha="center", va="center", fontsize=11)
ax_r.annotate("", xy=(RIGHT_RX - 0.10, MID_RY),
              xytext=(RIGHT_RX - 0.55, MID_RY),
              arrowprops=dict(arrowstyle="-|>", color=CHARCOAL, lw=1.5))
ax_r.text(RIGHT_RX - 0.65, MID_RY - 0.35, r"$I_3$",
          ha="center", va="center", fontsize=11)

# Sub-title without panel tag
ax_r.text(0.5, 1.00, "T-equivalent resistor network",
          transform=ax_r.transAxes, ha="center", va="top",
          fontsize=13, fontweight="bold", color=CHARCOAL)


# ════════════════════════════════════════════════════════════════════════
# Save + sync
# ════════════════════════════════════════════════════════════════════════
plt.tight_layout()
out_path = Path(__file__).resolve().parent / "fig_01_t_circuit.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white",
            pad_inches=0.15)
plt.close()
chapter_fig = (Path(__file__).resolve().parent.parent /
               "article" / "00_chapter_drafts" / "figures" / "fig_01_t_circuit.png")
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)
print(f"Saved  {out_path}")
print(f"Synced {chapter_fig}")
