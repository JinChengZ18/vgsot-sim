"""
Figure 2.2: layered behavioral model for probabilistic sMTJ computing.

The layout is a compact 2 x 2 journal-style composite figure. It preserves the
logical order from computational abstraction to device mechanisms while keeping
the aspect ratio suitable for thesis text pages.

Outputs:
  demo/fig_02_behavioral_layers.png
  article/00_chapter_drafts/figures/fig_02_behavioral_layers.png
"""
from __future__ import annotations

from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, Ellipse, FancyArrowPatch
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle


# Journal-style palette: restrained, high-contrast, and print friendly.
PURPLE = "#5A136F"
PURPLE_2 = "#8B5AA6"
PURPLE_TINT = "#F0E8F4"
BLUE = "#1F5EA8"
TEAL = "#087566"
AMBER = "#B87900"
RED = "#A71E35"
INK = "#202020"
MID = "#666666"
GRID = "#D8D3DC"
LINE = "#2D2D2D"
PANEL = "#FBFAFC"
WHITE = "#FFFFFF"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans"],
    "font.size": 9.2,
    "mathtext.fontset": "stixsans",
    "axes.unicode_minus": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})


fig, ax = plt.subplots(figsize=(8.35, 5.90))
fig.subplots_adjust(left=0.018, right=0.985, bottom=0.025, top=0.985)
ax.set_xlim(0, 10.0)
ax.set_ylim(0, 6.5)
ax.axis("off")


def arrow(x0: float, y0: float, x1: float, y1: float, *,
          color: str = PURPLE, lw: float = 1.05, ms: float = 8.5,
          rad: float = 0.0) -> None:
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle="-|>", mutation_scale=ms, lw=lw, color=color,
        shrinkA=1.2, shrinkB=1.2, connectionstyle=f"arc3,rad={rad}",
    ))


def text(x: float, y: float, s: str, *, fs: float = 8.2, color: str = INK,
         weight: str = "normal", ha: str = "center", va: str = "center",
         style: str = "normal") -> None:
    ax.text(x, y, s, fontsize=fs, color=color, fontweight=weight,
            ha=ha, va=va, fontstyle=style)


def panel(x: float, y: float, w: float, h: float,
          tag: str, title: str) -> tuple[float, float, float, float]:
    ax.add_patch(Rectangle((x, y), w, h, facecolor=PANEL,
                           edgecolor=LINE, lw=0.72))
    text(x + 0.18, y + h - 0.19, f"{tag} {title}",
         fs=10.0, weight="bold", ha="left")
    ax.plot([x + 0.14, x + w - 0.14], [y + h - 0.34, y + h - 0.34],
            color=GRID, lw=0.55)
    return x + 0.18, y + 0.13, w - 0.36, h - 0.47


def box(x: float, y: float, w: float, h: float, title: str,
        subtitle: str = "", *, fc: str = WHITE, fs: float = 8.8,
        sub_fs: float = 7.5, lw: float = 0.72) -> None:
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.018,rounding_size=0.035",
        facecolor=fc, edgecolor=LINE, lw=lw,
    ))
    if subtitle:
        text(x + w / 2, y + h * 0.66, title, fs=fs, weight="bold")
        text(x + w / 2, y + h * 0.31, subtitle, fs=sub_fs)
    else:
        text(x + w / 2, y + h / 2, title, fs=fs, weight="bold")


def draw_waveforms(x: float, y: float, w: float, h: float) -> None:
    ax.add_patch(Rectangle((x, y), w, h, facecolor=WHITE, edgecolor=LINE, lw=0.68))
    rows = [
        (r"$V_{\rm MTJ}$", TEAL, 0.74, 0.46, 0.76),
        (r"$I_{\rm SOT}$", PURPLE, 0.50, 0.28, 0.70),
        (r"$t_w$", BLUE, 0.27, 0.38, 0.66),
    ]
    for lab, color, frac, t0, t1 in rows:
        yy = y + frac * h
        x0, x1 = x + 0.30 * w, x + 0.93 * w
        base, top = yy - 0.055 * h, yy + 0.145 * h
        xs = [x0, x + t0 * w, x + t0 * w, x + t1 * w, x + t1 * w, x1]
        ys = [base, base, top, top, base, base]
        ax.plot(xs, ys, color=color, lw=1.35)
        ax.plot([x0, x1], [base, base], color=GRID, lw=0.48)
        text(x + 0.08 * w, yy + 0.01, lab, fs=8.2, color=color, ha="left")
    text(x + w / 2, y + 0.055 * h, "drives", fs=7.2, color=MID)


def draw_bitstream(x: float, y: float, w: float, h: float) -> None:
    bits = np.array([1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0])
    gap = 0.018
    cell = (w - (len(bits) - 1) * gap) / len(bits)
    for i, bit in enumerate(bits):
        xx = x + i * (cell + gap)
        ax.add_patch(Rectangle((xx, y), cell, h,
                               facecolor=PURPLE if bit else "#E9E7EB",
                               edgecolor=LINE, lw=0.35))
        text(xx + cell / 2, y + h / 2, str(bit), fs=6.3,
             color=WHITE if bit else INK, weight="bold")


def draw_sigmoid_plot(x: float, y: float, w: float, h: float) -> None:
    ax.add_patch(Rectangle((x, y), w, h, facecolor=WHITE, edgecolor=LINE, lw=0.68))
    px0, py0 = x + 0.34, y + 0.27
    pw, ph = w - 0.60, h - 0.62
    u = np.linspace(-4.0, 4.0, 320)
    psw = 1 / (1 + np.exp(-1.55 * u))
    X = px0 + (u - u.min()) / (u.max() - u.min()) * pw
    Y = py0 + psw * ph
    ax.add_patch(Rectangle((px0 + 0.43 * pw, py0), 0.16 * pw, ph,
                           facecolor=PURPLE_TINT, edgecolor="none"))
    ax.plot([px0, px0 + pw], [py0, py0], color=LINE, lw=0.65)
    ax.plot([px0, px0], [py0, py0 + ph], color=LINE, lw=0.65)
    ax.plot(X, Y, color=PURPLE, lw=1.85)
    uth = px0 + 0.50 * pw
    ax.plot([uth, uth], [py0, py0 + ph], color=GRID, lw=0.75, ls="--")
    ax.plot([uth - 0.36, uth + 0.48], [py0 + 0.35 * ph, py0 + 0.66 * ph],
            color=RED, lw=1.1)
    text(px0 - 0.08, py0 + ph / 2, r"$P_{\rm sw}$", fs=7.5, ha="right")
    text(px0 + pw / 2, py0 - 0.13, r"$u-u_{\rm th}$", fs=7.5)
    text(uth + 0.04, py0 + 0.09, r"$u_{\rm th}$", fs=7.0, color=MID, ha="left")
    text(uth + 0.26, py0 + 0.68 * ph, r"$\beta_s$", fs=7.5, color=RED, ha="left")
    text(uth - 0.02, py0 + 0.86 * ph, "critical\nregion",
         fs=6.7, color=PURPLE)
    ax.text(x + 0.26, y + h - 0.14,
            r"$P_{\rm sw}(u)\simeq[1+\exp[-\beta_s(u-u_{\rm th})]]^{-1}$",
            ha="left", va="top", fontsize=7.7, color=INK)


def draw_device_stack(x: float, y: float, w: float, h: float) -> None:
    hm = Polygon([
        (x + 0.10 * w, y + 0.20 * h),
        (x + 0.76 * w, y + 0.20 * h),
        (x + 0.90 * w, y + 0.33 * h),
        (x + 0.24 * w, y + 0.33 * h),
    ], closed=True, facecolor="#696969", edgecolor=LINE, lw=0.65)
    ax.add_patch(hm)
    text(x + 0.49 * w, y + 0.265 * h, "HM channel", fs=6.9, color=WHITE)
    arrow(x + 0.18 * w, y + 0.10 * h, x + 0.72 * w, y + 0.10 * h,
          color=INK, lw=0.75, ms=7)
    text(x + 0.45 * w, y + 0.015 * h, r"$I_{\rm SOT}$", fs=7.2)

    cx = x + 0.48 * w
    stack_w = 0.34 * w
    yy = y + 0.33 * h
    layers = [
        ("Free layer", "#4772A5", 0.19 * h),
        ("MgO", "#DCC783", 0.13 * h),
        ("Pinned layer", "#752A86", 0.19 * h),
    ]
    for name, color, hh in layers:
        ax.add_patch(Rectangle((cx - stack_w / 2, yy), stack_w, hh,
                               facecolor=color, edgecolor=LINE, lw=0.6))
        text(cx, yy + hh / 2, name, fs=6.6,
             color=WHITE if name != "MgO" else INK)
        yy += hh
    ax.add_patch(Ellipse((cx, yy), stack_w, 0.10 * h,
                         facecolor="#752A86", edgecolor=LINE, lw=0.6))

    arrow(cx + 0.23 * w, y + 0.36 * h, cx + 0.23 * w, y + 0.90 * h,
          color=PURPLE, lw=0.85, ms=7)
    text(cx + 0.27 * w, y + 0.63 * h, r"$V_{\rm MTJ}$", fs=7.1,
         color=PURPLE, ha="left")
    arrow(cx - 0.26 * w, y + 0.35 * h, cx - 0.26 * w, y + 0.86 * h,
          color=LINE, lw=0.75, ms=7)
    text(cx - 0.31 * w, y + 0.61 * h, "PMA", fs=7.1, weight="bold", ha="right")


def draw_thermal_icon(x: float, y: float, r: float) -> None:
    center = (x, y)
    ax.add_patch(Circle(center, r, facecolor="#FFF7E8", edgecolor=AMBER, lw=0.68))
    for deg in np.linspace(20, 340, 8):
        th = np.deg2rad(deg)
        arrow(x + 1.55 * r * np.cos(th), y + 1.55 * r * np.sin(th),
              x + 2.25 * r * np.cos(th), y + 2.25 * r * np.sin(th),
              color=AMBER, lw=0.55, ms=5.3)
    text(x, y, r"$H_{\rm th}$", fs=6.3, color=AMBER, weight="bold")


def draw_llg_icon(x: float, y: float, s: float) -> None:
    cx, cy = x + 0.42 * s, y + 0.50 * s
    ax.add_patch(Circle((cx, cy), 0.34 * s, facecolor="#F8F6FA",
                        edgecolor=GRID, lw=0.75))
    ax.add_patch(Arc((cx, cy), 0.58 * s, 0.36 * s, theta1=30, theta2=325,
                     color=PURPLE, lw=1.35))
    arrow(cx, cy, cx + 0.25 * s, cy + 0.25 * s, color=BLUE, lw=1.0, ms=7)
    arrow(cx, cy, cx - 0.08 * s, cy + 0.32 * s, color=TEAL, lw=1.0, ms=7)
    text(cx + 0.29 * s, cy + 0.28 * s, r"$\mathbf{m}$", fs=6.8, color=BLUE, ha="left")
    text(cx - 0.12 * s, cy + 0.36 * s, r"$\mathbf{H}_{\rm eff}$",
         fs=6.7, color=TEAL, ha="right")


def draw_barrier_plot(x: float, y: float, w: float, h: float) -> None:
    ax.add_patch(Rectangle((x, y), w, h, facecolor=WHITE, edgecolor=LINE, lw=0.68))
    text(x + w / 2, y + h - 0.16, "modulated barrier", fs=8.0, weight="bold")
    px0, py0 = x + 0.28, y + 0.27
    pw, ph = w - 0.50, h - 0.60
    q = np.linspace(-1.25, 1.25, 240)
    high = 0.82 * (q**2 - 1.0) ** 2
    low = 0.50 * (q**2 - 1.0) ** 2 + 0.05 * q
    lo, hi = min(high.min(), low.min()), max(high.max(), low.max())
    X = px0 + (q - q.min()) / (q.max() - q.min()) * pw
    Yh = py0 + (high - lo) / (hi - lo) * ph
    Yl = py0 + (low - lo) / (hi - lo) * ph
    ax.plot(X, Yh, color=MID, lw=1.0, ls=(0, (4, 2)))
    ax.plot(X, Yl, color=PURPLE, lw=1.65)
    mid = len(q) // 2
    arrow(X[mid] + 0.08, Yh[mid] - 0.02, X[mid] + 0.08, Yl[mid] + 0.03,
          color=RED, lw=0.8, ms=6)
    text(X[mid] + 0.17, (Yh[mid] + Yl[mid]) / 2, r"$\Delta E_b$",
         fs=6.8, color=RED, ha="left")
    text(x + w / 2, y + 0.10,
         r"$\Delta E_b=\Delta E_0-\eta_VV_{\rm MTJ}-\eta_II_{\rm SOT}$",
         fs=6.6)


def draw_escape_plot(x: float, y: float, w: float, h: float) -> None:
    ax.add_patch(Rectangle((x, y), w, h, facecolor=WHITE, edgecolor=LINE, lw=0.68))
    text(x + w / 2, y + h - 0.20, "Neel-Brown\nescape", fs=6.2, weight="bold")
    px0, py0 = x + 0.22, y + 0.38
    pw, ph = w - 0.40, h - 0.82
    q = np.linspace(-1.18, 1.18, 220)
    E = (q**2 - 1.0) ** 2
    X = px0 + (q - q.min()) / (q.max() - q.min()) * pw
    Y = py0 + (E - E.min()) / (E.max() - E.min()) * ph
    ax.plot(X, Y, color=LINE, lw=1.05)
    i0 = 36
    ax.add_patch(Circle((X[i0], Y[i0]), 0.035, facecolor=PURPLE,
                        edgecolor=WHITE, lw=0.45))
    arrow(X[i0] + 0.04, Y[i0] + 0.02, px0 + 0.60 * pw, py0 + 0.76 * ph,
          color=PURPLE, lw=1.0, ms=7, rad=-0.10)
    arrow(px0 + 0.60 * pw, py0 + 0.76 * ph, px0 + 0.82 * pw, py0 + 0.22 * ph,
          color=PURPLE, lw=1.0, ms=7, rad=-0.08)
    text(x + w / 2, y + 0.22, r"$\tau=\tau_0e^{\Delta_{\rm eff}}$", fs=6.1)
    text(x + w / 2, y + 0.08, r"$P_{\rm sw}=1-\exp(-t_w/\tau)$", fs=6.1)


# Balanced composite layout.
W, H = 4.72, 2.94
px1, px2 = 0.18, 5.10
py_top, py_bot = 3.28, 0.18

ca = panel(px1, py_top, W, H, "(a)", "Computational layer")
cb = panel(px2, py_top, W, H, "(b)", "Behavioral layer")
cc = panel(px1, py_bot, W, H, "(c)", "Physical-model layer")
cd = panel(px2, py_bot, W, H, "(d)", "Device layer")


# (a) Computational layer.
x, y, w, h = ca
draw_waveforms(x + 0.04, y + 1.08, 1.25, 1.27)
text(x + 1.58, y + 1.78, r"$(t_w,I_{\rm SOT},V_{\rm MTJ})$", fs=8.1)
arrow(x + 2.00, y + 1.78, x + 2.24, y + 1.78)
box(x + 2.26, y + 1.42, 0.86, 0.74, "map", r"$u=g(\cdot)$",
    fs=7.8, sub_fs=7.1)
arrow(x + 3.14, y + 1.78, x + 3.36, y + 1.78)
box(x + 3.38, y + 1.42, 0.88, 0.74, "sampler",
    r"$b_k\sim{\rm Bern}(P_{\rm sw})$", fs=7.6, sub_fs=6.2)
text(x + 2.73, y + 1.04, r"$P_{\rm sw}(u)$", fs=8.8, weight="bold")
arrow(x + 2.92, y + 1.10, x + 3.66, y + 1.37, rad=0.14)
draw_bitstream(x + 1.38, y + 0.45, 2.52, 0.31)
text(x + 2.64, y + 0.23, "programmable stochastic bit source", fs=6.8, color=MID)
box(x + 0.06, y + 0.30, 1.06, 0.62, "bit stream", r"$m\in\{0,1\}$",
    fs=7.2, sub_fs=6.5, fc="#F7F5F8")
arrow(x + 1.16, y + 0.61, x + 1.34, y + 0.61)


# (b) Behavioral layer.
x, y, w, h = cb
draw_sigmoid_plot(x + 0.06, y + 0.45, 2.62, 1.93)
box(x + 2.86, y + 1.66, 0.72, 0.58, "drive", r"$u,u_{\rm th}$",
    fs=7.2, sub_fs=6.5)
box(x + 3.72, y + 1.66, 0.68, 0.58, "slope", r"$\beta_s$",
    fs=7.2, sub_fs=6.8)
arrow(x + 3.60, y + 1.95, x + 3.70, y + 1.95, color=MID, lw=0.75, ms=6)
ax.add_patch(Rectangle((x + 2.80, y + 0.34), 1.54, 0.94,
                       facecolor=WHITE, edgecolor=LINE, lw=0.62))
text(x + 3.57, y + 1.12, "parameter map", fs=7.0, weight="bold")
for i, (lab, col) in enumerate([(r"$\Delta$", BLUE), (r"$\xi$", TEAL), (r"$\theta_{\rm SH}$", PURPLE)]):
    cx = x + 3.08 + i * 0.32
    ax.add_patch(Circle((cx, y + 0.79), 0.115, facecolor=col, edgecolor=WHITE, lw=0.45))
    text(cx, y + 0.79, lab, fs=5.3, color=WHITE, weight="bold")
arrow(x + 3.80, y + 0.79, x + 4.00, y + 0.79, color=MID, lw=0.65, ms=5.5)
text(x + 4.08, y + 0.79, r"$u_{\rm th},\beta_s$", fs=6.0, ha="left")
text(x + 3.57, y + 0.52, r"$f(\Delta,\xi,\theta_{\rm SH})$", fs=6.2)


# (c) Physical-model layer.
x, y, w, h = cc
box(x + 0.06, y + 1.58, 1.52, 0.66, "LLG + torques",
    r"$\dot{\bf m}=-\gamma\mu_0{\bf m}\times{\bf H}_{\rm eff}+\tau_{\rm SOT}+\tau_{\rm th}$",
    fs=7.2, sub_fs=5.9)
draw_llg_icon(x + 0.10, y + 0.45, 1.34)
arrow(x + 1.60, y + 1.48, x + 1.76, y + 1.48)
draw_barrier_plot(x + 1.76, y + 0.86, 1.36, 1.34)
arrow(x + 3.16, y + 1.48, x + 3.24, y + 1.48)
draw_escape_plot(x + 3.26, y + 0.86, 1.08, 1.34)
ax.add_patch(Rectangle((x + 1.34, y + 0.18), 2.96, 0.45,
                       facecolor="#F7F5F8", edgecolor=LINE, lw=0.62))
text(x + 2.82, y + 0.405,
     r"unified switching probability  $P_{\rm sw}(t_w,I_{\rm SOT},V_{\rm MTJ})$",
     fs=6.9)
arrow(x + 2.45, y + 0.84, x + 2.70, y + 0.64, color=PURPLE, lw=0.85, ms=6)
arrow(x + 3.84, y + 0.84, x + 3.70, y + 0.64, color=PURPLE, lw=0.85, ms=6)


# (d) Device layer.
x, y, w, h = cd
draw_device_stack(x + 0.00, y + 0.22, 2.22, 2.14)
box(x + 2.34, y + 1.66, 0.90, 0.56, "VCMA",
    r"$V_{\rm MTJ}$ across MgO", fs=7.2, sub_fs=6.2)
box(x + 3.38, y + 1.66, 0.78, 0.56, "SOT",
    r"$I_{\rm SOT}$ in HM", fs=7.2, sub_fs=6.2)
box(x + 2.52, y + 0.64, 1.28, 0.60, "thermal activation",
    r"$\mathbf{H}_{\rm th}(t)$", fs=7.0, sub_fs=6.4)
draw_thermal_icon(x + 4.18, y + 0.94, 0.155)
arrow(x + 2.23, y + 1.34, x + 2.50, y + 0.99, color=MID, lw=0.75, ms=6)
arrow(x + 2.24, y + 1.78, x + 2.32, y + 1.94, color=PURPLE, lw=0.85, ms=6)
arrow(x + 2.24, y + 0.89, x + 2.50, y + 0.89, color=AMBER, lw=0.85, ms=6)
text(x + 2.90, y + 0.22,
     "PMA free layer, VCMA barrier tuning,\nSOT drive, and thermal fluctuation",
     fs=6.4, color=MID)


out_path = Path(__file__).resolve().parent / "fig_02_behavioral_layers.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white", pad_inches=0.06)
plt.close(fig)

chapter_fig = (
    Path(__file__).resolve().parent.parent
    / "article" / "00_chapter_drafts" / "figures" / "fig_02_behavioral_layers.png"
)
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)

print(f"Saved  {out_path}")
print(f"Synced {chapter_fig}")
