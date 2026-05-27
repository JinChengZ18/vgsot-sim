"""
Figure 2.3: thermal and non-ideal effects in the sMTJ behavioral model.

Top-journal style revision: fewer words, larger typography, stronger visual
hierarchy, and compact 2 x 2 panel geometry.

Outputs:
  demo/Chapter02_local_03.png
  article/00_chapter_drafts/figs/Chapter02_local_03.png
"""
from __future__ import annotations

from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib.patches import Rectangle


PURPLE = "#5A136F"
PURPLE_TINT = "#F3EAF7"
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
    "font.size": 12.0,
    "mathtext.fontset": "stixsans",
    "axes.unicode_minus": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})


fig, ax = plt.subplots(figsize=(7.35, 5.05))
fig.subplots_adjust(left=0.018, right=0.985, bottom=0.025, top=0.985)
ax.set_xlim(0, 10)
ax.set_ylim(0, 6.0)
ax.axis("off")


def text(x: float, y: float, s: str, *, fs: float = 10.0,
         color: str = INK, weight: str = "normal",
         ha: str = "center", va: str = "center",
         style: str = "normal") -> None:
    ax.text(x, y, s, fontsize=fs, color=color, fontweight=weight,
            ha=ha, va=va, fontstyle=style)


def arrow(x0: float, y0: float, x1: float, y1: float, *,
          color: str = PURPLE, lw: float = 1.05, ms: float = 8.5,
          rad: float = 0.0) -> None:
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms,
        lw=lw, color=color, shrinkA=1.2, shrinkB=1.2,
        connectionstyle=f"arc3,rad={rad}",
    ))


def panel(x: float, y: float, w: float, h: float,
          tag: str, title: str) -> tuple[float, float, float, float]:
    ax.add_patch(Rectangle((x, y), w, h, facecolor=PANEL,
                           edgecolor=LINE, lw=0.70))
    text(x + 0.16, y + h - 0.17, f"{tag} {title}",
         fs=12.5, weight="bold", ha="left")
    ax.plot([x + 0.12, x + w - 0.12], [y + h - 0.31, y + h - 0.31],
            color=GRID, lw=0.55)
    return x + 0.14, y + 0.12, w - 0.28, h - 0.40


def box(x: float, y: float, w: float, h: float, title: str,
        subtitle: str = "", *, fc: str = WHITE, ec: str = LINE,
        fs: float = 9.2, sub_fs: float = 8.0) -> None:
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.014,rounding_size=0.030",
        facecolor=fc, edgecolor=ec, lw=0.68,
    ))
    if subtitle:
        text(x + w / 2, y + h * 0.64, title, fs=fs, weight="bold")
        text(x + w / 2, y + h * 0.30, subtitle, fs=sub_fs)
    else:
        text(x + w / 2, y + h / 2, title, fs=fs, weight="bold")


def draw_heat(x: float, y: float, w: float, h: float) -> None:
    text(x + 0.08, y + h - 0.18,
         r"$C_v\dot T=\lambda\nabla^2T+Q_{\rm J}$",
         fs=11.3, ha="left")
    cx, cy = x + 0.66, y + 0.90
    ax.add_patch(Circle((cx, cy), 0.24, facecolor="#FFF5E4",
                        edgecolor=AMBER, lw=0.8))
    for deg in np.linspace(20, 340, 8):
        th = np.deg2rad(deg)
        arrow(cx + 0.32 * np.cos(th), cy + 0.32 * np.sin(th),
              cx + 0.50 * np.cos(th), cy + 0.50 * np.sin(th),
              color=AMBER, lw=0.58, ms=5.5)
    text(cx, cy, r"$Q_{\rm J}$", fs=9.3, color=AMBER, weight="bold")

    px, py = x + 1.30, y + 0.28
    pw, ph = w - 1.55, h - 0.68
    ax.plot([px, px + pw], [py, py], color=LINE, lw=0.76)
    ax.plot([px, px], [py, py + ph], color=LINE, lw=0.76)
    t = np.linspace(0, 1, 220)
    hot = 1.0 - np.exp(-3.2 * t)
    cool = 0.62 * (1.0 - np.exp(-5.8 * t))
    ax.plot(px + t * pw, py + hot * ph, color=RED, lw=1.95)
    ax.plot(px + t * pw, py + cool * ph, color=AMBER, lw=1.35, ls=(0, (4, 2)))
    text(px - 0.08, py + ph / 2, r"$T(t)$", fs=10.5, ha="right")
    text(px + pw / 2, py - 0.12, r"$t$", fs=10.5)
    text(px + 0.65 * pw, py + 0.62 * ph, r"$\Delta T$", fs=9.8, color=RED)


def draw_channels(x: float, y: float, w: float, h: float) -> None:
    rows = [
        ("materials", r"$M_s,K_i,\eta$", BLUE, -0.30),
        ("transport", r"$R_{\rm MTJ}$", TEAL, 0.55),
        ("demag", r"$\mathbf{H}_{\rm D}$", AMBER, 0.05),
    ]
    for i, (name, eq, col, trend) in enumerate(rows):
        yy = y + h - 0.62 - i * 0.56
        box(x + 0.04, yy, 1.43, 0.42, name, eq, fs=9.4, sub_fs=9.0, ec=col)
        px, py = x + 1.75, yy + 0.09
        pw, ph = w - 1.95, 0.28
        ax.plot([px, px + pw], [py, py], color=GRID, lw=0.48)
        u = np.linspace(0, 1, 90)
        if name == "materials":
            curve = 0.80 + trend * u
        elif name == "transport":
            curve = 0.18 + trend * u**1.45
        else:
            curve = 0.48 + 0.18 * np.sin(2 * np.pi * u)
        ax.plot(px + u * pw, py + curve * ph, color=col, lw=1.45)
    text(x + w / 2, y + 0.16, r"$T$-corrected inputs", fs=9.5, color=MID)


def draw_llg(x: float, y: float, w: float, h: float) -> None:
    cx, cy = x + 0.78, y + 0.98
    ax.add_patch(Circle((cx, cy), 0.53, facecolor="#F7F5F9",
                        edgecolor=GRID, lw=0.82))
    ax.add_patch(Arc((cx, cy), 0.88, 0.54, theta1=25, theta2=325,
                     color=PURPLE, lw=1.65))
    arrow(cx, cy, cx + 0.36, cy + 0.36, color=BLUE, lw=1.05, ms=7.5)
    arrow(cx, cy, cx - 0.14, cy + 0.44, color=TEAL, lw=1.05, ms=7.5)
    text(cx + 0.42, cy + 0.38, r"$\mathbf{m}$", fs=9.2, color=BLUE, ha="left")
    text(cx - 0.18, cy + 0.47, r"$\mathbf{H}_{\rm eff}(T)$",
         fs=8.8, color=TEAL, ha="right")

    box(x + 1.58, y + 0.94, w - 1.68, 0.74, "stochastic LLG",
        r"$\dot{\bf m}=-\gamma\mu_0{\bf m}\times{\bf H}_{\rm eff}"
        r"+\tau_{\rm SOT}+\tau_{\rm th}$",
        fs=10.2, sub_fs=8.6, ec=PURPLE)
    box(x + 1.58, y + 0.28, 1.22, 0.46, r"$\Delta t$", "update",
        fs=9.7, sub_fs=8.0, fc="#F7F5F9")
    box(x + 3.00, y + 0.28, 1.22, 0.46, r"$t_w$", "sample",
        fs=9.7, sub_fs=8.0, fc="#F7F5F9")
    arrow(x + 2.82, y + 0.51, x + 2.98, y + 0.51, color=MID, lw=0.70, ms=6)


def draw_psw(x: float, y: float, w: float, h: float) -> None:
    px, py = x + 0.44, y + 0.28
    pw, ph = w - 0.78, h - 0.62
    ax.plot([px, px + pw], [py, py], color=LINE, lw=0.76)
    ax.plot([px, px], [py, py + ph], color=LINE, lw=0.76)
    u = np.linspace(-4, 4, 260)
    nominal = 1 / (1 + np.exp(-1.55 * u))
    shifted = 1 / (1 + np.exp(-0.95 * (u + 0.70)))
    X = px + (u - u.min()) / (u.max() - u.min()) * pw
    ax.plot(X, py + nominal * ph, color=INK, lw=1.85)
    ax.plot(X, py + shifted * ph, color=RED, lw=1.85, ls=(0, (4, 2)))
    text(px - 0.08, py + ph / 2, r"$P_{\rm sw}$", fs=10.3, ha="right")
    text(px + pw / 2, py - 0.12, r"$u$", fs=10.3)
    arrow(px + 0.30 * pw, py + 0.45 * ph, px + 0.55 * pw, py + 0.45 * ph,
          color=RED, lw=1.20, ms=7.5)
    text(px + 0.26 * pw, py + 0.34 * ph, r"$u_{\rm th}$ shift",
         fs=9.0, color=RED, ha="left", style="italic")
    ax.annotate("",
                xy=(px + 0.84 * pw, py + 0.26 * ph),
                xytext=(px + 0.84 * pw, py + 0.86 * ph),
                arrowprops=dict(arrowstyle="<->", color=RED, lw=1.20))
    text(px + 0.86 * pw, py + 0.56 * ph, r"$\beta_s$" + "\n" + "broadening",
         fs=9.0, color=RED, ha="left", style="italic")


W, H = 4.72, 2.62
px1, px2 = 0.18, 5.10
py_top, py_bot = 3.04, 0.16

pa = panel(px1, py_top, W, H, "(a)", "Self-heating")
pb = panel(px2, py_top, W, H, "(b)", "T-dependent channels")
pc = panel(px1, py_bot, W, H, "(c)", "Stochastic LLG")
pd = panel(px2, py_bot, W, H, "(d)", r"$P_{\rm sw}$ response")

draw_heat(*pa)
draw_channels(*pb)
draw_llg(*pc)
draw_psw(*pd)

arrow(2.55, 3.01, 2.55, 2.79, color=AMBER, lw=0.88, ms=6.5)
arrow(7.45, 3.01, 7.45, 2.79, color=TEAL, lw=0.88, ms=6.5)
arrow(4.90, 1.47, 5.08, 1.47, color=PURPLE, lw=0.92, ms=7.0)


out_path = Path(__file__).resolve().parent / "Chapter02_local_03.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white",
            pad_inches=0.04)
plt.close(fig)

chapter_fig = (
    Path(__file__).resolve().parent.parent
    / "article" / "00_chapter_drafts" / "figs" / "Chapter02_local_03.png"
)
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)

print(f"Saved  {out_path}")
print(f"Synced {chapter_fig}")
