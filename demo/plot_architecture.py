"""
Figure 2.8: vgsot-sim software architecture.

Journal-style revision: portrait-oriented layered architecture, reduced
verbal load, manuscript-scale typography, and orthogonal data flow.

Outputs:
  demo/Chapter02_local_08.png
  article/figs/Chapter02_local_08.png
"""
from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


PURPLE = "#5A136F"
PURPLE_TINT = "#F1E8F5"
BLUE = "#1F5EA8"
TEAL = "#087566"
AMBER = "#B87900"
INK = "#202020"
MID = "#666666"
GRID = "#D8D3DC"
LINE = "#2D2D2D"
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


fig, ax = plt.subplots(figsize=(6.25, 6.05))
fig.subplots_adjust(left=0.025, right=0.985, bottom=0.025, top=0.985)
ax.set_xlim(0, 7.0)
ax.set_ylim(0, 9.45)
ax.axis("off")


def text(x: float, y: float, s: str, *, fs: float = 10.0,
         color: str = INK, weight: str = "normal",
         ha: str = "center", va: str = "center",
         style: str = "normal") -> None:
    ax.text(x, y, s, fontsize=fs, color=color, fontweight=weight,
            ha=ha, va=va, fontstyle=style)


def arrow(x0: float, y0: float, x1: float, y1: float, *,
          color: str = INK, lw: float = 0.95, ms: float = 7.5,
          rad: float = 0.0) -> None:
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle="-|>", mutation_scale=ms, lw=lw, color=color,
        shrinkA=1.0, shrinkB=1.0, connectionstyle=f"arc3,rad={rad}",
    ))


def polyline(points: list[tuple[float, float]], *, color: str = INK,
             lw: float = 0.95, ms: float = 7.5) -> None:
    for (x0, y0), (x1, y1) in zip(points[:-2], points[1:-1]):
        ax.plot([x0, x1], [y0, y1], color=color, lw=lw,
                solid_capstyle="round")
    x0, y0 = points[-2]
    x1, y1 = points[-1]
    arrow(x0, y0, x1, y1, color=color, lw=lw, ms=ms)


def layer(x: float, y: float, w: float, h: float,
          title: str, note: str, fc: str) -> None:
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc,
                           edgecolor=LINE, lw=0.78))
    text(x + 0.18, y + h - 0.22, title, fs=12.8,
         weight="bold", ha="left")
    text(x + w - 0.18, y + h - 0.22, note, fs=9.6,
         color=MID, ha="right", style="italic")
    ax.plot([x + 0.16, x + w - 0.16], [y + h - 0.46, y + h - 0.46],
            color=GRID, lw=0.62)


def box(x: float, y: float, w: float, h: float, title: str,
        subtitle: str = "", *, fc: str = WHITE, ec: str = LINE,
        fs: float = 10.0, sub_fs: float = 8.6,
        lw: float = 0.78) -> None:
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.018,rounding_size=0.035",
        facecolor=fc, edgecolor=ec, lw=lw,
    ))
    if subtitle:
        text(x + w / 2, y + h * 0.64, title, fs=fs, weight="bold")
        text(x + w / 2, y + h * 0.30, subtitle, fs=sub_fs,
             style="italic")
    else:
        text(x + w / 2, y + h / 2, title, fs=fs, weight="bold")


def case_chip(x: float, y: float, w: float, h: float, label: str) -> None:
    ax.add_patch(Rectangle((x, y), w, h, facecolor=WHITE,
                           edgecolor=GRID, lw=0.62))
    text(x + w / 2, y + h / 2, label, fs=8.7)


LX, LW = 0.32, 6.36

# Layer 1: user-facing entry points.
layer(LX, 7.58, LW, 1.52, "User interfaces and I/O",
      "CLI / Python API", PURPLE_TINT)
box(0.58, 7.85, 1.24, 0.54, "install", "pip", fs=10.0, sub_fs=8.5)
box(2.18, 7.85, 1.74, 0.54, "CLI", "sweeps", fs=10.7,
    sub_fs=8.8, ec=PURPLE)
box(4.38, 7.85, 2.02, 0.54, "Python API", "scripts", fs=10.7,
    sub_fs=8.8, ec=PURPLE)

# Layer 2: cases and configuration.
layer(LX, 4.72, LW, 2.22, "Experiment scenarios and configuration",
      "cases / configs", "#F7F1FA")
ax.add_patch(FancyBboxPatch(
    (0.58, 5.03), 3.15, 1.15,
    boxstyle="round,pad=0.018,rounding_size=0.035",
    facecolor=WHITE, edgecolor=LINE, lw=0.78))
text(2.15, 5.92, "standard cases", fs=10.2, weight="bold")
case_labels = [
    "SOT baseline", "VCMA-assisted",
    "optimized VGSOT", "SER Monte Carlo",
]
for i, label in enumerate(case_labels):
    col = i % 2
    row = i // 2
    case_chip(0.78 + col * 1.39, 5.17 + (1 - row) * 0.32,
              1.18, 0.27, label)

box(4.18, 5.10, 2.22, 0.98, "configuration",
    r"$J,\ V,\ t_w,\ T$", fs=10.4, sub_fs=9.0, ec=TEAL)
arrow(3.76, 5.58, 4.14, 5.58, color=MID, lw=0.82, ms=7.0)

# Layer 3: core kernels.
layer(LX, 0.48, LW, 3.54, "Physical kernels and core modules",
      "single-step update loop", "#F3E7F7")
box(0.66, 2.76, 5.72, 0.66, "dynamic_switching",
    r"$\mathbf{m}_{n+1}=\mathcal{F}(\mathbf{m}_n,\mathbf{H}_{\rm eff},\mathbf{H}_{\rm th})$",
    fs=10.4, sub_fs=8.5, ec=PURPLE)

module_specs = [
    ("anisotropy", "field", BLUE),
    ("stochastic", "noise", AMBER),
    ("electronic", "transport", TEAL),
    ("tmr", "resistance", PURPLE),
]
mod_xs = [0.66, 2.10, 3.54, 4.98]
for x0, (title, sub, col) in zip(mod_xs, module_specs):
    box(x0, 1.50, 1.18, 0.82, title, sub, fs=9.3, sub_fs=8.0, ec=col)

# Top-to-bottom control flow.
polyline([(3.05, 7.84), (3.05, 7.22), (0.20, 7.22),
          (0.20, 5.58), (0.58, 5.58)],
         color=PURPLE, lw=0.90, ms=7.5)
polyline([(5.38, 7.84), (5.38, 7.22), (6.82, 7.22),
          (6.82, 5.59), (6.42, 5.59)],
         color=PURPLE, lw=0.90, ms=7.5)
polyline([(5.30, 5.10), (5.30, 4.36), (6.82, 4.36),
          (6.82, 3.09), (6.40, 3.09)],
         color=LINE, lw=0.86, ms=7.0)

# Kernel chain and feedback loop.
for x0, x1 in [(1.84, 2.08), (3.28, 3.52), (4.72, 4.96)]:
    arrow(x0, 1.91, x1, 1.91, color=PURPLE, lw=0.90, ms=7.0)
arrow(3.52, 2.76, 3.52, 2.34, color=LINE, lw=0.86, ms=7.0)
text(3.42, 1.05, "single-step data path",
     fs=8.9, color=PURPLE, style="italic")

# Inputs and outputs are kept inside the layer to avoid edge crowding.
box(0.66, 0.66, 1.18, 0.46, "inputs", r"$J,V,t_w,T$", fs=9.0,
    sub_fs=7.9, ec=LINE, fc="#FAFAFA")
box(4.98, 0.66, 1.40, 0.46, "outputs", r"$P_{\rm sw}$ / SER",
    fs=9.0, sub_fs=7.9, ec=LINE, fc="#FAFAFA")
arrow(1.25, 1.14, 1.25, 1.56, color=LINE, lw=0.82, ms=6.8)
arrow(5.57, 1.50, 5.57, 1.14, color=LINE, lw=0.82, ms=6.8)


out_path = Path(__file__).resolve().parent / "Chapter02_local_08.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white",
            pad_inches=0.05)
plt.close(fig)

chapter_fig = (
    Path(__file__).resolve().parent.parent
    / "article" / "figs" / "Chapter02_local_08.png"
)
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)

print(f"Saved  {out_path}")
print(f"Synced {chapter_fig}")
