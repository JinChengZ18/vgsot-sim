"""
3D 单位球面上的完整磁化矢量 m(t) 轨迹绘制。

利用 SimResult 已经记录的 theta(t) 与 phi(t)，重构
  m(t) = (sin theta cos phi, sin theta sin phi, cos theta)
并在单位球面上画出轨迹，配以时间彩色梯度（清华紫调色）、起止点标记，以及
m_x, m_y, m_z 三分量随时间的演化。明亮的现代风格 + 清华紫主题。

输出 Chapter02_local_10.png（本目录 + 同步到 article/00_chapter_drafts/figs/）。
"""
from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation

# ── Palette (THU purple family for accents, cool-warm for trajectory) ──
THU_DEEP, THU_MID, THU_SOFT, THU_PALE = "#660874", "#8B3A9E", "#A966BE", "#C99FD4"
THU_TINT, THU_LIGHT = "#EFE0F7", "#F7ECFB"
THU_GRID = "#DDD0E8"
CRIMSON, NAVY, TEAL, AMBER = "#A82038", "#1F5FA8", "#1A6B5A", "#C47A00"
CHARCOAL, NEAR_WHITE = "#2B2B2B", "#FFFFFF"
SPHERE_LINE = "#E8E0EE"   # very light lavender for sphere wireframe

# Use matplotlib's coolwarm (NAVY-blue → CRIMSON-red diverging) so the
# time-encoded trajectory matches the cool-warm palette used in
# Chapter02_local_11 and naturally maps "cold start at AP" → "warm end at P".
TRAJ_CMAP = plt.get_cmap("coolwarm")

plt.rcParams.update({
    "font.family"      : "sans-serif",
    "font.sans-serif"  : ["Arial", "Liberation Sans"],
    "font.size"        : 12,
    "axes.labelsize"   : 13,
    "axes.titlesize"   : 14,
    "axes.titlepad"    : 8,
    "legend.fontsize"  : 11,
    "xtick.labelsize"  : 11,
    "ytick.labelsize"  : 11,
    "mathtext.fontset" : "stix",
    "axes.linewidth"   : 0.9,
    "lines.linewidth"  : 1.7,
    "figure.dpi"       : 150,
    "savefig.dpi"      : 300,
})


# ── Simulation ──────────────────────────────────────────────────────────
cc = PhysicalConstantsConfig()
PULSE_NS  = 0.75
RELAX_NS  = 3.25
I_SOT_UA  = -2000
T_AMBIENT = 300.0

sim_end = int(round((PULSE_NS + RELAX_NS) * 1e-9 / cc.t_step))
mid1    = int(round(PULSE_NS * 1e-9 / cc.t_step))

np.random.seed(1)
res = run_piecewise_direct_excitation(
    sim_start_step=1, sim_mid1_step=mid1, sim_mid2_step=sim_end, sim_end_step=sim_end,
    pap=1,
    v_mtj_stage1=0.0, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
    i_sot_stage1=I_SOT_UA * 1e-6, i_sot_stage2=0.0, i_sot_stage3=0.0,
    estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1, estt_stage3=0, esot_stage3=1,
    vnv=0, non=1, r_sot_fl_dl=0.83,
    constants=cc, show_progress=False,
    enable_self_heating=True, T_ambient_K=T_AMBIENT,
)

t_ns  = res.time_s * 1e9
theta = res.theta
phi   = res.phi
mx    = np.sin(theta) * np.cos(phi)
my    = np.sin(theta) * np.sin(phi)
mz    = np.cos(theta)
norm  = np.sqrt(mx**2 + my**2 + mz**2)
mx, my, mz = mx / norm, my / norm, mz / norm


# ── Figure: balanced 3D sphere and time-series panels ──────────────────
fig = plt.figure(figsize=(12.0, 6.0), facecolor=NEAR_WHITE)
# Sphere panel marginally wider than time-series; gap tightened so the two
# panels feel visually anchored rather than separated.
gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0],
                      wspace=0.16, left=0.04, right=0.97,
                      top=0.94, bottom=0.10)

# ── Panel A — 3D sphere with trajectory (DOMINANT) ──────────────────────
ax3d = fig.add_subplot(gs[0, 0], projection="3d")
ax3d.set_facecolor(NEAR_WHITE)

# Unit sphere — wireframe ONLY (no surface), fully transparent interior
u = np.linspace(0, 2 * np.pi, 80)
v = np.linspace(0,     np.pi, 40)
xs = np.outer(np.cos(u), np.sin(v))
ys = np.outer(np.sin(u), np.sin(v))
zs = np.outer(np.ones_like(u), np.cos(v))

# Faint wireframe only — the inside is transparent so trajectory segments
# behind the sphere remain visible through it. Slightly stronger than the
# previous extra-faint setting since there is no surface backing it.
ax3d.plot_wireframe(xs, ys, zs, color=SPHERE_LINE, alpha=0.70,
                     rstride=6, cstride=12, linewidth=0.55)

# Equator highlight (thin reference circle)
phi_eq = np.linspace(0, 2 * np.pi, 200)
ax3d.plot(np.cos(phi_eq), np.sin(phi_eq), np.zeros_like(phi_eq),
          color=SPHERE_LINE, lw=0.9, alpha=0.95)

# Coordinate axes through origin
for direction, lab, col in [((1, 0, 0), r"$x$", CHARCOAL),
                             ((0, 1, 0), r"$y$", CHARCOAL),
                             ((0, 0, 1), r"$z$", CHARCOAL)]:
    ax3d.plot([0, 1.25 * direction[0]],
              [0, 1.25 * direction[1]],
              [0, 1.25 * direction[2]],
              color=col, lw=1.0)
    ax3d.text(1.40 * direction[0], 1.40 * direction[1], 1.40 * direction[2],
              lab, fontsize=13, fontweight="bold", color=col,
              ha="center", va="center")

# Two stable poles — colour-matched to coolwarm endpoints
ax3d.scatter([0], [0], [+1], s=110, color=TRAJ_CMAP(1.0), edgecolor=NEAR_WHITE,
             linewidth=1.5, zorder=4)
ax3d.scatter([0], [0], [-1], s=110, color=TRAJ_CMAP(0.0), edgecolor=NEAR_WHITE,
             linewidth=1.5, zorder=4)
ax3d.text(0.0, 0.18, 1.05, r"P ($m_z=+1$)", fontsize=10.5,
          color=CRIMSON, ha="left", va="bottom", fontweight="bold")
ax3d.text(0.0, 0.18, -1.10, r"AP ($m_z=-1$)", fontsize=10.5,
          color=NAVY, ha="left", va="top", fontweight="bold")

# Trajectory (time-coloured segments)
seg = np.linspace(0, 1, len(t_ns))
for i in range(len(t_ns) - 1):
    ax3d.plot([mx[i], mx[i + 1]],
              [my[i], my[i + 1]],
              [mz[i], mz[i + 1]],
              color=TRAJ_CMAP(seg[i]), lw=2.0, alpha=0.95)

# Start / end markers
ax3d.scatter([mx[0]],  [my[0]],  [mz[0]],  s=120, color=TRAJ_CMAP(0.0),
             edgecolor=CHARCOAL, linewidth=1.5, zorder=5)
ax3d.scatter([mx[-1]], [my[-1]], [mz[-1]], s=180, color=TRAJ_CMAP(1.0),
             edgecolor=CHARCOAL, linewidth=1.5, marker="*", zorder=6)

# Sphere occupies the full panel
ax3d.set_xlim(-1.05, 1.05)
ax3d.set_ylim(-1.05, 1.05)
ax3d.set_zlim(-1.05, 1.05)
ax3d.set_box_aspect((1, 1, 1))
ax3d.view_init(elev=22, azim=42)
ax3d.set_axis_off()

# Hide 3D pane backgrounds for cleaner modern look
ax3d.xaxis.pane.fill = False
ax3d.yaxis.pane.fill = False
ax3d.zaxis.pane.fill = False
ax3d.xaxis.pane.set_edgecolor("none")
ax3d.yaxis.pane.set_edgecolor("none")
ax3d.zaxis.pane.set_edgecolor("none")

ax3d.set_title(r"Magnetisation trajectory $\mathbf{m}(t)$ on the unit sphere",
               fontsize=13.5, fontweight="bold", color=CHARCOAL, pad=4)

# Colorbar placed INSIDE the sphere panel (left-bottom corner) so it is
# attached to the 3D plot rather than floating in the gap between panels.
cbax = ax3d.inset_axes([-0.02, 0.10, 0.028, 0.45])
mappable = cm.ScalarMappable(cmap=TRAJ_CMAP,
                             norm=plt.Normalize(vmin=t_ns.min(),
                                                 vmax=t_ns.max()))
cb = fig.colorbar(mappable, cax=cbax)
cb.set_label("time (ns)", fontsize=10.5, color=CHARCOAL, labelpad=4)
cb.ax.tick_params(labelsize=9.5, color=CHARCOAL)
cb.outline.set_edgecolor(CHARCOAL)
cb.outline.set_linewidth(0.5)


# ── Panel B — compact time-series components ───────────────────────────
ax2d = fig.add_subplot(gs[0, 1])
ax2d.set_facecolor(NEAR_WHITE)
ax2d.plot(t_ns, mx, color=NAVY,    lw=1.5, alpha=0.78, label=r"$m_x(t)$")
ax2d.plot(t_ns, my, color=TEAL,    lw=1.5, alpha=0.78, label=r"$m_y(t)$")
ax2d.plot(t_ns, mz, color=CRIMSON, lw=2.0, label=r"$m_z(t)$")
ax2d.axvline(PULSE_NS, color=CRIMSON, lw=1.2, ls=(0, (4, 2)))
ax2d.text(PULSE_NS + 0.04, 1.05, "pulse off",
          color=CRIMSON, fontsize=11, fontweight="bold",
          ha="left", va="top")
ax2d.axhline(0, color="gray", lw=0.5, ls=":")
ax2d.set_xlim(0, (PULSE_NS + RELAX_NS))
ax2d.set_ylim(-1.15, 1.15)
ax2d.set_xlabel("Time (ns)", fontsize=11.5)
ax2d.set_ylabel("magnetisation component", fontsize=11.5)
ax2d.set_title(rf"Cartesian components, $I_{{\rm SOT}} = {I_SOT_UA}\,\mu$A",
               fontsize=12, fontweight="bold", color=CHARCOAL, pad=6)
ax2d.grid(True, linestyle="--", linewidth=0.45, color=THU_GRID, alpha=0.7)
ax2d.legend(loc="center right", frameon=True,
            edgecolor=THU_PALE, framealpha=0.92)

for spine in ax2d.spines.values():
    spine.set_color(CHARCOAL)
    spine.set_linewidth(0.8)
ax2d.spines["top"].set_visible(False)
ax2d.spines["right"].set_visible(False)


# ── Save + sync ─────────────────────────────────────────────────────────
out_path = Path(__file__).resolve().parent / "Chapter02_local_10.png"
fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor=NEAR_WHITE,
            pad_inches=0.20)
plt.close(fig)
chapter_fig = (Path(__file__).resolve().parent.parent.parent /
               "article" / "figs" / "Chapter02_local_10.png")
chapter_fig.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(out_path, chapter_fig)
print(f"Saved  {out_path}")
print(f"Synced {chapter_fig}")
