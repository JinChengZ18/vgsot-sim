"""Journal-grade schematic panels for the norm-preserving Cayley integrator
section (explicit-omega step, convergence orders, Wong-Zakai drift, T_eff).

Four CLEAN panels (English-only, Arial, no baked-in panel letters or figure
numbers -- those are added later in the compositing PPT), plus one 2x2
preview composite for quick inspection only:

  integrator_panel_geometry.png   schematic: tangent Euler leaves the unit
                                  sphere, the Cayley step is an exact rotation;
                                  omega evaluation point (left endpoint vs
                                  midpoint) sets the deterministic order
  integrator_panel_order.png      MEASURED deterministic global convergence
                                  (data: result/sec_2_2_3_2/A/*.json)
  integrator_panel_drift.png      analytic stationary laws with/without the
                                  Wong-Zakai drift D*cot(theta) (annotated with
                                  the measured T_eff/T anchors)
  integrator_panel_teff.png       MEASURED equilibrium effective temperature,
                                  three integrators + dt->0 extrapolation
                                  (data: result/sec_2_2_3_2/B/teff_full_merged.json)
  integrator_panels_preview.png   2x2 preview (NOT for publication)

Every data point in the two measured panels comes from the committed,
reproducible experiment scripts (order_of_accuracy.py, teff_dtsweep.py +
merge_full_results.py); nothing is estimated or drawn by hand there.

Run:  PYTHONPATH=src python scripts/09_simulation_figures/fig_integrator_verification_panels.py
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "integrator_panels")
os.makedirs(OUT, exist_ok=True)

import sys
sys.path.insert(0, os.path.join(ROOT, "scripts", "11_boltzmann_teff"))


def _drift_hist_samples():
    """Equilibrium theta samples for the drift panel, from the SAME engine as
    experiment B (teff_dtsweep.run_trajectory), seeded and cached to JSON so
    re-renders are instant and the histogram is reproducible."""
    cache = os.path.join(OUT, "drift_hist_cache.json")
    if os.path.exists(cache):
        with open(cache, encoding="utf-8") as f:
            return json.load(f)

    import dataclasses
    import teff_dtsweep as tds  # scripts/11_boltzmann_teff (engine of exp B)
    from vgsot_sim.configs import PhysicalConstantsConfig
    from vgsot_sim.demag import demag_factors

    cc = PhysicalConstantsConfig()
    cc.h_ex_x = cc.h_ex_y = cc.h_ex_z = 0.0
    dt = 0.2e-12
    cc = dataclasses.replace(cc, t_step=dt)
    Nx, Ny, Nz = demag_factors(cc, mode="ellipsoid")
    Keff = cc.Ki / cc.tf - 0.5 * cc.u0 * cc.Ms ** 2 * (Nz - Nx)
    Delta = Keff * (cc.tf * cc.A1) / (cc.kb * 300.0)
    eq = tds._make_eq_sampler(Delta)
    n_steps = int(round(12e-9 / dt))
    n_burn = n_steps // 4
    out = {"Delta": Delta, "dt_ps": 0.2, "T_phys_ns": 12.0, "n_traj": 16}
    for integ, key in (("cayley_explicit", "cayley"),
                       ("euler_spherical", "euler")):
        xs = []
        for k in range(16):
            rng = np.random.default_rng(20260705 + 131 * k)
            mz = tds.run_trajectory(cc, dt, n_steps, eq(rng), rng, integ)
            th = np.arccos(np.clip(np.abs(mz[n_burn:]), -1.0, 1.0))
            xs.append(th[::25] * np.sqrt(Delta))   # thin: ~25 tau_int/traj
        out[key] = np.concatenate(xs).tolist()
    with open(cache, "w", encoding="utf-8") as f:
        json.dump(out, f)
    return out

# ── palette (scheme colours stay consistent across panels) ───────────────
PURPLE = "#660874"   # explicit-omega Cayley (published kernel)
TEAL   = "#1A6B5A"   # true implicit midpoint (cayley_midpoint)
NAVY   = "#1F5FA8"   # spherical-coordinate Euler
CRIMSON = "#A82038"  # defective / off-sphere behaviour
CHARCOAL = "#2B2B2B"
GRAY = "#8A8A8A"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12.5,
    "axes.titlepad": 8,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "mathtext.fontset": "stix",
    "axes.linewidth": 0.9,
    "axes.edgecolor": CHARCOAL,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "lines.linewidth": 1.7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

PANEL_SIZE = (4.9, 3.9)


# ══════════════════════════════════════════════════════════════════════════
# Panel 1 — geometry: tangent Euler vs Cayley rotation on the unit sphere
# ══════════════════════════════════════════════════════════════════════════
def panel_geometry(ax):
    # great-circle cross-section of the unit sphere
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color=GRAY, lw=1.1, zorder=1)
    ax.annotate(r"unit sphere  $|\mathbf{m}| = 1$",
                xy=(np.cos(3.75), np.sin(3.75)), xytext=(-1.83, -0.72),
                fontsize=9.5, color=GRAY, va="top",
                arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.8))
    # define the rotation form the whole figure is built on
    ax.text(-1.83, -0.16,
            "rotation form:\n"
            r"$\dot{\mathbf{m}} = \mathbf{w}\times\mathbf{m}$",
            fontsize=10, color=CHARCOAL, va="top")

    a0 = -0.10                     # angular position of m_n (right, near equator)
    dphi = 0.90                    # rotation angle of one (exaggerated) step
    m0 = np.array([np.cos(a0), np.sin(a0)])
    m1 = np.array([np.cos(a0 + dphi), np.sin(a0 + dphi)])

    # state vectors from the origin
    ax.annotate("", xy=m0, xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=CHARCOAL, lw=1.6))
    ax.annotate("", xy=m1, xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=PURPLE, lw=1.6))
    ax.text(m0[0] + 0.07, m0[1] - 0.20, r"$\mathbf{m}_n$", fontsize=13,
            color=CHARCOAL)
    ax.text(m1[0] - 0.34, m1[1] + 0.10, r"$\mathbf{m}_{n+1}$", fontsize=13,
            color=PURPLE)

    # naive tangent step (geometry illustration; distinct from the
    # (theta,phi)-chart Euler scheme of the stochastic panels)
    tang = np.array([-np.sin(a0), np.cos(a0)])
    me = m0 + dphi * tang
    ax.annotate("", xy=me, xytext=m0,
                arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.4,
                                linestyle="--", alpha=0.9))
    ax.plot(*me, marker="x", ms=9, mew=2.2, color=CRIMSON)
    ax.text(1.30, 1.44,
            "explicit tangent step\n"
            r"$\mathbf{m}_n + \Delta t\,(\mathbf{w}\times\mathbf{m}_n)$"
            "\n"
            r"$|\mathbf{m}| = \sqrt{1+\Delta t^2 |\mathbf{w}_{\!\perp}|^2} > 1$"
            "\n(norm drifts off the sphere)",
            fontsize=9.5, color=CRIMSON, ha="left", va="top")

    # Cayley rotation arc (on the circle)
    arc = np.linspace(a0, a0 + dphi, 60)
    ax.plot(np.cos(arc), np.sin(arc), color=PURPLE, lw=2.6, zorder=3)
    amid_arrow = a0 + 0.80 * dphi
    ax.annotate("",
                xy=(np.cos(amid_arrow + 0.02), np.sin(amid_arrow + 0.02)),
                xytext=(np.cos(amid_arrow - 0.02), np.sin(amid_arrow - 0.02)),
                arrowprops=dict(arrowstyle="-|>", color=PURPLE, lw=2.2))
    ax.text(-1.80, 1.58,
            "Cayley step = rotation through $\\Delta t\\,\\mathbf{w}$:\n"
            r"$\mathbf{m}_{n+1} = \mathbf{A}_n\mathbf{m}_n$"
            "\n"
            r"$|\mathbf{m}| = 1$ exactly, any $\Delta t$",
            fontsize=9.5, color=PURPLE, ha="left", va="top")

    # w evaluation points: left endpoint vs midpoint of the step
    amid = a0 + dphi / 2.0
    mm = np.array([np.cos(amid), np.sin(amid)])
    ax.plot(*m0, marker="o", ms=7, color=CHARCOAL, zorder=4)
    ax.plot(*mm, marker="o", ms=8, mfc="white", mec=TEAL, mew=2.0, zorder=5)
    ax.annotate(r"left-endpoint $\mathbf{w}(\mathbf{m}_n)$:  $O(\Delta t)$",
                xy=m0 + np.array([0.02, -0.02]), xytext=(0.05, -0.54),
                fontsize=9.5, color=CHARCOAL, ha="left",
                arrowprops=dict(arrowstyle="-", color=CHARCOAL, lw=0.8))
    ax.annotate(r"midpoint $\mathbf{w}$:  $O(\Delta t^2)$",
                xy=mm + np.array([0.02, -0.03]), xytext=(1.30, 0.16),
                fontsize=9.5, color=TEAL, ha="left",
                arrowprops=dict(arrowstyle="-", color=TEAL, lw=0.9,
                                connectionstyle="arc3,rad=0.15"))

    ax.set_xlim(-1.90, 3.15)
    ax.set_ylim(-1.24, 1.60)
    ax.set_aspect("equal", adjustable="box", anchor="N")  # top-anchored → title aligns
    ax.axis("off")
    ax.set_title("Geometric step and $\\mathbf{w}$ evaluation point")


# ══════════════════════════════════════════════════════════════════════════
# Panel 2 — measured deterministic global convergence (experiment A data)
# ══════════════════════════════════════════════════════════════════════════
def panel_order(ax):
    with open(os.path.join(ROOT, "result", "sec_2_2_3_2", "A",
                           "order_of_accuracy_results.json"),
              encoding="utf-8") as f:
        A = json.load(f)

    # label placement: dy stagger at the coarse (right) end so the two
    # first-order curves cannot collide
    series = [
        ("euler_spherical", NAVY, "^",
         "$(\\theta,\\phi)$-chart Euler", 1.0),
        ("published", PURPLE, "o",
         "Cayley, left-endpoint $\\mathbf{w}$", 0.20),
        ("midpoint", TEAL, "s", "Cayley, midpoint $\\mathbf{w}$", 1.0),
    ]
    for key, color, marker, label, dy in series:
        c = A["curves"][key]
        dt_ps = np.asarray(c["dts"], float) * 1e12
        err = np.asarray(c["errs"], float)
        order = np.argsort(dt_ps)
        dt_ps, err = dt_ps[order], err[order]
        p = c["p_asymptotic"]
        ax.loglog(dt_ps, err, marker=marker, ms=5.5, color=color,
                  mfc="white", mew=1.4, lw=1.6)
        # direct single-line label at the coarse end of each curve (no legend)
        ax.annotate(f"{label},  $p = {p:.2f}$",
                    xy=(dt_ps[-1], err[-1]),
                    xytext=(dt_ps[-1] * 1.3, err[-1] * dy),
                    fontsize=9.5, color=color, va="center")

    # slope guide triangles in guaranteed-free zones
    def guide(x0, x1, y0, power, label_dx=1.06):
        xx = np.array([x0, x1])
        yy = y0 * (xx / x0) ** power
        ax.plot([x0, x1, x1, x0], [y0, y0, yy[1], y0], color=GRAY, lw=0.9)
        ax.text(x1 * label_dx, np.sqrt(y0 * yy[1]), f"{power}",
                fontsize=9, color=GRAY, va="center")
    guide(0.35, 1.05, 2.2e-4, 1)
    guide(5.0, 15.0, 6.0e-7, 2)

    # bridge to the stochastic panels: the chart-Euler defect is NOT here
    ax.text(0.037, 2.6e-6,
            "$(\\theta,\\phi)$-chart Euler converges normally here —\n"
            "its defect appears only in equilibrium sampling",
            fontsize=9, color=NAVY, va="bottom")
    ax.text(0.037, 0.42, "fitted $p$: 6 finest $\\Delta t$ (asymptotic range)",
            fontsize=8.5, color=GRAY, va="top")

    ax.set_xlabel(r"$\Delta t$ (ps)")
    ax.set_ylabel(r"global error  $\left|\mathbf{m}(T) - \mathbf{m}_{\mathrm{ref}}(T)\right|$")
    ax.set_title("Deterministic global convergence (measured)")
    ax.set_xlim(0.03, 170)
    ax.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)


# ══════════════════════════════════════════════════════════════════════════
# Panel 3 — the missing Wong-Zakai drift and its stationary-law consequence
# ══════════════════════════════════════════════════════════════════════════
def panel_drift(ax):
    # dimensionless polar angle x = theta*sqrt(Delta); small-angle analytic
    # laws (sin(theta) ~ theta holds: the well confines theta to ~1/sqrt(48.5))
    x = np.linspace(0, 3.5, 400)
    p_ok = 2.0 * x * np.exp(-x ** 2)                    # Rayleigh (correct)
    p_bad = 2.0 / np.sqrt(np.pi) * np.exp(-x ** 2)      # Jacobian lost

    # measured equilibrium samples from the experiment-B engine (cached)
    H = _drift_hist_samples()
    bins = np.linspace(0, 3.5, 46)
    for key, color in (("cayley", PURPLE), ("euler", NAVY)):
        ax.hist(np.asarray(H[key], float), bins=bins, density=True,
                histtype="step", lw=1.3, color=color, alpha=0.85)

    ax.plot(x, p_ok, color=CHARCOAL, lw=2.0)
    ax.plot(x, p_bad, color=NAVY, lw=2.0, ls="--")

    ax.annotate("with noise-induced drift $D\\cot\\theta$:\n"
                "$p(x)\\propto x\\,e^{-x^2}$  (Boltzmann,"
                " $\\sin\\theta\\!\\approx\\!\\theta$)\n"
                "sampled by both Cayley steps:"
                " $T_{\\mathrm{eff}}/T = 0.979$",
                xy=(1.10, 0.61), xytext=(1.62, 0.50),
                fontsize=9.5, color=CHARCOAL, va="top", ha="left",
                arrowprops=dict(arrowstyle="-", color=CHARCOAL, lw=0.8))
    ax.annotate("drift term missing in the $(\\theta,\\phi)$-chart Euler:\n"
                "Jacobian lost, $p(x)\\propto e^{-x^2}$\n"
                "measured: $T_{\\mathrm{eff}}/T = 0.492 \\approx 1/2$",
                xy=(0.30, 1.03), xytext=(1.18, 1.46),
                fontsize=9.5, color=NAVY, va="top", ha="left",
                arrowprops=dict(arrowstyle="-", color=NAVY, lw=0.9))
    ax.text(3.42, 0.90,
            "histograms: sampled equilibrium\n"
            "($\\Delta t = 0.2$ ps, 12 ns, 16 trajectories)",
            fontsize=8.5, color=GRAY, ha="right", va="top")

    ax.set_xlabel(r"scaled polar angle  $x = \theta\sqrt{\Delta}$"
                  r"   ($\Delta = E_b/k_BT = 48.5$)")
    ax.set_ylabel(r"stationary density  $p(x)$")
    ax.set_title("Stochastic drift in the $(\\theta,\\phi)$ chart")
    ax.set_xlim(0, 3.5)
    ax.set_ylim(0, 1.55)
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)


# ══════════════════════════════════════════════════════════════════════════
# Panel 4 — measured equilibrium effective temperature (experiment B data)
# ══════════════════════════════════════════════════════════════════════════
def panel_teff(ax):
    with open(os.path.join(ROOT, "result", "sec_2_2_3_2", "B",
                           "teff_full_merged.json"), encoding="utf-8") as f:
        B = json.load(f)

    series = [
        ("cayley_explicit", PURPLE, "o",
         "Cayley, left-endpoint $\\mathbf{w}$"),
        ("cayley_true_midpoint", TEAL, "s", "Cayley, midpoint $\\mathbf{w}$"),
        ("euler_spherical", NAVY, "^", "$(\\theta,\\phi)$-chart Euler"),
    ]
    label_y = {"cayley_explicit": 1.075, "cayley_true_midpoint": 0.885,
               "euler_spherical": 0.60}
    star_ms = {"cayley_explicit": 13, "cayley_true_midpoint": 10,
               "euler_spherical": 13}
    # de-collide the two near-identical Cayley series (points AND dt->0
    # stars coincide): small horizontal jitter + distinct dot patterns
    x_off = {"cayley_explicit": -0.004, "cayley_true_midpoint": +0.004,
             "euler_spherical": 0.0}
    fit_ls = {"cayley_explicit": (0, (1, 1.6)),
              "cayley_true_midpoint": (0, (3, 1.6)),
              "euler_spherical": (0, (1, 1.6))}
    for key, color, marker, label in series:
        rows = sorted(B["sweep"][key], key=lambda r: r["dt_ps"])
        dt = np.array([r["dt_ps"] for r in rows]) + x_off[key]
        y = np.array([r["moment"]["point"] for r in rows])
        e = np.array([r["moment"]["sigma"] for r in rows])
        # markers + error bars only (no connecting line: points are
        # independent runs, the dotted line is the weighted linear fit)
        ax.errorbar(dt, y, yerr=e, marker=marker, ms=5.5, color=color,
                    mfc="white", mew=1.4, ls="none", capsize=2.5)
        ex = B["extrapolation"][key]["moment"]
        xx = np.linspace(0, 0.42, 30)
        ax.plot(xx, ex["intercept"] + ex["slope"] * xx, ls=fit_ls[key],
                lw=1.1, color=color, alpha=0.8)
        ax.plot(x_off[key], ex["intercept"], marker="*", ms=star_ms[key],
                color=color, mec=CHARCOAL, mew=0.5, zorder=5, clip_on=False)
        ax.text(0.165, label_y[key],
                f"{label},  $\\Delta t\\to 0$: "
                f"${ex['intercept']:.4f}\\pm{ex['intercept_sigma']:.4f}$",
                fontsize=9.5, color=color, va="center")

    ax.axhline(1.0, color=GRAY, lw=1.0)
    ax.text(0.008, 1.014, "Boltzmann", fontsize=9, color=GRAY, ha="left")
    ax.axhline(0.5, color=GRAY, lw=1.0, ls="--")
    ax.text(0.008, 0.522, r"$T/2$", fontsize=9, color=GRAY, ha="left")
    ax.text(0.008, 0.435,
            "error bars: block-bootstrap $1\\sigma$ (48 trajectories);"
            "  $\\star$ = weighted $\\Delta t\\to 0$ fit",
            fontsize=8.5, color=GRAY, ha="left", va="bottom")

    ax.set_xlabel(r"$\Delta t$ (ps)")
    ax.set_ylabel(r"$T_{\mathrm{eff}}/T$")
    ax.set_title("Equilibrium effective temperature (measured, 60 ns)")
    ax.set_xlim(-0.015, 0.43)
    ax.set_ylim(0.42, 1.13)
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)


# ══════════════════════════════════════════════════════════════════════════
PANELS = [
    ("integrator_panel_geometry.png", panel_geometry),
    ("integrator_panel_order.png", panel_order),
    ("integrator_panel_drift.png", panel_drift),
    ("integrator_panel_teff.png", panel_teff),
]

for fname, fn in PANELS:
    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    fn(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print("wrote", os.path.join(OUT, fname))

# ── article-ready 2x2 composite (no banner; titles top-aligned via a fixed
# GridSpec + panel_geometry's anchor='N'; panel letters are added later in
# PPT per the figure convention, so none are baked here) ─────────────────
fig, axes = plt.subplots(2, 2, figsize=(PANEL_SIZE[0] * 2, PANEL_SIZE[1] * 2),
                         gridspec_kw=dict(hspace=0.34, wspace=0.40,
                                          left=0.065, right=0.90,
                                          top=0.955, bottom=0.07))
for (fname, fn), ax in zip(PANELS, axes.ravel()):
    fn(ax)
# panel letters (a)-(d) at each panel's top-left, in figure coordinates so
# they clear the panel titles (this figure is auto-composited, so the letters
# are baked here rather than added in a separate compositing pass)
for ax, letter in zip(axes.ravel(), "abcd"):
    ax.annotate(f"({letter})", xy=(0.0, 1.0), xytext=(-0.11, 1.13),
                xycoords="axes fraction", fontsize=15, fontweight="bold",
                color=CHARCOAL, va="top", ha="left")
composite = os.path.join(OUT, "integrator_panels_composite.png")
# bbox_inches='tight' expands the canvas to include the right-column direct
# labels (which extend past the axes); relative title alignment is preserved
fig.savefig(composite, facecolor="white", bbox_inches="tight", pad_inches=0.08)
plt.close(fig)
print("wrote", composite)

# mirror into the manuscript figure folder as the next chapter figure
import shutil
dest = os.path.join(ROOT, "article", "figs", "Chapter02_local_07.png")
shutil.copy(composite, dest)
print("mirrored", dest)
