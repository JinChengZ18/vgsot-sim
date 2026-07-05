"""Merge the two halves of the Experiment-B 60 ns full run into the final
thesis artifact.

The original --full run was killed by a system reboot after completing all
four cayley_explicit cells (values transcribed to cayley_from_log.json; the
run is fully seeded so they are exact). The remainder (euler_spherical,
cayley_true_midpoint) completed as tag=full_rest2. This script re-does the
weighted Richardson dt->0 fit for the cayley_explicit rows and emits the
merged three-integrator figure + verdict JSON.

Run:  PYTHONPATH=src python scripts/11_boltzmann_teff/merge_full_results.py
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.join("result", "sec_2_2_3_2", "B")


def wfit(x, y, s):
    W = 1.0 / np.asarray(s) ** 2
    X = np.vstack([np.ones(len(x)), x]).T
    WX = X * W[:, None]
    cov = np.linalg.inv(X.T @ WX)
    a, b = cov @ (WX.T @ np.asarray(y))
    sa, sb = np.sqrt(np.diag(cov))
    return dict(intercept=float(a), intercept_sigma=float(sa),
                slope=float(b), slope_sigma=float(sb),
                intercept_dev_sigma=float((a - 1.0) / sa),
                slope_pos_sigma=float(b / sb))


with open(os.path.join(OUT, "cayley_from_log.json"), encoding="utf-8") as f:
    cay = json.load(f)
with open(os.path.join(OUT, "teff_dtsweep_full_rest2.json"), encoding="utf-8") as f:
    rest = json.load(f)

sweep = {"cayley_explicit": cay["sweep"]}
sweep.update(rest["sweep"])

extrap = dict(rest["extrapolation"])
x = [r["dt_ps"] for r in cay["sweep"]]
extrap["cayley_explicit"] = {
    est: wfit(x, [r[est]["point"] for r in cay["sweep"]],
              [r[est]["sigma"] for r in cay["sweep"]])
    for est in ("moment", "slope")
}

pub = extrap["cayley_explicit"]["moment"]
verdict = dict(
    published_integrator="cayley_explicit",
    sc5_false=bool(abs(pub["intercept_dev_sigma"]) >= 3.0
                   or pub["slope_pos_sigma"] >= 3.0),
    intercept=f"{pub['intercept']:.4f}±{pub['intercept_sigma']:.4f} "
              f"({pub['intercept_dev_sigma']:+.1f}σ from 1)",
    slope=f"{pub['slope']:+.4f}±{pub['slope_sigma']:.4f}/ps "
          f"({pub['slope_pos_sigma']:+.1f}σ from 0)",
    euler_intercept=f"{extrap['euler_spherical']['moment']['intercept']:.4f}"
                    f"±{extrap['euler_spherical']['moment']['intercept_sigma']:.4f} "
                    f"({extrap['euler_spherical']['moment']['intercept_dev_sigma']:+.1f}σ)",
    midpoint_intercept=f"{extrap['cayley_true_midpoint']['moment']['intercept']:.4f}"
                       f"±{extrap['cayley_true_midpoint']['moment']['intercept_sigma']:.4f}",
)

merged = dict(meta=dict(source=["cayley_from_log.json (seeded, exact)",
                                "teff_dtsweep_full_rest2.json"],
                        T_phys_ns=60.0, n_traj=48, n_boot=2000),
              sweep=sweep, extrapolation=extrap, verdict=verdict)
with open(os.path.join(OUT, "teff_full_merged.json"), "w", encoding="utf-8") as f:
    json.dump(merged, f, indent=2)

colors = {"cayley_explicit": "#c0392b", "euler_spherical": "#2980b9",
          "cayley_true_midpoint": "#27ae60"}
labels = {"cayley_explicit": "Cayley explicit-$\\omega$ (published)",
          "euler_spherical": "Euler spherical (missing Wong–Zakai drift)",
          "cayley_true_midpoint": "Cayley true midpoint (control)"}

fig, ax = plt.subplots(figsize=(7.2, 5.0))
for integ, rows in sweep.items():
    xs = [r["dt_ps"] for r in rows]
    ys = [r["moment"]["point"] for r in rows]
    es = [r["moment"]["sigma"] for r in rows]
    ax.errorbar(xs, ys, yerr=es, marker="o", capsize=3,
                color=colors[integ], label=labels[integ])
    e = extrap[integ]["moment"]
    xx = np.linspace(0, 0.42, 40)
    ax.plot(xx, e["intercept"] + e["slope"] * xx, ":", color=colors[integ], alpha=0.7)
    ax.scatter([0], [e["intercept"]], marker="*", s=130, color=colors[integ],
               zorder=5, edgecolor="k", linewidth=0.5)
ax.axhline(1.0, color="k", lw=1, alpha=0.6)
ax.axhline(0.5, color="k", lw=0.8, ls="--", alpha=0.4)
ax.text(0.405, 0.503, "T/2", fontsize=9, va="bottom", ha="right", alpha=0.7)
ax.set_xlabel(r"$\Delta t$ (ps)")
ax.set_ylabel(r"$T_{\mathrm{eff}}/T$  (moment-inversion estimator)")
ax.set_title("Equilibrium effective temperature, 60 ns, "
             r"$\Delta=48.5$ ($\star$ = $\Delta t\to0$ extrapolation)")
ax.legend(fontsize=9, loc="center right")
ax.grid(alpha=0.3)
ax.set_xlim(left=-0.012)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "teff_full_merged.png"), dpi=150)

print(json.dumps(verdict, indent=2, ensure_ascii=False))
print("wrote teff_full_merged.{json,png}")
