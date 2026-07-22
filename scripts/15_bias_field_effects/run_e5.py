"""E5 — direction asymmetry vs reference-layer stray field.

Section 2.3.4 attributes the measured direction asymmetry (P->AP shows a
back-hopping plateau from ~940 mV and two-segment transitions; AP->P is
clean; hysteresis intercepts differ by ~30 mV) to an uncompensated
reference-layer stray field, with "可能提示" wording and no simulation
support. This experiment injects a z stray field H_z in {0, 25, 50, 100} Oe
and scans BOTH switching directions, testing which of the three measured
features the macrospin picture reproduces:
  (i)  threshold asymmetry between directions growing with H_z;
  (ii) direction-selective over-drive plateau (back-hopping) suppression /
       enhancement;
  (iii) two-segment transitions (not expected from a single macrospin).

Directions at the fig 2.11 workpoint (0.75 ns, Cayley, self-heating ON):
  p2ap: pap=1 (start m_z ~ -1 = P), drive -I, target m_z = +1.
  ap2p: pap=0 (start m_z ~ +1 = AP), drive +I, target m_z = -1.

Modes:
  python run_e5.py --direction p2ap --hstray 50    # -> e5_p2ap_h50.json
  python run_e5.py --analyze                       # -> e5_results.json + panels
  python run_e5.py --smoke
"""
from __future__ import annotations

import argparse
import glob
import json
import time
from dataclasses import replace
from datetime import datetime
from pathlib import Path

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation

HERE = Path(__file__).resolve().parent
BASE_SEED = 20260719

PULSE_NS, RELAX_NS = 0.75, 3.25
ITH0_UA = 1160.0
HSTRAY_OE = (0.0, 25.0, 50.0, 100.0)
MULT = np.array([0.85, 0.90, 0.94, 0.97, 1.00, 1.03, 1.06, 1.10, 1.16, 1.25, 1.35])
OE_TO_APM = 1000.0 / (4 * np.pi)


def _trial_seed(seed, i_sot, trial):
    return (int(abs(seed)) * 0x9E3779B1 + int(round(i_sot * 1e9)) * 0x85EBCA6B
            + int(trial) * 0xC2B2AE35) & 0x7FFFFFFF


def mc_psw(grid_uA, trials, seed, *, direction="p2ap", hstray_oe=0.0):
    cc = replace(PhysicalConstantsConfig(), h_ex_z=hstray_oe * OE_TO_APM)
    mid1 = int(round(PULSE_NS * 1e-9 / cc.t_step))
    end = mid1 + int(round(RELAX_NS * 1e-9 / cc.t_step))
    pap, sgn, target = (1, -1.0, +1.0) if direction == "p2ap" else (0, +1.0, -1.0)
    psw = []
    t0 = time.time()
    for i_uA in grid_uA:
        i_sot = sgn * i_uA * 1e-6
        ok = 0
        for tr in range(trials):
            rng = np.random.default_rng(_trial_seed(seed, i_sot, tr))
            res = run_piecewise_direct_excitation(
                sim_start_step=1, sim_mid1_step=mid1, sim_mid2_step=end,
                sim_end_step=end, pap=pap,
                v_mtj_stage1=0.0, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
                i_sot_stage1=i_sot, i_sot_stage2=0.0, i_sot_stage3=0.0,
                estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1,
                estt_stage3=0, esot_stage3=1,
                vnv=0, non=1, r_sot_fl_dl=0.83, show_progress=False,
                constants=cc, enable_self_heating=True, T_ambient_K=300.0,
                integrator="cayley", rng=rng,
            )
            if abs(float(res.mz[end]) - target) <= 0.2:
                ok += 1
        psw.append(ok / trials)
        print(f"    {direction} Hz={hstray_oe:.0f}Oe |I|={i_uA:.0f}uA -> {psw[-1]:.3f}",
              flush=True)
    return dict(grid_uA=list(map(float, grid_uA)), psw=psw, trials=trials,
                direction=direction, hstray_oe=hstray_oe,
                R_W=PhysicalConstantsConfig().R_W, wall_s=time.time() - t0,
                timestamp=datetime.now().isoformat(timespec="seconds"))


def crossing(grid_uA, P, trials):
    I = np.asarray(grid_uA, float); P = np.asarray(P, float)
    o = np.argsort(I); I, P = I[o], P[o]
    for i in range(len(I) - 1):
        if (P[i] - 0.5) * (P[i + 1] - 0.5) <= 0 and P[i] != P[i + 1]:
            slope = (P[i + 1] - P[i]) / (I[i + 1] - I[i])
            return float(I[i] + (0.5 - P[i]) / slope), float(np.sqrt(0.25 / trials) / abs(slope))
    return float("nan"), float("nan")


def analyze():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = []
    for f in sorted(glob.glob(str(HERE / "e5_*.json"))):
        if "results" in f:
            continue
        runs.append(json.loads(Path(f).read_text(encoding="utf-8")))
    if not runs:
        raise SystemExit("No e5 shard JSONs found.")
    cc = PhysicalConstantsConfig()

    rows = []
    for d in runs:
        ith, sig = crossing(d["grid_uA"], d["psw"], d["trials"])
        P = np.asarray(d["psw"])
        rows.append(dict(direction=d["direction"], hstray_oe=d["hstray_oe"],
                         ith_uA=ith, ith_sigma_uA=sig,
                         plateau_band=[float(P[-3:].min()), float(P[-3:].max())]))
    rows.sort(key=lambda r: (r["hstray_oe"], r["direction"]))
    asym = []
    for h in sorted({r["hstray_oe"] for r in rows}):
        a = next((r for r in rows if r["hstray_oe"] == h and r["direction"] == "ap2p"), None)
        p = next((r for r in rows if r["hstray_oe"] == h and r["direction"] == "p2ap"), None)
        if a and p and np.isfinite(a["ith_uA"]) and np.isfinite(p["ith_uA"]):
            asym.append(dict(hstray_oe=h,
                             dI_uA=float(a["ith_uA"] - p["ith_uA"]),
                             dV_mV=float((a["ith_uA"] - p["ith_uA"]) * 1e-6 * cc.R_W * 1e3)))

    summary = dict(rows=rows, asymmetry=asym,
                   measured_intercept_diff_mV=30.0,
                   timestamp=datetime.now().isoformat(timespec="seconds"))
    (HERE / "e5_results.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    CHARCOAL, NAVY, CRIMSON = "#2B2B2B", "#1F5FA8", "#A82038"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans"],
        "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13.5,
        "mathtext.fontset": "stix", "axes.linewidth": 0.9,
        "axes.edgecolor": CHARCOAL, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "figure.dpi": 150, "savefig.dpi": 300,
    })

    # Panel: P_sw(I) for both directions across stray fields.
    hs = sorted({r["hstray_oe"] for r in rows})
    fig, axes = plt.subplots(1, len(hs), figsize=(3.4 * len(hs), 4.2), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, h in zip(axes, hs):
        for d in runs:
            if d["hstray_oe"] != h:
                continue
            c, lbl = (NAVY, "P" + r"$\rightarrow$" + "AP") if d["direction"] == "p2ap" \
                else (CRIMSON, "AP" + r"$\rightarrow$" + "P")
            ax.plot(d["grid_uA"], d["psw"], "o-", ms=4, lw=1.2, color=c, label=lbl)
        ax.axhline(0.5, color="gray", lw=0.5, ls=":")
        ax.set_title(rf"$H_z$ = {h:.0f} Oe", fontsize=12)
        ax.set_xlabel(r"$|I_{\mathrm{SOT}}|$ ($\mu$A)")
        ax.set_ylim(-0.04, 1.06)
    axes[0].set_ylabel(r"$P_{\mathrm{sw}}$")
    axes[0].legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "e5_direction_asymmetry.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print("=" * 74)
    print("E5 summary (0.75 ns, Cayley, self-heating ON)")
    print("=" * 74)
    print(f"{'H_z (Oe)':>9} {'dir':>6} {'Ith (uA)':>10} {'plateau':>13}")
    for r in rows:
        band = f"{r['plateau_band'][0]:.2f}-{r['plateau_band'][1]:.2f}"
        print(f"{r['hstray_oe']:>9.0f} {r['direction']:>6} {r['ith_uA']:>10.0f} {band:>13}")
    for a in asym:
        print(f"asymmetry at H_z={a['hstray_oe']:.0f} Oe: dI={a['dI_uA']:+.0f} uA "
              f"(dV={a['dV_mV']:+.1f} mV; measured intercept diff ~30 mV)")
    print("Wrote e5_results.json, e5_direction_asymmetry.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--direction", choices=("p2ap", "ap2p"), default=None)
    ap.add_argument("--hstray", type=float, default=0.0)
    ap.add_argument("--trials", type=int, default=150)
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        for d in ("p2ap", "ap2p"):
            out = mc_psw([1160.0, 1400.0], 8, BASE_SEED + 5900, direction=d, hstray_oe=0.0)
            print(f"[smoke] {d} H_z=0:", out["psw"])
        return
    if args.analyze:
        analyze()
        return
    if args.direction is None:
        ap.error("choose --direction/--hstray, --analyze, or --smoke")
    grid = np.round(ITH0_UA * MULT)
    print(f"[E5] {args.direction} H_z={args.hstray:g} Oe: {len(grid)} pts x "
          f"{args.trials} trials", flush=True)
    out = mc_psw(grid, args.trials,
                 BASE_SEED + 5000 + int(args.hstray) + (500 if args.direction == "ap2p" else 0),
                 direction=args.direction, hstray_oe=args.hstray)
    p = HERE / f"e5_{args.direction}_h{args.hstray:g}.json"
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[E5] wrote {p.name} (wall {out['wall_s']/60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()
