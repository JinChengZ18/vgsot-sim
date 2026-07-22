"""E4 — symmetry-breaking bias-field calibre check.

The experiment measured at H_x = 200 Oe (section 2.3.1) while the simulator's
default exchange bias is -50 Oe along -y (section 2.2.1.2 footnote), a 4x gap
the chapter never discloses; maintenance docs record that 200 Oe produces a
~0.4 SER tail attributed to channels not in the model (Neel/DMI, roughness).
This experiment maps both knobs:

1. Field magnitude: P_sw(I) and the 50% threshold at |H_ex| in
   {0, 25, 50, 100, 150, 200, 300} Oe along the code convention -y
   (perpendicular to sigma_SH = -x). Quantifies how much of the 50<->200 Oe
   difference the effective theta_SH calibration absorbed, shows the H = 0
   random-bit limit, and reproduces (or not) the documented 200 Oe
   degradation with committed data.
2. Angle tolerance: rotate H_ex (50 Oe) in-plane away from perpendicular by
   0-45 deg and probe P_sw at 1160/1400 uA - quantifying the footnote's
   "deviation from perpendicular collapses SER to ~0.5" claim into a
   usable alignment window.

Modes:
  python run_e4.py --field 200        # one shard per field -> e4_h200.json
  python run_e4.py --angles           # -> e4_angles.json
  python run_e4.py --analyze          # -> e4_results.json + panels
  python run_e4.py --smoke
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
FIELDS_OE = (0.0, 25.0, 50.0, 100.0, 150.0, 200.0, 300.0)
FIELD_MULT = np.array([0.55, 0.65, 0.75, 0.83, 0.90, 0.95, 1.00, 1.05, 1.12, 1.20, 1.30])
ANGLES_DEG = (0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0)
ANGLE_CURRENTS_UA = (1160.0, 1400.0)
OE_TO_APM = 1000.0 / (4 * np.pi)


def _trial_seed(seed, i_sot, trial):
    return (int(abs(seed)) * 0x9E3779B1 + int(round(i_sot * 1e9)) * 0x85EBCA6B
            + int(trial) * 0xC2B2AE35) & 0x7FFFFFFF


def mc_psw(grid_uA, trials, seed, *, h_oe=50.0, angle_deg=0.0):
    """P_sw over an I_SOT grid with in-plane bias field of magnitude h_oe,
    rotated angle_deg away from the perpendicular (-y) convention toward
    the sigma_SH axis (-x)."""
    h = h_oe * OE_TO_APM
    ang = np.deg2rad(angle_deg)
    cc = replace(PhysicalConstantsConfig(),
                 h_ex_y=-h * np.cos(ang), h_ex_x=-h * np.sin(ang))
    mid1 = int(round(PULSE_NS * 1e-9 / cc.t_step))
    end = mid1 + int(round(RELAX_NS * 1e-9 / cc.t_step))
    psw = []
    t0 = time.time()
    for i_uA in grid_uA:
        i_sot = -i_uA * 1e-6
        ok = 0
        for tr in range(trials):
            rng = np.random.default_rng(_trial_seed(seed, i_sot, tr))
            res = run_piecewise_direct_excitation(
                sim_start_step=1, sim_mid1_step=mid1, sim_mid2_step=end,
                sim_end_step=end, pap=1,
                v_mtj_stage1=0.0, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
                i_sot_stage1=i_sot, i_sot_stage2=0.0, i_sot_stage3=0.0,
                estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1,
                estt_stage3=0, esot_stage3=1,
                vnv=0, non=1, r_sot_fl_dl=0.83, show_progress=False,
                constants=cc, enable_self_heating=True, T_ambient_K=300.0,
                integrator="cayley", rng=rng,
            )
            if abs(float(res.mz[end]) - 1.0) <= 0.2:
                ok += 1
        psw.append(ok / trials)
        print(f"    H={h_oe:.0f}Oe ang={angle_deg:g} |I|={i_uA:.0f}uA -> {psw[-1]:.3f}",
              flush=True)
    return dict(grid_uA=list(map(float, grid_uA)), psw=psw, trials=trials,
                h_oe=h_oe, angle_deg=angle_deg, R_W=PhysicalConstantsConfig().R_W,
                wall_s=time.time() - t0,
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
    for f in sorted(glob.glob(str(HERE / "e4_h*.json"))):
        runs.append(json.loads(Path(f).read_text(encoding="utf-8")))
    runs.sort(key=lambda d: d["h_oe"])
    if not runs:
        raise SystemExit("No e4_h*.json found.")
    cc = PhysicalConstantsConfig()

    rows = []
    for d in runs:
        ith, sig = crossing(d["grid_uA"], d["psw"], d["trials"])
        P = np.asarray(d["psw"])
        rows.append(dict(h_oe=d["h_oe"], ith_uA=ith, ith_sigma_uA=sig,
                         p_max=float(P.max()),
                         p_top_band=[float(P[-4:].min()), float(P[-4:].max())]))
    ith50 = next((r["ith_uA"] for r in rows if r["h_oe"] == 50.0), float("nan"))
    ith200 = next((r["ith_uA"] for r in rows if r["h_oe"] == 200.0), float("nan"))

    angles = None
    p_ang = HERE / "e4_angles.json"
    if p_ang.exists():
        angles = json.loads(p_ang.read_text(encoding="utf-8"))

    summary = dict(
        fields=rows,
        vth50_mV=float(ith50 * 1e-6 * cc.R_W * 1e3) if np.isfinite(ith50) else None,
        vth200_mV=float(ith200 * 1e-6 * cc.R_W * 1e3) if np.isfinite(ith200) else None,
        ratio_200_50=float(ith200 / ith50) if np.isfinite(ith200) and np.isfinite(ith50) else None,
        angles=angles,
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    (HERE / "e4_results.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    CHARCOAL, NAVY, CRIMSON, TEAL, AMBER = "#2B2B2B", "#1F5FA8", "#A82038", "#1A6B5A", "#C47A00"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans"],
        "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13.5,
        "mathtext.fontset": "stix", "axes.linewidth": 0.9,
        "axes.edgecolor": CHARCOAL, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "figure.dpi": 150, "savefig.dpi": 300,
    })

    # Panel 1: P_sw(I) family across field magnitudes.
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    cmap = plt.get_cmap("viridis")
    for i, d in enumerate(runs):
        c = cmap(i / max(1, len(runs) - 1) * 0.9)
        ax.plot(d["grid_uA"], d["psw"], "o-", ms=4.2, lw=1.25, color=c,
                label=rf"$|H_{{\mathrm{{ex}}}}|$ = {d['h_oe']:.0f} Oe")
    ax.axhline(0.5, color="gray", lw=0.6, ls=":")
    ax.set_xlabel(r"$|I_{\mathrm{SOT}}|$ ($\mu$A)")
    ax.set_ylabel(r"$P_{\mathrm{sw}}$")
    ax.set_title("Switching transition vs symmetry-breaking field magnitude")
    ax.set_ylim(-0.04, 1.06)
    ax.legend(fontsize=8.6, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "e4_field_family.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Panel 2: angle tolerance at 50 Oe.
    if angles:
        fig, ax = plt.subplots(figsize=(5.8, 4.5))
        colors = {1160.0: NAVY, 1400.0: CRIMSON}
        for i_uA in ANGLE_CURRENTS_UA:
            xs, ys = [], []
            for rec in angles["records"]:
                if rec["i_uA"] == i_uA:
                    xs.append(rec["angle_deg"]); ys.append(rec["psw"])
            ax.plot(xs, ys, "o-", ms=6, lw=1.5, color=colors.get(i_uA, CHARCOAL),
                    label=rf"$|I_{{\mathrm{{SOT}}}}|$ = {i_uA:.0f} $\mu$A")
        ax.axhline(0.5, color="gray", lw=0.6, ls=":")
        ax.set_xlabel(r"bias-field angle from $\perp\hat\sigma_{\mathrm{SH}}$ (deg)")
        ax.set_ylabel(r"$P_{\mathrm{sw}}$")
        ax.set_title("Alignment tolerance of deterministic switching (50 Oe)")
        ax.set_ylim(-0.04, 1.06)
        ax.legend(fontsize=9.5)
        fig.tight_layout()
        fig.savefig(HERE / "e4_angle_tolerance.png", bbox_inches="tight", facecolor="white")
        plt.close(fig)

    print("=" * 72)
    print("E4 summary (0.75 ns, AP->P, Cayley, self-heating ON, vnv=0)")
    print("=" * 72)
    print(f"{'H (Oe)':>8} {'Ith (uA)':>10} {'p_max':>7} {'top band':>15}")
    for r in rows:
        band = f"{r['p_top_band'][0]:.2f}-{r['p_top_band'][1]:.2f}"
        print(f"{r['h_oe']:>8.0f} {r['ith_uA']:>10.0f} {r['p_max']:>7.2f} {band:>15}")
    if summary["ratio_200_50"]:
        print(f"Ith(200 Oe)/Ith(50 Oe) = {summary['ratio_200_50']:.3f}")
    if angles:
        print("angle probe:", ", ".join(
            f"{rec['angle_deg']:g}deg@{rec['i_uA']:.0f}uA:{rec['psw']:.2f}"
            for rec in angles["records"]))
    print("Wrote e4_results.json + panels")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--field", type=float, default=None)
    ap.add_argument("--trials", type=int, default=150)
    ap.add_argument("--angles", action="store_true")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        out = mc_psw([900.0, 1160.0, 1400.0], 8, BASE_SEED + 4900, h_oe=50.0)
        print(json.dumps({k: out[k] for k in ("grid_uA", "psw", "wall_s")}, indent=2))
        out2 = mc_psw([1160.0], 8, BASE_SEED + 4901, h_oe=50.0, angle_deg=30.0)
        print("[smoke] 30deg probe:", out2["psw"])
        return
    if args.field is not None:
        H = args.field
        grid = np.round(ITH0_UA * FIELD_MULT)
        print(f"[E4] field {H:g} Oe: {len(grid)} pts x {args.trials} trials", flush=True)
        out = mc_psw(grid, args.trials, BASE_SEED + 4000 + int(H), h_oe=H)
        p = HERE / f"e4_h{H:g}.json"
        p.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[E4] wrote {p.name} (wall {out['wall_s']/60:.1f} min)", flush=True)
        return
    if args.angles:
        records = []
        for ang in ANGLES_DEG:
            for i_uA in ANGLE_CURRENTS_UA:
                out = mc_psw([i_uA], 200, BASE_SEED + 4500 + int(ang * 10), h_oe=50.0,
                             angle_deg=ang)
                records.append(dict(angle_deg=ang, i_uA=i_uA, psw=out["psw"][0],
                                    trials=200))
        p = HERE / "e4_angles.json"
        p.write_text(json.dumps(dict(h_oe=50.0, records=records), indent=2),
                     encoding="utf-8")
        print(f"[E4] wrote {p.name}", flush=True)
        return
    if args.analyze:
        analyze()
        return
    ap.error("choose one of --field / --angles / --analyze / --smoke")


if __name__ == "__main__":
    main()
