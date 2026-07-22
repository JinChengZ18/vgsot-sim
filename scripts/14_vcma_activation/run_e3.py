"""E3 — activate the VCMA branch: I_th(V_MTJ) map + dual-pulse energy check.

The platform is named vgsot-sim, section 2.1.3 builds the SOT-VCMA joint
drive model, and section 2.2.4 promises a dual-pulse scenario "used to
quantify the energy-optimization benefit" — yet no chapter figure ever
exercises the VCMA path. This experiment runs it, in two configurations:

1. Simultaneous assist: V_MTJ held during the 0.75 ns SOT write stage
   (vnv=1, v_mtj_stage1=V, 0 during relaxation). Output: I_th(V_MTJ) for
   V in {-0.8, -0.4, 0, +0.4, +0.8} V, compared against the analytic
   linear-VCMA scaling I_c0(V)/I_c0(0) = 1 - dKv(V)/K_U^eff with
   dKv = beta*V/(t_ox*t_f).
2. Sequential dual pulse (the 2.2.4 preset's literal description): stage1
   VCMA-only pre-pulse (1 ns, +0.8 V, I_SOT=0), stage2 SOT-only write
   (0.75 ns), stage3 relaxation. Without a persistent barrier effect the
   pre-pulse should buy nothing — an honest test of the preset's premise.

Write energy at the p = 0.5 workpoint is compared in --analyze from the
extracted thresholds: E = I_th^2 R_W t_w + V_MTJ^2/R_P t_V.

Modes (one detached process per V_MTJ value parallelizes):
  python run_e3.py --vmtj 0.8         # -> e3_v0.8.json
  python run_e3.py --dualpulse        # -> e3_dualpulse.json
  python run_e3.py --analyze          # -> e3_results.json + panels
  python run_e3.py --smoke
"""
from __future__ import annotations

import argparse
import glob
import json
import time
from datetime import datetime
from pathlib import Path

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation

HERE = Path(__file__).resolve().parent
BASE_SEED = 20260719

PULSE_NS = 0.75
RELAX_NS = 3.25
PRECHARGE_NS = 1.0            # VCMA-only pre-pulse length (sequential mode)
ITH0_UA = 1160.0              # simulated V_MTJ=0 threshold (fig 2.11 / E1 anchor)
VMTJ_LIST = (-0.8, -0.4, 0.0, 0.4, 0.8)
# Wide sign-dependent windows anchored on I_th(0) rather than on the analytic
# static scaling: whether the dynamic 0.75 ns threshold follows the static
# K_U^eff scaling is exactly what the experiment measures, so the windows
# must not presume it.
POS_MULT = np.array([0.48, 0.56, 0.64, 0.72, 0.80, 0.87, 0.93, 0.98, 1.02, 1.06])
NEG_MULT = np.array([0.97, 1.03, 1.09, 1.16, 1.24, 1.32, 1.42, 1.52, 1.64])
ZERO_MULT = np.array([0.85, 0.90, 0.94, 0.97, 1.00, 1.03, 1.06, 1.10, 1.16])


def window_for(v_mtj: float) -> np.ndarray:
    m = POS_MULT if v_mtj > 0 else NEG_MULT if v_mtj < 0 else ZERO_MULT
    return np.round(ITH0_UA * m)


def vcma_scale(cc: PhysicalConstantsConfig, v_mtj: float) -> float:
    """Analytic linear-VCMA threshold scaling I_c0(V)/I_c0(0) = 1 - dKv/K_U^eff.

    K_U^eff is reconstructed from the same ingredients the kernel uses:
    interface term Ki/tf minus the ellipsoid shape term."""
    from vgsot_sim.demag import demag_factors

    nx, ny, nz = demag_factors(cc, mode="ellipsoid")
    k_shape = 0.5 * cc.u0 * cc.Ms**2 * (nz - 0.5 * (nx + ny))
    k_ueff = cc.Ki / cc.tf - k_shape
    d_kv = cc.beta * v_mtj / (cc.tox * cc.tf)
    return 1.0 - d_kv / k_ueff


def _trial_seed(seed, i_sot, trial):
    return (int(abs(seed)) * 0x9E3779B1 + int(round(i_sot * 1e9)) * 0x85EBCA6B
            + int(trial) * 0xC2B2AE35) & 0x7FFFFFFF


def mc_psw(grid_uA, trials, seed, *, v_mtj_write=0.0, precharge_ns=0.0,
           v_mtj_precharge=0.0):
    """P_sw over an I_SOT grid with stage-wise V_MTJ.

    Stage layout: [precharge (VCMA only)] -> [write: I_SOT (+ V_MTJ if
    simultaneous)] -> [relax: all off]. vnv=1 throughout so any nonzero
    stage V_MTJ acts on the barrier; estt=0 keeps STT off (SOT-only drive,
    same as every chapter workpoint).
    """
    cc = PhysicalConstantsConfig()
    pre = int(round(precharge_ns * 1e-9 / cc.t_step))
    # Stage layout is FIXED: stage1 = precharge (VCMA only; a single idle
    # step when there is no precharge), stage2 = write, stage3 = relax.
    mid1 = max(pre, 1)
    mid2 = mid1 + int(round(PULSE_NS * 1e-9 / cc.t_step))
    end = mid2 + int(round(RELAX_NS * 1e-9 / cc.t_step))
    psw = []
    t0 = time.time()
    for i_uA in grid_uA:
        i_sot = -i_uA * 1e-6
        ok = 0
        for tr in range(trials):
            rng = np.random.default_rng(_trial_seed(seed, i_sot, tr))
            res = run_piecewise_direct_excitation(
                sim_start_step=1, sim_mid1_step=mid1, sim_mid2_step=mid2,
                sim_end_step=end, pap=1,
                v_mtj_stage1=(v_mtj_precharge if pre > 0 else 0.0),
                v_mtj_stage2=v_mtj_write,
                v_mtj_stage3=0.0,
                i_sot_stage1=0.0,
                i_sot_stage2=i_sot,
                i_sot_stage3=0.0,
                estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1,
                estt_stage3=0, esot_stage3=1,
                vnv=1, non=1, r_sot_fl_dl=0.83, show_progress=False,
                constants=cc, enable_self_heating=True, T_ambient_K=300.0,
                integrator="cayley", rng=rng,
            )
            if abs(float(res.mz[end]) - 1.0) <= 0.2:
                ok += 1
        psw.append(ok / trials)
        print(f"    |I|={i_uA:.0f} uA -> P_sw={psw[-1]:.3f}", flush=True)
    return dict(grid_uA=list(map(float, grid_uA)), psw=psw, trials=trials,
                R_W=cc.R_W, wall_s=time.time() - t0,
                timestamp=datetime.now().isoformat(timespec="seconds"))


def crossing(grid_uA, P, trials):
    """Raw 50% crossing in uA + binomial-propagated sigma (same calibre as E1)."""
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

    cc = PhysicalConstantsConfig()
    raw = []
    for f in sorted(glob.glob(str(HERE / "e3_v*.json"))):
        raw.append(json.loads(Path(f).read_text(encoding="utf-8")))
    if not raw:
        raise SystemExit("No e3_v*.json found.")
    # Merge files sharing a V_MTJ (base window + low-current supplements).
    from collections import defaultdict
    groups = defaultdict(list)
    for d in raw:
        groups[d["v_mtj"]].append(d)
    runs = []
    for v in sorted(groups):
        ds = groups[v]
        I = np.concatenate([np.asarray(d["grid_uA"], float) for d in ds])
        P = np.concatenate([np.asarray(d["psw"], float) for d in ds])
        o = np.argsort(I)
        runs.append(dict(v_mtj=v, grid_uA=I[o].tolist(), psw=P[o].tolist(),
                         trials=ds[0]["trials"]))
    V = np.array([d["v_mtj"] for d in runs])
    cross = [crossing(d["grid_uA"], d["psw"], d["trials"]) for d in runs]
    ith = np.array([c[0] for c in cross])
    ith_sig = np.array([c[1] for c in cross])
    scale_sim = ith / ith[np.argmin(np.abs(V))]
    scale_ana = np.array([vcma_scale(cc, v) for v in V])

    dual = None
    p_dual = HERE / "e3_dualpulse.json"
    if p_dual.exists():
        dual = json.loads(p_dual.read_text(encoding="utf-8"))
        dual_ith, dual_sig = crossing(dual["grid_uA"], dual["psw"], dual["trials"])
    else:
        dual_ith = dual_sig = float("nan")

    # Write energy at the p=0.5 workpoint. SOT channel + MTJ leakage during
    # the stages where V_MTJ is applied (R_P = RA / A_elec).
    A_elec = np.pi / 4 * cc.D_elec**2
    R_P = cc.RA / A_elec        # RA is in Ohm*m^2 -> ~5 kOhm at D_elec = 65 nm
    tw = PULSE_NS * 1e-9
    def energy(i_uA, v_mtj, t_v):
        return (i_uA * 1e-6) ** 2 * cc.R_W * tw + v_mtj**2 / R_P * t_v
    e_rows = []
    for v, i in zip(V, ith):
        if np.isfinite(i):
            e_rows.append(dict(v_mtj=float(v), ith_uA=float(i),
                               E_pJ=float(energy(i, v, tw) * 1e12)))
    i0 = float(ith[np.argmin(np.abs(V))])
    e0 = energy(i0, 0.0, 0.0) * 1e12

    summary = dict(
        v_mtj=V.tolist(), ith_uA=ith.tolist(), ith_sigma_uA=ith_sig.tolist(),
        scale_sim=scale_sim.tolist(), scale_analytic=scale_ana.tolist(),
        dualpulse=dict(ith_uA=dual_ith, sigma_uA=dual_sig,
                       precharge_ns=PRECHARGE_NS,
                       raw=dual) if dual else None,
        energy_rows=e_rows, E_sot_only_pJ=float(e0), R_P_ohm=float(R_P),
        R_W_ohm=float(cc.R_W), beta_vcma=cc.beta,
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    (HERE / "e3_results.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    CHARCOAL, NAVY, CRIMSON, TEAL, AMBER = "#2B2B2B", "#1F5FA8", "#A82038", "#1A6B5A", "#C47A00"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans"],
        "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13.5,
        "mathtext.fontset": "stix", "axes.linewidth": 0.9,
        "axes.edgecolor": CHARCOAL, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "figure.dpi": 150, "savefig.dpi": 300,
    })

    # Panel 1: I_th vs V_MTJ, sim vs analytic linear-VCMA scaling.
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    vv = np.linspace(-0.9, 0.9, 200)
    ax.plot(vv, [vcma_scale(cc, v) * i0 for v in vv], "--", color=TEAL, lw=1.7,
            label=r"analytic  $I_{c0}(V)\propto K_U^{\mathrm{eff}}-\beta_{\mathrm{VCMA}}V/(t_{\mathrm{ox}}t_f)$")
    ax.errorbar(V, ith, yerr=np.maximum(ith_sig, 1.0), fmt="o", ms=7.5, color=CRIMSON,
                mfc="white", mew=1.6, capsize=3, lw=1.4,
                label="vgsot-sim, simultaneous assist")
    if dual and np.isfinite(dual_ith):
        ax.errorbar([0.8], [dual_ith], yerr=[max(dual_sig, 1.0)], fmt="D", ms=8,
                    color=AMBER, mfc="white", mew=1.6, capsize=3,
                    label=f"sequential dual pulse (+0.8 V, {PRECHARGE_NS:g} ns pre)")
    ax.set_xlabel(r"$V_{\mathrm{MTJ}}$ (V)")
    ax.set_ylabel(r"$|I_{\mathrm{th}}|$ ($\mu$A)")
    ax.set_title("VCMA modulation of the 0.75 ns switching threshold")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "e3_ith_vs_vmtj.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Panel 2: the P_sw(I) family.
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    cmap = plt.get_cmap("coolwarm")
    for d in runs:
        c = cmap((d["v_mtj"] + 0.8) / 1.6)
        ax.plot(d["grid_uA"], d["psw"], "o-", ms=4.5, lw=1.3, color=c,
                label=rf"$V_{{\mathrm{{MTJ}}}}$ = {d['v_mtj']:+.1f} V")
    ax.axhline(0.5, color="gray", lw=0.6, ls=":")
    ax.set_xlabel(r"$|I_{\mathrm{SOT}}|$ ($\mu$A)")
    ax.set_ylabel(r"$P_{\mathrm{sw}}$")
    ax.set_title(r"$P_{\mathrm{sw}}(I_{\mathrm{SOT}})$ under simultaneous VCMA bias")
    ax.set_ylim(-0.04, 1.06)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "e3_psw_family.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print("=" * 70)
    print("E3 summary (0.75 ns, AP->P, Cayley, self-heating ON, vnv=1)")
    print("=" * 70)
    print(f"{'V_MTJ (V)':>10} {'Ith (uA)':>10} {'sim scale':>10} {'analytic':>9} {'E (pJ)':>8}")
    for v, i, ss, sa in zip(V, ith, scale_sim, scale_ana):
        e = energy(i, v, tw) * 1e12 if np.isfinite(i) else float("nan")
        print(f"{v:>10.1f} {i:>10.0f} {ss:>10.3f} {sa:>9.3f} {e:>8.3f}")
    if dual:
        print(f"sequential dual pulse: Ith = {dual_ith:.0f} uA "
              f"(V=0 reference {i0:.0f} uA) -> benefit "
              f"{(1 - dual_ith/i0)*100:+.1f}%")
    print(f"R_P = {R_P/1e3:.2f} kOhm, R_W = {cc.R_W:.0f} Ohm, "
          f"E_SOT-only(p=0.5) = {e0:.3f} pJ")
    print("Wrote e3_results.json, e3_ith_vs_vmtj.png, e3_psw_family.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vmtj", type=float, default=None)
    ap.add_argument("--trials", type=int, default=150)
    ap.add_argument("--dualpulse", action="store_true")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    cc = PhysicalConstantsConfig()
    if args.smoke:
        print("[smoke] V_MTJ=+0.8, 3 points, 8 trials")
        grid = np.array([700.0, 900.0, 1100.0])
        out = mc_psw(grid, 8, BASE_SEED + 3900, v_mtj_write=0.8)
        print(json.dumps({k: out[k] for k in ("grid_uA", "psw", "wall_s")}, indent=2))
        print(f"[smoke] analytic scale(+0.8 V) = {vcma_scale(cc, 0.8):.3f}")
        return
    if args.vmtj is not None:
        v = args.vmtj
        grid = window_for(v)
        print(f"[E3] V_MTJ={v:+.1f} V: {len(grid)} pts x {args.trials} trials "
              f"(analytic scale {vcma_scale(cc, v):.3f})", flush=True)
        out = mc_psw(grid, args.trials, BASE_SEED + 3000 + int(round(v * 10)),
                     v_mtj_write=v)
        out["v_mtj"] = v
        p = HERE / f"e3_v{v:g}.json"
        p.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[E3] wrote {p.name} (wall {out['wall_s']/60:.1f} min)", flush=True)
        return
    if args.dualpulse:
        grid = window_for(0.0)
        print(f"[E3] sequential dual pulse: +0.8 V x {PRECHARGE_NS} ns pre, "
              f"{len(grid)} pts x {args.trials} trials", flush=True)
        out = mc_psw(grid, args.trials, BASE_SEED + 3500,
                     v_mtj_write=0.0, precharge_ns=PRECHARGE_NS, v_mtj_precharge=0.8)
        p = HERE / "e3_dualpulse.json"
        p.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[E3] wrote {p.name} (wall {out['wall_s']/60:.1f} min)", flush=True)
        return
    if args.analyze:
        analyze()
        return
    ap.error("choose one of --vmtj / --dualpulse / --analyze / --smoke")


if __name__ == "__main__":
    main()
