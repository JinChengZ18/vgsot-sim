"""E2 — overlay the measured Sigmoid on the simulated P_sw(V) and decompose eta_c.

Fig 2.11 aligns simulation and experiment at the 50% threshold only; the
full measured Device-A AP->P Sigmoid (fig 2.14b) was never overlaid on the
simulated transition. This experiment (i) densely samples the simulated
transition (>=500 trials/point) and extracts beta_s^sLLG with a confidence
interval, placing the macrospin engine on the slope ladder
beta^NB (7.9 V^-1) -> beta_s^sLLG (?) -> beta^meas (44.6 V^-1), i.e. how
much of the C2C narrowing factor eta_c = 5.6 the macrospin dynamics
explains; and (ii) confronts the simulated over-drive back-hopping plateau
(~0.8, fig 2.11) with the measured clean saturation to 1.0 by ~1020 mV,
including a relaxation-window sensitivity check that separates judgement
convention from device physics.

Measured reference (mirrors scripts/07_process_variability/sigmoid_fig.py
DATA[("A", "AP->P")], 50 write-read cycles per point, H_x = 200 Oe):
sigmoid fit V_th = 894.0 mV, k = 22.43 mV, beta_s = 44.6 V^-1 (Table 2.6).
NB reference at tau0 = 1 ns: Delta = 4.91, V_c0 = 857 mV (Table 2.8).

Modes (parallel-friendly):
  python run_e2.py --dense-shard 0..3     # 13-point transition grid, 125 trials each
  python run_e2.py --overdrive            # 6 points 1250-1600 uA, 200 trials
  python run_e2.py --relaxwin 3.25|10|50  # window sensitivity at 1400/2000 uA
  python run_e2.py --analyze              # pool -> fits + figures + e2_results.json
  python run_e2.py --smoke
"""
from __future__ import annotations

import argparse
import glob
import json
import time
from datetime import datetime
from pathlib import Path

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig, SerSotNoVcmaThermalConfig
from vgsot_sim.ser_cases import ser_sot_no_vcma_thermal
from vgsot_sim.analysis.sigmoid_fit import sigmoid4p, wilson
from vgsot_sim.analysis import nb_fit

HERE = Path(__file__).resolve().parent
BASE_SEED = 20260719

# ── Measured Device A, AP->P, 0.75 ns (V_SOT mV, P_sw; N = 50 each) ──────
MEAS_V_MV = np.array([800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020])
MEAS_P = np.array([0.000, 0.020, 0.100, 0.180, 0.340, 0.500, 0.820, 0.840,
                   0.900, 0.980, 0.940, 1.000])
MEAS_N = 50          # per-shot records hold 50 cycles per point, not 100
MEAS_FIT = dict(Vth_mV=894.0, k_mV=22.43, beta=44.6)     # Table 2.6
NB_REF = dict(Delta=4.91, Vc0=0.857, tau0_ns=1.0)        # Table 2.8
BETA_NB_ANALYTIC = nb_fit.beta_nb_analytic(NB_REF["Delta"], NB_REF["Vc0"])  # ~7.94

PULSE_NS = 0.75
DENSE_GRID_UA = np.round(np.linspace(1085, 1235, 13))
OVERDRIVE_UA = np.array([1250.0, 1300.0, 1360.0, 1420.0, 1500.0, 1600.0])
RELAX_CURRENTS_UA = np.array([1400.0, 2000.0])


def run_scan(grid_uA, trials, seed, relax_ns=3.25):
    cc = PhysicalConstantsConfig()
    mid1 = int(round(PULSE_NS * 1e-9 / cc.t_step))
    end = mid1 + int(round(relax_ns * 1e-9 / cc.t_step))
    cfg = SerSotNoVcmaThermalConfig(
        i_sot_list=tuple(-np.asarray(grid_uA, float) * 1e-6), trials=trials,
        sim_start_step=1, sim_mid1_step=mid1, sim_end_step=end,
        pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
        target_mz=1.0, failure_tol=0.2, constants=cc,
    )
    t0 = time.time()
    res = ser_sot_no_vcma_thermal(
        cfg, show_progress=False, enable_self_heating=True, T_ambient_K=300.0,
        seed=seed, rng_mode="generator", integrator="cayley",
    )
    return dict(
        grid_uA=np.asarray(grid_uA, float).tolist(),
        V_mV=(np.asarray(grid_uA, float) * 1e-6 * cc.R_W * 1e3).tolist(),
        psw=res.psw.tolist(), trials=trials, relax_ns=relax_ns, seed=seed,
        R_W=cc.R_W, theta_SH=cc.theta_SH, integrator="cayley",
        self_heating=True, wall_s=time.time() - t0,
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )


def analyze():
    from scipy.optimize import curve_fit
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ── Pool dense shards ────────────────────────────────────────────────
    shards = []
    for f in sorted(glob.glob(str(HERE / "e2_dense_shard*.json"))):
        shards.append(json.loads(Path(f).read_text(encoding="utf-8")))
    if not shards:
        raise SystemExit("No e2_dense_shard*.json found.")
    grid = np.array(shards[0]["grid_uA"])
    V_mV = np.array(shards[0]["V_mV"])
    n_tot = sum(s["trials"] for s in shards)
    succ = np.zeros_like(grid, dtype=float)
    for s in shards:
        assert s["grid_uA"] == shards[0]["grid_uA"], "shard grids differ"
        succ += np.array(s["psw"]) * s["trials"]
    psw = succ / n_tot

    # 4-parameter logistic in volts -> beta = 1/k (V^-1)
    Vv = V_mV / 1e3
    p0 = (0.0, 1.0, float(Vv[np.argmin(np.abs(psw - 0.5))]), 0.02)
    bounds = ([-0.05, 0.3, Vv.min() - 0.05, 1e-3], [0.30, 1.10, Vv.max() + 0.05, 0.20])
    popt, pcov = curve_fit(sigmoid4p, Vv, psw, p0=p0, bounds=bounds, maxfev=30000)
    perr = np.sqrt(np.diag(pcov))
    y0, L, vth_sim, k_sim = (float(x) for x in popt)
    beta_sim = 1.0 / k_sim
    beta_sim_sigma = float(perr[3] / k_sim**2)
    ss_res = float(np.sum((psw - sigmoid4p(Vv, *popt))**2))
    ss_tot = float(np.sum((psw - psw.mean())**2))
    r2 = 1.0 - ss_res / ss_tot

    # ── Overdrive + relaxation-window data ───────────────────────────────
    over = None
    p_over = HERE / "e2_overdrive.json"
    if p_over.exists():
        over = json.loads(p_over.read_text(encoding="utf-8"))
    relax = []
    for f in sorted(glob.glob(str(HERE / "e2_relaxwin_*.json"))):
        relax.append(json.loads(Path(f).read_text(encoding="utf-8")))

    eta_meas = MEAS_FIT["beta"] / BETA_NB_ANALYTIC
    eta_sim = beta_sim / BETA_NB_ANALYTIC
    share = (np.log(eta_sim) / np.log(eta_meas)) if eta_meas > 1 else float("nan")

    summary = dict(
        dense=dict(grid_uA=grid.tolist(), V_mV=V_mV.tolist(), psw=psw.tolist(),
                   trials_per_point=n_tot, n_shards=len(shards)),
        sigmoid_sim=dict(y0=y0, L=L, Vth_mV=vth_sim * 1e3, k_mV=k_sim * 1e3,
                         beta=beta_sim, beta_sigma=beta_sim_sigma, R2=r2),
        slope_ladder=dict(beta_NB=BETA_NB_ANALYTIC, beta_sLLG=beta_sim,
                          beta_meas=MEAS_FIT["beta"],
                          eta_c_meas=eta_meas, eta_c_sim=eta_sim,
                          log_share_explained=share),
        overdrive=over, relaxwin=relax,
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    (HERE / "e2_results.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    def yerr_pair(P, lo, hi):
        """Non-negative (lower, upper) error-bar pair; clips the float-eps
        overshoot of the Wilson bounds at the degenerate p = 0/1 edges."""
        P = np.asarray(P, float)
        return [np.maximum(P - lo, 0.0), np.maximum(hi - P, 0.0)]

    # ── Style (repo academic convention; clean panels, no baked labels) ──
    CHARCOAL, NAVY, CRIMSON, TEAL, AMBER = "#2B2B2B", "#1F5FA8", "#A82038", "#1A6B5A", "#C47A00"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans"],
        "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13.5,
        "mathtext.fontset": "stix", "axes.linewidth": 0.9,
        "axes.edgecolor": CHARCOAL, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "figure.dpi": 150, "savefig.dpi": 300,
    })

    # Panel 1: measured vs simulated transition + NB reference.
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    lo_m, hi_m = wilson(MEAS_P, MEAS_N)
    ax.errorbar(MEAS_V_MV, MEAS_P, yerr=yerr_pair(MEAS_P, lo_m, hi_m), fmt="s",
                ms=6.5, color=TEAL, mfc="white", mew=1.5, capsize=2.5, lw=1.1,
                label=f"measured (N = {MEAS_N}/point)")
    vv = np.linspace(0.78, 1.26, 400)
    ax.plot(vv * 1e3, sigmoid4p(vv * 1e3, 0.0, 1.0, MEAS_FIT["Vth_mV"], MEAS_FIT["k_mV"]),
            "--", color=TEAL, lw=1.5,
            label=rf"measured fit  $\beta_s$ = {MEAS_FIT['beta']:.1f} V$^{{-1}}$")
    lo_s, hi_s = wilson(psw, n_tot)
    ax.errorbar(V_mV, psw, yerr=yerr_pair(psw, lo_s, hi_s), fmt="o", ms=6,
                color=CRIMSON, mfc="white", mew=1.5, capsize=2.5, lw=1.1,
                label=f"vgsot-sim (N = {n_tot}/point)")
    ax.plot(vv * 1e3, sigmoid4p(vv, *popt), "-", color=CRIMSON, lw=1.6,
            label=rf"sim fit  $\beta_s$ = {beta_sim:.1f} $\pm$ {1.96*beta_sim_sigma:.1f} V$^{{-1}}$")
    if over is not None:
        oV = np.array(over["V_mV"]); oP = np.array(over["psw"])
        lo_o, hi_o = wilson(oP, over["trials"])
        ax.errorbar(oV, oP, yerr=yerr_pair(oP, lo_o, hi_o), fmt="o", ms=6, color=CRIMSON,
                    mfc=CRIMSON, mew=0, capsize=2.5, lw=1.1, alpha=0.75,
                    label="vgsot-sim over-drive")
    ax.plot(vv * 1e3, nb_fit.psw_nb(vv, PULSE_NS, NB_REF["Delta"], NB_REF["Vc0"]),
            ":", color="gray", lw=1.6,
            label=rf"Neel-Brown, no C2C corr. ($\beta$ = {BETA_NB_ANALYTIC:.1f} V$^{{-1}}$)")
    ax.axhline(0.5, color="gray", lw=0.5, ls=":")
    ax.set_xlabel(r"$|V_{\mathrm{SOT}}|$ (mV)")
    ax.set_ylabel(r"$P_{\mathrm{sw}}$")
    ax.set_title("Measured vs simulated switching transition, 0.75 ns, P" + r"$\rightarrow$" + "AP")
    ax.set_xlim(790, 1260)
    ax.set_ylim(-0.04, 1.06)
    ax.legend(fontsize=8.8, loc="lower right")
    fig.tight_layout()
    fig.savefig(HERE / "e2_overlay.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Panel 2: relaxation-window sensitivity of the over-drive plateau.
    if relax:
        fig, ax = plt.subplots(figsize=(5.6, 4.4))
        colors = {1400.0: NAVY, 2000.0: AMBER}
        for i_uA in RELAX_CURRENTS_UA:
            ws, ps, ns_ = [], [], []
            for r in sorted(relax, key=lambda d: d["relax_ns"]):
                if i_uA in r["grid_uA"]:
                    j = r["grid_uA"].index(i_uA)
                    ws.append(r["relax_ns"]); ps.append(r["psw"][j]); ns_.append(r["trials"])
            if not ws:
                continue
            ps = np.array(ps); ns_ = np.array(ns_)
            lo, hi = wilson(ps, ns_)
            ax.errorbar(ws, ps, yerr=yerr_pair(ps, lo, hi), fmt="o-", ms=6.5, lw=1.5,
                        capsize=3, color=colors.get(i_uA, CHARCOAL),
                        label=rf"$|I_{{\mathrm{{SOT}}}}|$ = {i_uA:.0f} $\mu$A")
        ax.set_xscale("log")
        ax.set_xticks([3.25, 10, 50])
        ax.set_xticklabels(["3.25", "10", "50"])
        ax.set_xlabel("post-pulse relaxation window (ns)")
        ax.set_ylabel(r"$P_{\mathrm{sw}}$")
        ax.set_title("Over-drive plateau vs judgement window")
        ax.set_ylim(0.5, 1.05)
        ax.legend(fontsize=9.5)
        fig.tight_layout()
        fig.savefig(HERE / "e2_relaxwin.png", bbox_inches="tight", facecolor="white")
        plt.close(fig)

    # ── Console summary ──────────────────────────────────────────────────
    print("=" * 72)
    print("E2 summary (0.75 ns, AP->P, Cayley, self-heating ON)")
    print("=" * 72)
    print(f"pooled dense scan : {len(grid)} points x {n_tot} trials")
    print(f"sim sigmoid       : Vth = {vth_sim*1e3:.1f} mV, k = {k_sim*1e3:.1f} mV, "
          f"beta = {beta_sim:.1f} +/- {1.96*beta_sim_sigma:.1f} V^-1 (R2 = {r2:.3f})")
    print(f"slope ladder      : NB {BETA_NB_ANALYTIC:.1f} -> sLLG {beta_sim:.1f} "
          f"-> measured {MEAS_FIT['beta']:.1f} V^-1")
    print(f"eta_c             : measured {eta_meas:.2f}, macrospin-explained {eta_sim:.2f} "
          f"(log-share {share*100:.0f}%)")
    if over is not None:
        print("overdrive         : " + ", ".join(
            f"{v:.0f}mV:{p:.2f}" for v, p in zip(over["V_mV"], over["psw"])))
    for r in sorted(relax, key=lambda d: d["relax_ns"]):
        print(f"relax window {r['relax_ns']:>5} ns: " + ", ".join(
            f"{u:.0f}uA:{p:.2f}" for u, p in zip(r["grid_uA"], r["psw"])))
    print("Wrote e2_results.json, e2_overlay.png" + (", e2_relaxwin.png" if relax else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dense-shard", type=int, default=None)
    ap.add_argument("--trials", type=int, default=125)
    ap.add_argument("--overdrive", action="store_true")
    ap.add_argument("--relaxwin", type=float, default=None)
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        out = run_scan(DENSE_GRID_UA[::6], trials=4, seed=BASE_SEED + 2900)
        print(json.dumps({k: out[k] for k in ("grid_uA", "psw", "wall_s")}, indent=2))
        print(f"[smoke] OK, {out['wall_s']/(len(DENSE_GRID_UA[::6])*4):.2f} s/trial")
        return
    if args.dense_shard is not None:
        k = args.dense_shard
        print(f"[E2] dense shard {k}: {len(DENSE_GRID_UA)} pts x {args.trials} trials", flush=True)
        out = run_scan(DENSE_GRID_UA, trials=args.trials, seed=BASE_SEED + 2000 + k)
        p = HERE / f"e2_dense_shard{k}.json"
        p.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[E2] wrote {p.name} (wall {out['wall_s']/60:.1f} min)", flush=True)
        return
    if args.overdrive:
        print(f"[E2] overdrive: {len(OVERDRIVE_UA)} pts x 200 trials", flush=True)
        out = run_scan(OVERDRIVE_UA, trials=200, seed=BASE_SEED + 2100)
        (HERE / "e2_overdrive.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[E2] wrote e2_overdrive.json (wall {out['wall_s']/60:.1f} min)", flush=True)
        return
    if args.relaxwin is not None:
        w = args.relaxwin
        trials = 100 if w >= 50 else 200
        print(f"[E2] relaxwin {w} ns: {len(RELAX_CURRENTS_UA)} currents x {trials} trials", flush=True)
        out = run_scan(RELAX_CURRENTS_UA, trials=trials,
                       seed=BASE_SEED + 2200 + int(w * 10), relax_ns=w)
        p = HERE / f"e2_relaxwin_{w:g}.json"
        p.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[E2] wrote {p.name} (wall {out['wall_s']/60:.1f} min)", flush=True)
        return
    if args.analyze:
        analyze()
        return
    ap.error("choose one of --dense-shard / --overdrive / --relaxwin / --analyze / --smoke")


if __name__ == "__main__":
    main()
