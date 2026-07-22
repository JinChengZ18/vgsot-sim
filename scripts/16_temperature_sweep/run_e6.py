"""E6 — ambient-temperature sweep of the probabilistic write interface.

The chapter builds the full temperature machinery (Bloch M_s(T),
Callen-Callen K_i(T), eta(T), self-heating feedback) but exercises it only
at the single 300 K workpoint (fig 2.6/2.11 self-heating comparisons).
This experiment sweeps the ambient temperature T0 in {250, 275, 300, 325,
350, 400} K with self-heating ON and extracts the interface's temperature
coefficients: dV_th/dT (primary, from the raw 50% crossing) and the
per-temperature logistic slope beta_s (secondary; wide CI at 200
trials/point). These are the drift numbers downstream chapters need for
robustness arguments about the Bernoulli interface.

Modes:
  python run_e6.py --temp 325       # one shard per temperature -> e6_t325.json
  python run_e6.py --analyze        # -> e6_results.json + panels
  python run_e6.py --smoke
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

HERE = Path(__file__).resolve().parent
BASE_SEED = 20260719

PULSE_NS, RELAX_NS = 0.75, 3.25
ITH0_UA = 1160.0
TEMPS_K = (250.0, 275.0, 300.0, 325.0, 350.0, 400.0)
MULT = np.array([0.70, 0.78, 0.85, 0.91, 0.96, 1.00, 1.04, 1.09, 1.15, 1.22, 1.30])


def run_temp(T_K: float, trials: int):
    cc = PhysicalConstantsConfig()
    mid1 = int(round(PULSE_NS * 1e-9 / cc.t_step))
    end = mid1 + int(round(RELAX_NS * 1e-9 / cc.t_step))
    grid_uA = np.round(ITH0_UA * MULT)
    cfg = SerSotNoVcmaThermalConfig(
        i_sot_list=tuple(-grid_uA * 1e-6), trials=trials,
        sim_start_step=1, sim_mid1_step=mid1, sim_end_step=end,
        pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
        target_mz=1.0, failure_tol=0.2, constants=cc,
    )
    t0 = time.time()
    res = ser_sot_no_vcma_thermal(
        cfg, show_progress=False, enable_self_heating=True, T_ambient_K=T_K,
        seed=BASE_SEED + 6000 + int(T_K), rng_mode="generator", integrator="cayley",
    )
    return dict(T_K=T_K, grid_uA=grid_uA.tolist(),
                V_mV=(grid_uA * 1e-6 * cc.R_W * 1e3).tolist(),
                psw=res.psw.tolist(), trials=trials, R_W=cc.R_W,
                wall_s=time.time() - t0,
                timestamp=datetime.now().isoformat(timespec="seconds"))


def crossing(V, P, trials):
    V = np.asarray(V, float); P = np.asarray(P, float)
    o = np.argsort(V); V, P = V[o], P[o]
    for i in range(len(V) - 1):
        if (P[i] - 0.5) * (P[i + 1] - 0.5) <= 0 and P[i] != P[i + 1]:
            slope = (P[i + 1] - P[i]) / (V[i + 1] - V[i])
            return float(V[i] + (0.5 - P[i]) / slope), float(np.sqrt(0.25 / trials) / abs(slope))
    return float("nan"), float("nan")


def analyze():
    from scipy.optimize import curve_fit
    from vgsot_sim.analysis.sigmoid_fit import sigmoid4p
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = []
    for f in sorted(glob.glob(str(HERE / "e6_t*.json"))):
        runs.append(json.loads(Path(f).read_text(encoding="utf-8")))
    runs.sort(key=lambda d: d["T_K"])
    if len(runs) < 3:
        raise SystemExit("Need >=3 temperature shards.")

    T = np.array([d["T_K"] for d in runs])
    vth = np.array([crossing(d["V_mV"], d["psw"], d["trials"])[0] for d in runs])
    vth_sig = np.array([crossing(d["V_mV"], d["psw"], d["trials"])[1] for d in runs])

    betas, beta_sigs = [], []
    for d in runs:
        Vv = np.array(d["V_mV"]) / 1e3
        P = np.array(d["psw"])
        try:
            p0 = (0.0, 1.0, float(Vv[np.argmin(np.abs(P - 0.5))]), 0.03)
            bounds = ([-0.05, 0.3, Vv.min() - 0.05, 1e-3], [0.30, 1.10, Vv.max() + 0.05, 0.3])
            popt, pcov = curve_fit(sigmoid4p, Vv, P, p0=p0, bounds=bounds, maxfev=30000)
            betas.append(1.0 / popt[3])
            beta_sigs.append(float(np.sqrt(np.diag(pcov))[3] / popt[3] ** 2))
        except Exception:
            betas.append(float("nan")); beta_sigs.append(float("nan"))
    betas = np.array(betas); beta_sigs = np.array(beta_sigs)

    ok = np.isfinite(vth)
    slope, intercept = np.polyfit(T[ok], vth[ok], 1)   # mV per K

    summary = dict(
        T_K=T.tolist(), vth_mV=vth.tolist(), vth_sigma_mV=vth_sig.tolist(),
        beta_V=betas.tolist(), beta_sigma_V=beta_sigs.tolist(),
        dVth_dT_mV_per_K=float(slope), vth_fit_intercept_mV=float(intercept),
        trials_per_point=runs[0]["trials"],
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    (HERE / "e6_results.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    CHARCOAL, NAVY, CRIMSON, TEAL = "#2B2B2B", "#1F5FA8", "#A82038", "#1A6B5A"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans"],
        "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13.5,
        "mathtext.fontset": "stix", "axes.linewidth": 0.9,
        "axes.edgecolor": CHARCOAL, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "figure.dpi": 150, "savefig.dpi": 300,
    })

    fig, ax = plt.subplots(figsize=(6.2, 4.7))
    ax.errorbar(T[ok], vth[ok], yerr=np.maximum(vth_sig[ok], 1.0), fmt="o", ms=7,
                color=CRIMSON, mfc="white", mew=1.6, capsize=3, lw=1.3,
                label="vgsot-sim 50% threshold")
    tt = np.linspace(240, 410, 100)
    ax.plot(tt, slope * tt + intercept, "--", color=TEAL, lw=1.6,
            label=rf"linear fit  {slope:+.2f} mV/K")
    ax.set_xlabel(r"ambient temperature $T_0$ (K)")
    ax.set_ylabel(r"$|V_{\mathrm{th}}|$ (mV)")
    ax.set_title("Write-threshold temperature drift, self-heating ON")
    ax.legend(fontsize=9.5)
    fig.tight_layout()
    fig.savefig(HERE / "e6_vth_vs_T.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.7))
    cmap = plt.get_cmap("coolwarm")
    for i, d in enumerate(runs):
        c = cmap(i / max(1, len(runs) - 1))
        ax.plot(d["V_mV"], d["psw"], "o-", ms=4.2, lw=1.25, color=c,
                label=rf"$T_0$ = {d['T_K']:.0f} K")
    ax.axhline(0.5, color="gray", lw=0.6, ls=":")
    ax.set_xlabel(r"$|V_{\mathrm{SOT}}|$ (mV)")
    ax.set_ylabel(r"$P_{\mathrm{sw}}$")
    ax.set_title(r"$P_{\mathrm{sw}}(V)$ vs ambient temperature")
    ax.set_ylim(-0.04, 1.06)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "e6_psw_family.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print("=" * 68)
    print("E6 summary (0.75 ns, AP->P, Cayley, self-heating ON)")
    print("=" * 68)
    print(f"{'T0 (K)':>8} {'Vth (mV)':>10} {'sigma':>7} {'beta (V^-1)':>12} {'95%CI':>8}")
    for t, v, s, b, bs in zip(T, vth, vth_sig, betas, beta_sigs):
        print(f"{t:>8.0f} {v:>10.1f} {s:>7.1f} {b:>12.1f} {1.96*bs:>8.1f}")
    print(f"dVth/dT = {slope:+.3f} mV/K over {T[ok].min():.0f}-{T[ok].max():.0f} K")
    print("Wrote e6_results.json, e6_vth_vs_T.png, e6_psw_family.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--temp", type=float, default=None)
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        cc = PhysicalConstantsConfig()
        mid1 = int(round(PULSE_NS * 1e-9 / cc.t_step))
        end = mid1 + int(round(RELAX_NS * 1e-9 / cc.t_step))
        for Tk in (250.0, 400.0):
            cfg = SerSotNoVcmaThermalConfig(
                i_sot_list=(-1160e-6,), trials=8,
                sim_start_step=1, sim_mid1_step=mid1, sim_end_step=end,
                pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
                target_mz=1.0, failure_tol=0.2, constants=cc,
            )
            res = ser_sot_no_vcma_thermal(cfg, show_progress=False,
                                          enable_self_heating=True, T_ambient_K=Tk,
                                          seed=BASE_SEED + 6900, rng_mode="generator",
                                          integrator="cayley")
            print(f"[smoke] T0={Tk:.0f} K @1160 uA -> P_sw={res.psw[0]:.2f}")
        return
    if args.temp is not None:
        print(f"[E6] T0={args.temp:g} K: {len(MULT)} pts x {args.trials} trials", flush=True)
        out = run_temp(args.temp, args.trials)
        p = HERE / f"e6_t{args.temp:g}.json"
        p.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[E6] wrote {p.name} (wall {out['wall_s']/60:.1f} min)", flush=True)
        return
    if args.analyze:
        analyze()
        return
    ap.error("choose one of --temp / --analyze / --smoke")


if __name__ == "__main__":
    main()
