"""E1 — pulse-width transferability of the single-point theta_SH calibration.

The chapter calibrates theta_SH = 0.066 against ONE experimental anchor
(Device A, P->AP, t_w = 0.75 ns, V_th = 894 mV, Cayley integrator,
self-heating ON). This experiment asks whether that single-point
calibration *transfers*: simulate V_th^sim(t_w) at the pulse widths the
experiment actually measured (0.75 / 1 / 2 / 5 ns, plus 1.5 / 3 ns
fill-ins), overlay the experimental log-linear law, then invert the
simulated V_th(ln t_w) with the same Neel-Brown machinery used on the
measured data (analysis.nb_fit) and compare (Delta_pulse, V_c0, b,
tau_ret) against the experimental Table-2.8 values.

A zero-drive free-evolution run at the calibrated (deep-barrier)
parameter set bounds the *actual* retention of the simulated device at
|H_ex| = 50 Oe (model default) and 200 Oe (experimental bias): zero
flips over the observed window shows tau_ret^NB ~ 172 ns is a
linear-barrier extrapolation artifact, not a physical dwell time.

Experimental references (article/chapter02.md, 2.3.2-2.3.3, Device A, P->AP):
  V(t_w) = a - b ln(t_w/ns),  a = 0.795 V, b = 0.175 V  (Table 2.7 precision;
  prose rounds to 0.79/0.18);  sigmoid V_th(0.75 ns) = 894 mV;
  NB inversion at tau0 = 1 ns: Delta = 4.91, V_c0 = 857 mV, tau_ret = 135 ns.

Modes (composable for parallel execution, one process per width):
  python run_e1.py --widths 0.75            # MC for one width -> e1_w0.75.json
  python run_e1.py --freerun-only           # zero-drive runs   -> e1_freerun.json
  python run_e1.py --analyze                # merge e1_w*.json  -> figures + e1_results.json
  python run_e1.py --smoke                  # tiny end-to-end sanity pass
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

from vgsot_sim.configs import PhysicalConstantsConfig, SerSotNoVcmaThermalConfig
from vgsot_sim.ser_cases import ser_sot_no_vcma_thermal
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation
from vgsot_sim.analysis import nb_fit
from vgsot_sim.analysis.sigmoid_fit import sigmoid4p

HERE = Path(__file__).resolve().parent

# ── Experimental references (Device A, P->AP) ────────────────────────────
EXP_A, EXP_B = 0.795, 0.175          # V(t) = a - b ln(t/ns), Table-2.7 precision
EXP_WIDTHS_NS = (0.75, 1.0, 2.0, 5.0)  # pulse widths actually measured
EXP_VTH_SIGMOID_075 = 0.894          # V, sigmoid 50% point at 0.75 ns
EXP_NB = dict(Delta=4.91, Vc0=0.857, tau_ret_ns=135.0, tau0_ns=1.0)

# Simulated 50% point at the calibration anchor (fig 2.11 inset, amber curve).
SIM_ITH_075_UA = 1160.0

RELAX_NS = 3.25                      # same relaxation tail as fig 2.11
MULTIPLIERS = np.array([0.85, 0.90, 0.94, 0.97, 1.00, 1.03, 1.06, 1.10, 1.16])
BASE_SEED = 20260719


def exp_line(tw_ns):
    """Experimental log-linear law |V_th|(t_w), volts."""
    return EXP_A - EXP_B * np.log(np.asarray(tw_ns, dtype=float))


def center_uA(tw_ns: float) -> float:
    """Scan-window center: scale the known simulated 0.75-ns threshold by
    the experimental log-linear shape (assumption only affects window
    placement, not the extracted V_th)."""
    return SIM_ITH_075_UA * float(exp_line(tw_ns) / exp_line(0.75))


def vth_from_curve(V, P):
    """50% crossing, two ways: 4-parameter logistic fit (headline, with
    1-sigma from the covariance) and direct monotone interpolation
    (cross-check). Returns dict."""
    from scipy.optimize import curve_fit

    V = np.asarray(V, float)
    P = np.asarray(P, float)
    order = np.argsort(V)
    V, P = V[order], P[order]

    p0 = (0.0, 1.0, float(V[np.argmin(np.abs(P - 0.5))]), 0.02)
    bounds = ([-0.05, 0.3, V.min() - 0.05, 1e-3],
              [0.30, 1.10, V.max() + 0.05, 0.20])
    popt, pcov = curve_fit(sigmoid4p, V, P, p0=p0, bounds=bounds, maxfev=20000)
    y0, L, vth_fit, k = (float(x) for x in popt)
    vth_sigma = float(np.sqrt(np.diag(pcov))[2])
    # The logistic midpoint parameter is where P = y0 + L/2; report also the
    # literal P = 0.5 crossing of the fitted curve for transparency.
    vth_p50_fit = float(vth_fit + k * np.log(L / (0.5 - y0) - 1.0)) \
        if (0.0 < 0.5 - y0 < L) else float("nan")

    # Direct interpolation of the raw points across 0.5.
    vth_interp = float("nan")
    for i in range(len(V) - 1):
        if (P[i] - 0.5) * (P[i + 1] - 0.5) <= 0 and P[i] != P[i + 1]:
            vth_interp = float(V[i] + (0.5 - P[i]) * (V[i + 1] - V[i]) / (P[i + 1] - P[i]))
            break
    return dict(vth_fit=vth_fit, vth_sigma=vth_sigma, vth_p50_fit=vth_p50_fit,
                vth_interp=vth_interp, y0=y0, L=L, k=k)


def run_width(tw_ns: float, trials: int, multipliers=MULTIPLIERS, seed_offset: int = 0,
              center_override_uA: float | None = None):
    """Monte-Carlo P_sw(V) transition scan at one pulse width.

    `center_override_uA` re-centers the scan window; needed at 3/5 ns where
    the simulated threshold sits well above the experimental-line-scaled
    prediction (the transfer deviation under test) and the default window
    misses the transition."""
    cc = PhysicalConstantsConfig()
    mid1 = int(round(tw_ns * 1e-9 / cc.t_step))
    end = mid1 + int(round(RELAX_NS * 1e-9 / cc.t_step))
    grid_uA = np.round((center_override_uA or center_uA(tw_ns)) * multipliers)
    cfg = SerSotNoVcmaThermalConfig(
        i_sot_list=tuple(-grid_uA * 1e-6), trials=trials,
        sim_start_step=1, sim_mid1_step=mid1, sim_end_step=end,
        pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
        target_mz=1.0, failure_tol=0.2, constants=cc,
    )
    t0 = time.time()
    res = ser_sot_no_vcma_thermal(
        cfg, show_progress=False, enable_self_heating=True, T_ambient_K=300.0,
        seed=BASE_SEED + seed_offset, rng_mode="generator", integrator="cayley",
    )
    wall_s = time.time() - t0
    V = grid_uA * 1e-6 * cc.R_W          # |V_SOT| in volts
    psw = res.psw
    out = dict(
        tw_ns=tw_ns, grid_uA=grid_uA.tolist(), V=V.tolist(), psw=psw.tolist(),
        trials=trials, wall_s=wall_s, R_W=cc.R_W, theta_SH=cc.theta_SH,
        t_step=cc.t_step, relax_ns=RELAX_NS, integrator="cayley",
        rng_mode="generator", self_heating=True, seed=BASE_SEED + seed_offset,
        vth=vth_from_curve(V, psw),
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    return out


def run_freerun(seg_ns: float, n_seeds: int, fields_Oe=(50.0, 200.0)):
    """Zero-drive free evolution at the calibrated deep-barrier parameter
    set; count +/-0.5-hysteresis m_z crossings (expected: zero)."""
    results = []
    for H_Oe in fields_Oe:
        cc = replace(PhysicalConstantsConfig(), t_step=4e-12,
                     h_ex_y=-H_Oe * 1000.0 / (4 * np.pi))
        steps = int(round(seg_ns * 1e-9 / cc.t_step))
        for s in range(n_seeds):
            rng = np.random.default_rng(BASE_SEED + 1000 + int(H_Oe) * 10 + s)
            t0 = time.time()
            res = run_piecewise_direct_excitation(
                sim_start_step=1, sim_mid1_step=steps, sim_mid2_step=steps,
                sim_end_step=steps, pap=1,
                v_mtj_stage1=0.0, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
                i_sot_stage1=0.0, i_sot_stage2=0.0, i_sot_stage3=0.0,
                estt_stage1=0, esot_stage1=0, estt_stage2=0, esot_stage2=0,
                estt_stage3=0, esot_stage3=0,
                vnv=0, non=1, r_sot_fl_dl=0.83, show_progress=False,
                constants=cc, enable_self_heating=False, integrator="cayley",
                rng=rng,
            )
            mz = np.asarray(res.mz)
            # Hysteresis-threshold flip counting (same |m_z| > 0.5 convention
            # as the 2.4.2 dwell extraction).
            state, crossings = -1, 0
            for v in mz[:: max(1, steps // 500_000)]:
                if state == -1 and v > 0.5:
                    state, crossings = +1, crossings + 1
                elif state == +1 and v < -0.5:
                    state, crossings = -1, crossings + 1
            results.append(dict(
                H_Oe=H_Oe, seed_idx=s, seg_ns=seg_ns, t_step=cc.t_step,
                crossings=crossings, max_mz=float(mz.max()), min_mz=float(mz.min()),
                wall_s=time.time() - t0,
            ))
            print(f"  freerun H={H_Oe:.0f} Oe seed {s}: {seg_ns:.0f} ns, "
                  f"crossings={crossings}, max_mz={mz.max():+.3f}", flush=True)
    return results


def analyze(out_json: Path, fig_prefix: str):
    """Merge per-width JSONs + freerun JSON -> NB inversion, figures, summary."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    raw = []
    for f in sorted(glob.glob(str(HERE / "e1_w*.json"))):
        with open(f, encoding="utf-8") as fh:
            raw.append(json.load(fh))
    # Merge files sharing a pulse width (base window + _hi supplements).
    from collections import defaultdict
    groups = defaultdict(list)
    for d in raw:
        groups[d["tw_ns"]].append(d)
    widths = []
    for t in sorted(groups):
        ds = groups[t]
        V = np.concatenate([np.asarray(d["V"], float) for d in ds])
        P = np.concatenate([np.asarray(d["psw"], float) for d in ds])
        o = np.argsort(V)
        widths.append(dict(tw_ns=t, V=V[o].tolist(), psw=P[o].tolist(),
                           trials=ds[0]["trials"], vth=ds[0]["vth"],
                           n_files=len(ds)))
    if len(widths) < 3:
        raise SystemExit(f"Need >=3 widths for the log-linear fit, found {len(widths)}.")

    def crossing_and_sigma(V, P, n_trials):
        """Raw P = 0.5 crossing + binomial-noise-propagated 1-sigma.

        This is the like-for-like counterpart of the experimental threshold
        (hysteresis jump / sigmoid 50% point). The 4-parameter logistic
        midpoint parameter is NOT comparable here: the over-drive
        back-hopping plateau compresses its L below 1 and drags the
        midpoint parameter tens of mV below the actual 50% crossing.
        """
        V = np.asarray(V, float); P = np.asarray(P, float)
        o = np.argsort(V); V, P = V[o], P[o]
        for i in range(len(V) - 1):
            if (P[i] - 0.5) * (P[i + 1] - 0.5) <= 0 and P[i] != P[i + 1]:
                slope = (P[i + 1] - P[i]) / (V[i + 1] - V[i])
                vc = V[i] + (0.5 - P[i]) / slope
                sig = np.sqrt(0.25 / n_trials) / abs(slope)
                return float(vc), float(sig)
        return float("nan"), float("nan")

    tw = np.array([d["tw_ns"] for d in widths])
    cross = [crossing_and_sigma(d["V"], d["psw"], d["trials"]) for d in widths]
    vth = np.array([c[0] for c in cross])          # headline: raw 50% crossing
    vth_sig = np.array([c[1] for c in cross])
    vth_itp = np.array([d["vth"]["vth_fit"] for d in widths])  # secondary: logistic midpoint

    sim_fit = nb_fit.fit_direction(tw, vth, tau0=1.0)
    dev_line = vth - exp_line(tw)

    freerun = []
    fr_path = HERE / "e1_freerun.json"
    if fr_path.exists():
        freerun = json.loads(fr_path.read_text(encoding="utf-8"))

    summary = dict(
        widths=widths, freerun=freerun,
        sim_nb=dict(a=sim_fit.a, b=sim_fit.b, Delta=sim_fit.Delta,
                    Vc0=sim_fit.Vc0, tau_ret_ns=sim_fit.tau_ret_ns, tau0_ns=1.0),
        exp_nb=EXP_NB, exp_line=dict(a=EXP_A, b=EXP_B),
        vth_definition="raw P=0.5 crossing (like-for-like with experimental threshold); "
                       "logistic midpoint kept as secondary",
        vth_crossing_mV=(vth * 1e3).tolist(),
        vth_logistic_mid_mV=(vth_itp * 1e3).tolist(),
        deviation_vs_exp_line_mV=(dev_line * 1e3).tolist(),
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # ── Figure style (repo academic convention: Arial, English, no baked
    #    panel letters / figure numbers) ─────────────────────────────────
    CHARCOAL, NAVY, CRIMSON, TEAL, AMBER = "#2B2B2B", "#1F5FA8", "#A82038", "#1A6B5A", "#C47A00"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans"],
        "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13.5,
        "mathtext.fontset": "stix", "axes.linewidth": 0.9,
        "axes.edgecolor": CHARCOAL, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "figure.dpi": 150, "savefig.dpi": 300,
    })

    # Panel 1: V_th vs t_w (log axis) — sim points vs experimental law.
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    tt = np.geomspace(0.6, 6.5, 200)
    ax.plot(tt, exp_line(tt) * 1e3, "--", color=TEAL, lw=1.8,
            label=r"experiment log-linear fit  $V=a-b\,\ln t_w$")
    ax.plot(EXP_WIDTHS_NS, exp_line(EXP_WIDTHS_NS) * 1e3, "s", ms=8, mfc="white",
            mec=TEAL, mew=1.6, label="experimental pulse widths")
    # The simulator was calibrated to the sigmoid 50% anchor (894 mV), which
    # itself sits ~6% above the hysteresis-jump log-linear line at 0.75 ns —
    # mark it so definitional offset is not misread as model deviation.
    ax.plot([0.75], [EXP_VTH_SIGMOID_075 * 1e3], "*", ms=13, color=AMBER,
            mec=CHARCOAL, mew=0.6, zorder=5,
            label="sigmoid 50% anchor (calibration target)")
    ax.errorbar(tw, vth * 1e3, yerr=np.maximum(vth_sig, 1e-4) * 1e3, fmt="o",
                ms=7, color=CRIMSON, mfc="white", mew=1.6, capsize=3, lw=1.4,
                label=r"vgsot-sim ($\theta_{\mathrm{SH}}$ calibrated at 0.75 ns only)")
    aa, bb = sim_fit.a, sim_fit.b
    ax.plot(tt, (aa - bb * np.log(tt)) * 1e3, "-", color=CRIMSON, lw=1.2, alpha=0.7)
    ax.set_xscale("log")
    ax.set_xticks([0.75, 1, 1.5, 2, 3, 5])
    ax.set_xticklabels(["0.75", "1", "1.5", "2", "3", "5"])
    ax.set_xlabel(r"pulse width $t_w$ (ns)")
    ax.set_ylabel(r"$|V_{\mathrm{th}}|$ (mV)")
    ax.set_title("Threshold-voltage transfer across pulse widths")
    txt = (f"sim:  b = {bb*1e3:.0f} mV,  $\\Delta$ = {sim_fit.Delta:.2f},  "
           f"$V_{{c0}}$ = {sim_fit.Vc0*1e3:.0f} mV\n"
           f"exp:  b = {EXP_B*1e3:.0f} mV,  $\\Delta$ = {EXP_NB['Delta']:.2f},  "
           f"$V_{{c0}}$ = {EXP_NB['Vc0']*1e3:.0f} mV   ($\\tau_0$ = 1 ns)")
    ax.text(0.03, 0.06, txt, transform=ax.transAxes, fontsize=10,
            va="bottom", ha="left",
            bbox=dict(fc="white", ec="#CCCCCC", alpha=0.9))
    ax.legend(fontsize=9.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(HERE / f"{fig_prefix}_vth_vs_tw.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Panel 2: the P_sw(V) family behind the thresholds.
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    cmap = plt.get_cmap("viridis")
    for i, d in enumerate(widths):
        c = cmap(i / max(1, len(widths) - 1) * 0.85)
        ax.plot(np.array(d["V"]) * 1e3, d["psw"], "o-", ms=4.5, lw=1.3, color=c,
                label=rf"$t_w$ = {d['tw_ns']:g} ns")
    ax.axhline(0.5, color="gray", lw=0.6, ls=":")
    ax.set_xlabel(r"$|V_{\mathrm{SOT}}|$ (mV)")
    ax.set_ylabel(r"$P_{\mathrm{sw}}$")
    ax.set_title(r"Simulated $P_{\mathrm{sw}}(V)$ transitions, self-heating ON")
    ax.set_ylim(-0.04, 1.06)
    ax.legend(fontsize=9.5)
    fig.tight_layout()
    fig.savefig(HERE / f"{fig_prefix}_psw_family.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # ── Console summary ──────────────────────────────────────────────────
    print("=" * 74)
    print("E1 summary  (theta_SH = 0.066, Cayley, self-heating ON, P->AP)")
    print("=" * 74)
    print(f"{'t_w (ns)':>9} {'Vth_sim (mV)':>13} {'Vth_expline (mV)':>17} {'dev':>7}")
    for t, v, dv in zip(tw, vth, dev_line):
        print(f"{t:>9g} {v*1e3:>13.1f} {exp_line(t)*1e3:>17.1f} {dv/exp_line(t)*100:>6.1f}%")
    print(f"\nlog-linear fit  : a = {sim_fit.a*1e3:.0f} mV, b = {sim_fit.b*1e3:.0f} mV "
          f"(exp {EXP_A*1e3:.0f}/{EXP_B*1e3:.0f} mV)")
    print(f"NB inversion    : Delta = {sim_fit.Delta:.2f}, Vc0 = {sim_fit.Vc0*1e3:.0f} mV, "
          f"tau_ret = {sim_fit.tau_ret_ns:.0f} ns  (exp {EXP_NB['Delta']:.2f} / "
          f"{EXP_NB['Vc0']*1e3:.0f} mV / {EXP_NB['tau_ret_ns']:.0f} ns, tau0 = 1 ns)")
    if freerun:
        tot = {}
        for r in freerun:
            k = r["H_Oe"]
            tot.setdefault(k, [0.0, 0])
            tot[k][0] += r["seg_ns"]
            tot[k][1] += r["crossings"]
        for H, (ns, cr) in sorted(tot.items()):
            exp_flips = ns / EXP_NB["tau_ret_ns"]
            print(f"free evolution  : H = {H:.0f} Oe, {ns/1e3:.1f} us observed, "
                  f"{cr} flips (tau_ret = 135 ns would predict ~{exp_flips:.0f})")
    print(f"\nWrote {out_json.name}, {fig_prefix}_vth_vs_tw.png, {fig_prefix}_psw_family.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--widths", type=str, default=None,
                    help="comma-separated pulse widths in ns to run MC for")
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--center-ua", type=float, default=None,
                    help="override scan-window center (uA) for supplementary scans")
    ap.add_argument("--suffix", type=str, default="",
                    help="output filename suffix, e.g. _hi for a re-centered supplement")
    ap.add_argument("--freerun-only", action="store_true")
    ap.add_argument("--freerun-seg-ns", type=float, default=2000.0)
    ap.add_argument("--freerun-seeds", type=int, default=5)
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--smoke", action="store_true",
                    help="tiny sanity pass: 1 width x 3 points x 8 trials + 40 ns freerun")
    args = ap.parse_args()

    if args.smoke:
        print("[smoke] MC 0.75 ns, 3 points, 8 trials")
        out = run_width(0.75, trials=8, multipliers=np.array([0.94, 1.00, 1.06]))
        (HERE / "e1_smoke.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(json.dumps({k: out[k] for k in ("grid_uA", "psw", "wall_s")}, indent=2))
        print("[smoke] freerun 40 ns x 1 seed x 2 fields")
        fr = run_freerun(seg_ns=40.0, n_seeds=1)
        print(f"[smoke] OK; MC wall {out['wall_s']:.1f} s -> "
              f"{out['wall_s']/ (3*8):.2f} s/trial at 4 ns span")
        return

    if args.freerun_only:
        fr = run_freerun(seg_ns=args.freerun_seg_ns, n_seeds=args.freerun_seeds)
        (HERE / "e1_freerun.json").write_text(json.dumps(fr, indent=2), encoding="utf-8")
        print(f"Wrote e1_freerun.json ({len(fr)} segments)")
        return

    if args.analyze:
        analyze(HERE / "e1_results.json", "e1")
        return

    if not args.widths:
        ap.error("choose one of --widths / --freerun-only / --analyze / --smoke")
    for i, w in enumerate(float(x) for x in args.widths.split(",")):
        print(f"[E1] MC width {w:g} ns, {args.trials} trials x {len(MULTIPLIERS)} points",
              flush=True)
        out = run_width(w, trials=args.trials,
                        seed_offset=int(round(w * 100)) + (7 if args.center_ua else 0),
                        center_override_uA=args.center_ua)
        p = HERE / f"e1_w{w:g}{args.suffix}.json"
        p.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[E1] wrote {p.name}  (wall {out['wall_s']/60:.1f} min, "
              f"Vth_fit = {out['vth']['vth_fit']*1e3:.1f} mV)", flush=True)


if __name__ == "__main__":
    main()
