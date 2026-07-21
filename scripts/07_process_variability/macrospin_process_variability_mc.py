"""Macrospin Psw Monte Carlo with PDK process-mismatch samples.

The NB variability figure averages over Delta only.  This script pushes the
same PDK/Brinkman mismatch budget into the macrospin solver so the magnetic
parameters and the SOT voltage-to-current conversion vary per sampled device.
The axes stay in terminal SOT voltage: the main sweep mirrors Chapter02_local_11's broad
current span through the nominal R_W, while each sampled device keeps that
voltage fixed and therefore sees its own process-shifted I_SOT.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from vgsot_sim.analysis.variability import (
    PDKBudgetInputs,
    cv_delta_budget,
    sample_macrospin_process_constants,
)
from vgsot_sim.configs import PhysicalConstantsConfig, SerSotNoVcmaThermalConfig
from vgsot_sim.ser_cases import ser_sot_no_vcma_thermal


PULSE_NS = 0.75
TOTAL_NS = 4.0
MEASURED_VTH_MV = 894.0

# Follow Chapter02_local_11: a broad-spectrum main sweep plus a dense threshold
# sweep.  These nominal current landmarks are converted to terminal voltages
# once using the baseline channel resistance before process MC is evaluated.
WIDE_I_SOT_UA = np.array([300, 600, 900, 1100, 1300, 1500, 1800,
                          2200, 2800, 3500], dtype=float)
INSET_I_SOT_UA = np.array([800, 900, 1000, 1080, 1120, 1140, 1160, 1180,
                           1220, 1280, 1350, 1400], dtype=float)

PDK_INPUTS = PDKBudgetInputs(
    CV_RP=0.07,
    CV_RSOT=0.07,
    CV_TMR=0.04,
    CV_TOX=0.003,
    CV_PHI=0.005,
    CV_TF=0.003,
    CV_MS=0.02,
    T_OX=1.0,
    PHI_BAR=0.6,
    KAPPA=10.25,
)
PROCESS_Z_KEYS = ("D", "tf", "Ms", "Ki_proxy", "Rsot", "tox", "phi")

NAVY = "#1F5FA8"
CRIMSON = "#A82038"
TEAL = "#1A6B5A"
THU_PALE = "#C99FD4"
THU_GRID = "#DDD0E8"
CHARCOAL = "#2B2B2B"
NEAR_WHITE = "#FFFFFF"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Arial", "Liberation Serif"],
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "axes.titlepad": 8,
    "legend.fontsize": 10,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "mathtext.fontset": "stix",
    "axes.linewidth": 0.9,
    "axes.edgecolor": CHARCOAL,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "legend.frameon": True,
    "legend.framealpha": 0.94,
    "legend.edgecolor": THU_PALE,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--devices", type=int, default=16,
                        help="number of D2D macrospin parameter samples")
    parser.add_argument("--thermal-trials", type=int, default=32,
                        help="thermal MC trials per sampled device and voltage")
    parser.add_argument("--seed", type=int, default=7,
                        help="master seed for device mismatch and thermal noise")
    parser.add_argument("--integrator", choices=("euler_spherical", "cayley"),
                        default="cayley",
                        help="macrospin time-stepper; Cayley is the calibrated kernel "
                             "(theta_SH = 0.066 was recalibrated against it)")
    parser.add_argument("--disable-self-heating", action="store_true",
                        help="run without the self-heating feedback loop")
    return parser.parse_args()


def wilson_interval(p, n, z=1.96):
    p = np.clip(np.asarray(p, dtype=float), 0.0, 1.0)
    denom = 1.0 + z ** 2 / n
    centre = (p + z ** 2 / (2.0 * n)) / denom
    half = z * np.sqrt(p * (1.0 - p) / n + z ** 2 / (4.0 * n ** 2)) / denom
    return centre - half, centre + half


def run_psw_curve(constants, voltages, trials, *, seed, self_heating, integrator):
    pulse_end = int(round(PULSE_NS * 1e-9 / constants.t_step))
    sim_end = int(round(TOTAL_NS * 1e-9 / constants.t_step))
    i_sot_list = tuple(-np.asarray(voltages, dtype=float) / constants.R_W)
    cfg = SerSotNoVcmaThermalConfig(
        i_sot_list=i_sot_list,
        trials=trials,
        sim_start_step=1,
        sim_mid1_step=pulse_end,
        sim_end_step=sim_end,
        pap=1,
        non=1,
        vnv=0,
        v_mtj=0.0,
        r_sot_fl_dl=0.83,
        target_mz=1.0,
        failure_tol=0.2,
        constants=constants,
    )
    return ser_sot_no_vcma_thermal(
        cfg,
        show_progress=False,
        enable_self_heating=self_heating,
        T_ambient_K=300.0,
        seed=seed,
        rng_mode="generator",
        integrator=integrator,
    ).psw


def make_process_samples(base, n_devices, rng):
    """Use antithetic z-score pairs so small D2D populations stay centered."""
    z_rows = []
    for _ in range(n_devices // 2):
        z = rng.normal(size=len(PROCESS_Z_KEYS))
        z_rows.extend((z, -z))
    if n_devices % 2:
        z_rows.append(np.zeros(len(PROCESS_Z_KEYS)))
    return [
        sample_macrospin_process_constants(
            base,
            PDK_INPUTS,
            rng=rng,
            z_scores=dict(zip(PROCESS_Z_KEYS, z_row)),
        )
        for z_row in z_rows
    ]


def print_scale_summary(samples):
    keys = ("D", "tf", "Ms", "Ki_proxy", "Rsot", "RA")
    print("  Sampled scale factors (mean +/- std):")
    for key in keys:
        values = np.array([sample.scales[key] for sample in samples], dtype=float)
        ddof = 1 if values.size > 1 else 0
        print(f"    {key:>8s}: {values.mean():.4f} +/- {values.std(ddof=ddof):.4f}")
    r_w = np.array([sample.constants.R_W for sample in samples], dtype=float)
    ddof = 1 if r_w.size > 1 else 0
    print(f"    {'R_W (ohm)':>8s}: {r_w.mean():.2f} +/- {r_w.std(ddof=ddof):.2f}")


def nominal_voltage_sweep(i_sot_ua, constants):
    """Convert Chapter02_local_11-style |I_SOT| landmarks into baseline |V_SOT| values."""
    return np.asarray(i_sot_ua, dtype=float) * 1e-6 * constants.R_W


def split_sweeps(values, n_wide):
    return values[:n_wide], values[n_wide:]


def plot_band_curve(ax, x_mv, psw, n_trials, *, color, fmt, label,
                    markersize=6, lw=1.9):
    lo, hi = wilson_interval(psw, n_trials)
    ax.fill_between(x_mv, lo, hi, color=color, alpha=0.14, lw=0)
    ax.plot(x_mv, psw, fmt, color=color, lw=lw, markersize=markersize,
            markerfacecolor="white", markeredgewidth=1.5, label=label)


def plot_curves(wide_voltages, baseline_wide, mean_wide,
                inset_voltages, baseline_inset, mean_inset, args, budget):
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    wide_mv = wide_voltages * 1e3
    inset_mv = inset_voltages * 1e3
    pooled_trials = args.devices * args.thermal_trials

    plot_band_curve(
        ax, wide_mv, baseline_wide, pooled_trials,
        color=NAVY, fmt="o-", label="nominal macrospin",
    )
    plot_band_curve(
        ax, wide_mv, mean_wide, pooled_trials,
        color=CRIMSON, fmt="s--",
        label=rf"PDK mismatch mean "
              rf"($\mathrm{{CV}}_\Delta={budget.CV_Delta*100:.1f}\%$)",
    )

    ax.axvline(MEASURED_VTH_MV, color=TEAL, lw=1.4, ls=(0, (4, 2)),
               label=r"measured $V_{\rm th}=894$ mV")
    ax.axhline(0.5, color=CHARCOAL, lw=0.6, ls=":", alpha=0.65)
    ax.set_xlabel(r"$|V_{\rm SOT}|$ (mV)")
    ax.set_ylabel(r"$P_{\rm sw}\equiv P(\mathrm{switch})$")
    ax.set_title(
        rf"vgsot-sim process-variability Monte-Carlo $P_{{\rm sw}}$ sweep, "
        rf"$t_p={PULSE_NS:g}\,\mathrm{{ns}}$, $V_{{\rm MTJ}}=0$"
    )
    ax.set_xlim(0, wide_mv.max() + 40)
    ax.set_ylim(-0.05, 1.08)
    ax.grid(True, color=THU_GRID, ls="--", lw=0.55, alpha=0.75)
    ax.legend(loc="upper left", framealpha=0.95)

    inset = ax.inset_axes([0.43, 0.10, 0.55, 0.50])
    plot_band_curve(
        inset, inset_mv, baseline_inset, pooled_trials,
        color=NAVY, fmt="o-", label=None, markersize=5.3, lw=1.7,
    )
    plot_band_curve(
        inset, inset_mv, mean_inset, pooled_trials,
        color=CRIMSON, fmt="s--", label=None, markersize=5.1, lw=1.7,
    )
    inset.axvline(MEASURED_VTH_MV, color=TEAL, lw=1.2, ls=(0, (4, 2)))
    inset.axhline(0.5, color=CHARCOAL, lw=0.55, ls=":", alpha=0.65)
    inset.set_xlabel(r"$|V_{\rm SOT}|$ (mV)", fontsize=10, labelpad=2)
    inset.set_ylabel(r"$P_{\rm sw}$", fontsize=10, labelpad=2)
    inset.set_title("Threshold-region sigmoid (zoom)", fontsize=10.5, pad=4)
    inset.set_xlim(inset_mv.min() - 12, inset_mv.max() + 12)
    inset.set_ylim(-0.05, 1.08)
    inset.tick_params(axis="both", labelsize=9)
    inset.grid(True, color=THU_GRID, ls="--", lw=0.45, alpha=0.65)

    fig.tight_layout()
    out_path = Path(__file__).resolve().parent / "Chapter02_local_17.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor=NEAR_WHITE)
    plt.close(fig)
    print(f"Saved {out_path}")


def level_crossing(voltages_v, psw, level=0.5):
    """Raw P_sw = `level` crossing (V) by linear interpolation, plus the
    binomial-propagated 1-sigma from the local slope. Fit-free, so it stays
    meaningful for the wafer-mean curve, which does not saturate inside the
    inset window."""
    V = np.asarray(voltages_v, dtype=float)
    P = np.asarray(psw, dtype=float)
    order = np.argsort(V)
    V, P = V[order], P[order]
    for i in range(len(V) - 1):
        if (P[i] - level) * (P[i + 1] - level) <= 0 and P[i] != P[i + 1]:
            slope = (P[i + 1] - P[i]) / (V[i + 1] - V[i])
            return float(V[i] + (level - P[i]) / slope), float(abs(slope))
    return float("nan"), float("nan")


def logistic_slope(voltages_v, psw, y0_min=0.0):
    """4-parameter logistic slope beta = 1/k (V^-1) with a 95% half-width.

    `y0_min` is exposed because the wafer-mean curve does not saturate inside
    the window: with a negative baseline allowed, y0 pins to its bound and
    beta becomes bound-controlled (4.7-9.1 V^-1 across the usual conventions).
    The default y0_min = 0 is the physical floor (P_sw >= 0) and is the only
    setting for which the fit leaves every parameter interior.
    """
    from scipy.optimize import curve_fit
    from vgsot_sim.analysis.sigmoid_fit import sigmoid4p

    V = np.asarray(voltages_v, dtype=float)
    P = np.asarray(psw, dtype=float)
    p0 = (0.0, 0.75, float(V[np.argmin(np.abs(P - 0.4))]), 0.03)
    bounds = ([y0_min, 0.2, V.min() - 0.1, 1e-3], [0.30, 1.10, V.max() + 0.3, 0.6])
    popt, pcov = curve_fit(sigmoid4p, V, P, p0=p0, bounds=bounds, maxfev=60000)
    resid = P - sigmoid4p(V, *popt)
    r2 = 1.0 - float(np.sum(resid ** 2) / np.sum((P - P.mean()) ** 2))
    beta = 1.0 / popt[3]
    return dict(beta=float(beta), beta_ci95=float(1.96 * np.sqrt(np.diag(pcov))[3] / popt[3] ** 2),
                y0=float(popt[0]), L=float(popt[1]), Vth_mV=float(popt[2] * 1e3),
                k_mV=float(popt[3] * 1e3), R2=r2, y0_min=y0_min)


def summarize_broadening(inset_voltages, baseline_inset, mean_inset,
                         device_inset, pooled_trials, args):
    """Quantify D2D broadening without leaning on a saturation-dependent fit.

    Headline metric is the fit-free 0.25 -> 0.50 voltage span; the logistic
    slopes are reported alongside with their bound sensitivity, and the
    per-device crossings give the threshold scatter that drives the effect.
    """
    def span(psw):
        v25, _ = level_crossing(inset_voltages, psw, 0.25)
        v50, _ = level_crossing(inset_voltages, psw, 0.50)
        return v25, v50, v50 - v25

    v25_n, v50_n, w_n = span(baseline_inset)
    v25_m, v50_m, w_m = span(mean_inset)
    sig_n = np.sqrt(0.25 / pooled_trials) / level_crossing(inset_voltages, baseline_inset)[1]
    sig_m = np.sqrt(0.25 / pooled_trials) / level_crossing(inset_voltages, mean_inset)[1]

    dev_cross = [level_crossing(inset_voltages, row)[0] for row in device_inset]
    finite = [c for c in dev_cross if np.isfinite(c)]

    def safe_fit(psw, y0_min):
        try:
            return logistic_slope(inset_voltages, psw, y0_min=y0_min)
        except Exception as exc:                      # degenerate pilot runs
            return dict(error=f"{type(exc).__name__}: {exc}")

    fits = {}
    for tag, y0_min in (("physical_y0ge0", 0.0), ("e2_convention_y0ge-0.05", -0.05)):
        fits[tag] = dict(nominal=safe_fit(baseline_inset, y0_min),
                         wafer_mean=safe_fit(mean_inset, y0_min))

    z = 1.96
    wilson_hw = lambda n: z * np.sqrt(0.25 / n) / (1 + z ** 2 / n)
    summary = dict(
        devices=args.devices, thermal_trials=args.thermal_trials,
        pooled_trials=pooled_trials, integrator=args.integrator, seed=args.seed,
        inset_voltages_mV=(np.asarray(inset_voltages) * 1e3).tolist(),
        psw_nominal=list(map(float, baseline_inset)),
        psw_wafer_mean=list(map(float, mean_inset)),
        crossing_nominal_mV=v50_n * 1e3, crossing_nominal_sigma_mV=float(sig_n * 1e3),
        crossing_wafer_mV=v50_m * 1e3, crossing_wafer_sigma_mV=float(sig_m * 1e3),
        span25to50_nominal_mV=w_n * 1e3, span25to50_wafer_mV=w_m * 1e3,
        span_ratio=float(w_m / w_n) if w_n else float("nan"),
        logistic_fits=fits,
        wilson_halfwidth_pooled=float(wilson_hw(pooled_trials)),
        device_crossings_mV=[float(c * 1e3) for c in dev_cross],
        device_crossing_mean_mV=float(np.mean(finite) * 1e3) if finite else None,
        device_crossing_std_mV=float(np.std(finite, ddof=1) * 1e3) if len(finite) > 1 else None,
        device_crossings_resolved=len(finite),
    )

    print()
    print("  D2D broadening summary (threshold-inset window):")
    print(f"    50% crossing      : nominal {summary['crossing_nominal_mV']:.1f} "
          f"+/- {summary['crossing_nominal_sigma_mV']:.1f} mV | wafer mean "
          f"{summary['crossing_wafer_mV']:.1f} +/- {summary['crossing_wafer_sigma_mV']:.1f} mV")
    print(f"    0.25->0.50 span   : nominal {summary['span25to50_nominal_mV']:.1f} mV | "
          f"wafer mean {summary['span25to50_wafer_mV']:.1f} mV "
          f"(x{summary['span_ratio']:.2f} broadening, fit-free)")
    for tag, pair in fits.items():
        if "error" in pair["nominal"] or "error" in pair["wafer_mean"]:
            print(f"    logistic [{tag}]: fit failed on this run")
            continue
        print(f"    logistic [{tag}]: nominal beta={pair['nominal']['beta']:.1f}"
              f"+/-{pair['nominal']['beta_ci95']:.1f} (R2={pair['nominal']['R2']:.4f}) | "
              f"wafer beta={pair['wafer_mean']['beta']:.1f}"
              f"+/-{pair['wafer_mean']['beta_ci95']:.1f} (R2={pair['wafer_mean']['R2']:.4f})")
    print(f"    Wilson 95% half-width at p=0.5, n={pooled_trials}: "
          f"+/-{summary['wilson_halfwidth_pooled']*100:.1f}%  "
          f"(nominal curve only; the wafer mean is a {args.devices}-device cluster sample)")
    if summary["device_crossing_std_mV"] is not None:
        print(f"    per-device 50% crossings: {summary['device_crossings_resolved']}"
              f"/{args.devices} resolved, mean {summary['device_crossing_mean_mV']:.0f} mV, "
              f"std {summary['device_crossing_std_mV']:.0f} mV")
    return summary


def main():
    args = parse_args()
    if args.devices < 1 or args.thermal_trials < 1:
        raise ValueError("--devices and --thermal-trials must be positive")
    t_start = time.time()

    base = PhysicalConstantsConfig()
    budget = cv_delta_budget(PDK_INPUTS)
    sample_rng = np.random.default_rng(args.seed)
    samples = make_process_samples(base, args.devices, sample_rng)
    self_heating = not args.disable_self_heating
    wide_voltages = nominal_voltage_sweep(WIDE_I_SOT_UA, base)
    inset_voltages = nominal_voltage_sweep(INSET_I_SOT_UA, base)
    voltages = np.concatenate((wide_voltages, inset_voltages))

    print("=" * 78)
    print("  Macrospin process-variability Monte Carlo")
    print("=" * 78)
    print(f"  PDK/Brinkman CV(Delta) reference = {budget.CV_Delta*100:.2f}%")
    print(f"  Wide terminal-voltage sweep    = {wide_voltages[0]*1e3:.0f}.."
          f"{wide_voltages[-1]*1e3:.0f} mV")
    print(f"  Threshold voltage inset        = {inset_voltages[0]*1e3:.0f}.."
          f"{inset_voltages[-1]*1e3:.0f} mV")
    print("  Process MC control variable    = fixed terminal |V_SOT|")
    print(f"  Self-heating                   = {self_heating}")
    print(f"  Integrator                     = {args.integrator}")
    print_scale_summary(samples)

    pooled_trials = args.devices * args.thermal_trials
    baseline_psw = run_psw_curve(
        base,
        voltages,
        pooled_trials,
        seed=args.seed + 17,
        self_heating=self_heating,
        integrator=args.integrator,
    )

    device_psw = []
    for idx, sample in enumerate(samples):
        curve = run_psw_curve(
            sample.constants,
            voltages,
            args.thermal_trials,
            seed=args.seed + 1009 * (idx + 1),
            self_heating=self_heating,
            integrator=args.integrator,
        )
        device_psw.append(curve)
    device_psw = np.asarray(device_psw, dtype=float)
    mean_psw = device_psw.mean(axis=0)
    baseline_wide, baseline_inset = split_sweeps(baseline_psw, len(wide_voltages))
    mean_wide, mean_inset = split_sweeps(mean_psw, len(wide_voltages))

    print()
    print("  Wide Psw table (nominal vs process-MC wafer mean):")
    for voltage, p_nom, p_mc in zip(wide_voltages, baseline_wide, mean_wide):
        print(f"    {voltage*1e3:6.0f} mV : nominal={p_nom:.3f}  process_MC={p_mc:.3f}")
    print()
    print("  Threshold-inset Psw table (nominal vs process-MC wafer mean):")
    for voltage, p_nom, p_mc in zip(inset_voltages, baseline_inset, mean_inset):
        print(f"    {voltage*1e3:6.0f} mV : nominal={p_nom:.3f}  process_MC={p_mc:.3f}")
    print("=" * 78)

    summary = summarize_broadening(inset_voltages, baseline_inset, mean_inset,
                                   device_psw[:, len(wide_voltages):],
                                   pooled_trials, args)
    summary["elapsed_s"] = time.time() - t_start
    print(f"  Wall-clock elapsed: {summary['elapsed_s']/3600:.2f} h")
    print("=" * 78)
    out_json = Path(__file__).resolve().parent / "macrospin_variability_summary.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Saved {out_json}")

    plot_curves(wide_voltages, baseline_wide, mean_wide,
                inset_voltages, baseline_inset, mean_inset, args, budget)


if __name__ == "__main__":
    main()
