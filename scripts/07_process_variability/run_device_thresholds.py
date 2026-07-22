"""Per-device threshold scatter of the 图2.17 mismatch ensemble.

The 图2.17 run averages 16 mismatched devices into one wafer-mean curve and
never exposes the individual devices, so two statements in §2.3.5 rest on
inference rather than measurement: that the wafer mean's 50% crossing moves
up because high-threshold samples delay saturation, and the long-standing
promise to look at device-to-device uniformity.

This script reconstructs the SAME 16 devices (identical sampler, seed and
antithetic pairing as macrospin_process_variability_mc.py) but scans each
one individually over a threshold window, so the D2D threshold distribution
is measured directly. It is deliberately cheap: only the transition window
is swept, not the full 22-point curve.

Note on scope: per-device *slope* is not extractable at this trial budget —
the 图2.17 ensemble spends 32 trials per device per point, and even the
512-trial pooled curve carries a +/-5 V^-1 slope interval. Only the
threshold (a level crossing) is well determined here.

Modes:
  python run_device_thresholds.py --device 0     # one shard -> dev_thr_d0.json
  python run_device_thresholds.py --analyze      # -> device_threshold_scatter.json
  python run_device_thresholds.py --smoke
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
from vgsot_sim.analysis.variability import sample_macrospin_process_constants

import importlib.util

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "_fig17", HERE / "macrospin_process_variability_mc.py")
_fig17 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fig17)

PULSE_NS, TOTAL_NS = _fig17.PULSE_NS, _fig17.TOTAL_NS
N_DEVICES, SEED, TRIALS = 16, 7, 48
PLATEAU_GATE = 0.5      # a device whose top-of-scan P_sw stays below this has
                        # no usable deterministic window; its crossing is noise
# Wider than the 图2.17 inset: individual devices scatter well outside the
# window that brackets the wafer mean.
SCAN_UA = np.array([700, 800, 900, 1000, 1080, 1160, 1240, 1330, 1450, 1600], float)


def rebuild_devices():
    """Reproduce the figure's 16 mismatch samples exactly (same sampler,
    same antithetic z-pairs, same master seed)."""
    base = PhysicalConstantsConfig()
    rng = np.random.default_rng(SEED)
    return _fig17.make_process_samples(base, N_DEVICES, rng)


def run_device(idx: int, trials: int = TRIALS):
    samples = rebuild_devices()
    sample = samples[idx]
    cc = sample.constants
    mid1 = int(round(PULSE_NS * 1e-9 / cc.t_step))
    end = int(round(TOTAL_NS * 1e-9 / cc.t_step))
    # Fixed terminal voltage is the control variable, exactly as in 图2.17:
    # each device sees the nominal voltage grid through its own R_W.
    v_grid = SCAN_UA * 1e-6 * PhysicalConstantsConfig().R_W
    cfg = SerSotNoVcmaThermalConfig(
        i_sot_list=tuple(-v_grid / cc.R_W), trials=trials,
        sim_start_step=1, sim_mid1_step=mid1, sim_end_step=end,
        pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
        target_mz=1.0, failure_tol=0.2, constants=cc,
    )
    t0 = time.time()
    res = ser_sot_no_vcma_thermal(
        cfg, show_progress=False, enable_self_heating=True, T_ambient_K=300.0,
        seed=SEED + 1009 * (idx + 1), rng_mode="generator", integrator="cayley",
    )
    return dict(
        device=idx, V_mV=(v_grid * 1e3).tolist(), psw=res.psw.tolist(),
        trials=trials, R_W=cc.R_W, Ki=cc.Ki, Ms=cc.Ms, tf=cc.tf,
        D_elec=cc.D_elec, scales={k: float(v) for k, v in sample.scales.items()},
        wall_s=time.time() - t0,
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )


def crossing(V_mV, P, trials, level=0.5):
    V = np.asarray(V_mV, float); P = np.asarray(P, float)
    o = np.argsort(V); V, P = V[o], P[o]
    for i in range(len(V) - 1):
        if (P[i] - level) * (P[i + 1] - level) <= 0 and P[i] != P[i + 1]:
            s = (P[i + 1] - P[i]) / (V[i + 1] - V[i])
            return float(V[i] + (level - P[i]) / s), float(np.sqrt(0.25 / trials) / abs(s))
    return float("nan"), float("nan")


def analyze():
    shards = []
    for f in sorted(glob.glob(str(HERE / "dev_thr_d*.json"))):
        shards.append(json.loads(Path(f).read_text(encoding="utf-8")))
    shards.sort(key=lambda d: d["device"])
    if not shards:
        raise SystemExit("No dev_thr_d*.json found.")

    V = np.asarray(shards[0]["V_mV"], float)
    rows, spans = [], []
    for d in shards:
        assert np.allclose(d["V_mV"], V), "shards must share the voltage grid"
        v50, sig = crossing(d["V_mV"], d["psw"], d["trials"])
        v25, _ = crossing(d["V_mV"], d["psw"], d["trials"], level=0.25)
        P = np.asarray(d["psw"], float)
        span = v50 - v25 if np.isfinite(v50) and np.isfinite(v25) else float("nan")
        spans.append(span)
        rows.append(dict(device=d["device"], vth_mV=v50, sigma_mV=sig,
                         span25to50_mV=span, plateau=float(P[-3:].max()),
                         R_W=d["R_W"], Ki_scale=d["scales"].get("Ki_proxy"),
                         Rsot_scale=d["scales"].get("Rsot")))

    def stats(vals):
        v = np.asarray([x for x in vals if np.isfinite(x)], float)
        if v.size == 0:
            return None
        return dict(n=int(v.size), mean_mV=float(v.mean()),
                    std_mV=float(v.std(ddof=1)) if v.size > 1 else None,
                    cv=float(v.std(ddof=1) / v.mean()) if v.size > 1 else None,
                    min_mV=float(v.min()), max_mV=float(v.max()))

    vth = np.array([r["vth_mV"] for r in rows], float)
    # A device whose curve never rises above 0.5 near the top of the scan has
    # no usable deterministic window; its nominal crossing is a noise artifact.
    # Both populations are reported so the exclusion's effect is visible.
    usable = np.array([r["plateau"] > PLATEAU_GATE for r in rows])

    # Direct average of the measured per-device curves. This is a model-free
    # cross-check of the 图2.17 wafer mean: same 16 devices, independent
    # voltage grid and trial count.
    P_all = np.array([d["psw"] for d in shards], float)
    mean_curve = P_all.mean(axis=0)
    mv50, _ = crossing(V, mean_curve, TRIALS * len(shards))
    mv25, _ = crossing(V, mean_curve, TRIALS * len(shards), level=0.25)

    summary = dict(
        devices=len(rows), trials_per_point=TRIALS,
        plateau_gate=PLATEAU_GATE,
        excluded_devices=[r["device"] for r, u in zip(rows, usable) if not u],
        scan_mV=V.tolist(), rows=rows,
        vth_all=stats(vth), vth_usable=stats(vth[usable]),
        span25to50_usable=stats(np.asarray(spans, float)[usable]),
        mean_curve_psw=mean_curve.tolist(),
        mean_curve_v25_mV=mv25, mean_curve_v50_mV=mv50,
        mean_curve_span_mV=(mv50 - mv25) if np.isfinite(mv50) and np.isfinite(mv25) else None,
        nominal_vth_mV=896.1,          # 图2.17 nominal raw crossing
        nominal_span_mV=51.5,          # 图2.17 nominal 0.25->0.50 span
        fig217_wafer_v50_mV=965.5, fig217_wafer_span_mV=198.2,
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    (HERE / "device_threshold_scatter.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print("=" * 72)
    print("Per-device 50% thresholds (0.75 ns, AP->P, Cayley, self-heating ON)")
    print("=" * 72)
    for r, u in zip(rows, usable):
        v, s = r["vth_mV"], r["span25to50_mV"]
        vs = f"{v:7.1f}" if np.isfinite(v) else "    n/a"
        ss = f"{s:5.1f}" if np.isfinite(s) else "  n/a"
        flag = "" if u else "   <- no usable plateau, crossing not meaningful"
        print(f"  device {r['device']:2d}: Vth = {vs} mV  span = {ss} mV  "
              f"(R_W = {r['R_W']:6.1f} ohm, plateau {r['plateau']:.2f}){flag}")
    for tag, st in (("all devices", summary["vth_all"]),
                    (f"plateau > {PLATEAU_GATE}", summary["vth_usable"])):
        if st and st["std_mV"]:
            print(f"\n  threshold [{tag}] n={st['n']}: mean {st['mean_mV']:.1f} mV, "
                  f"std {st['std_mV']:.1f} mV (CV {st['cv']*100:.1f}%), "
                  f"range {st['min_mV']:.0f}-{st['max_mV']:.0f} mV")
    sp = summary["span25to50_usable"]
    if sp and sp["std_mV"]:
        print(f"  per-device 0.25->0.50 span: mean {sp['mean_mV']:.1f} mV, "
              f"std {sp['std_mV']:.1f} mV, range {sp['min_mV']:.0f}-{sp['max_mV']:.0f} mV "
              f"(nominal single device: {summary['nominal_span_mV']:.1f} mV)")
    print(f"\n  direct average of the {len(rows)} measured curves: "
          f"V50 = {summary['mean_curve_v50_mV']:.1f} mV, "
          f"span = {summary['mean_curve_span_mV']:.1f} mV")
    print(f"  图2.17 wafer mean (independent grid, 32 trials): "
          f"V50 = {summary['fig217_wafer_v50_mV']:.1f} mV, "
          f"span = {summary['fig217_wafer_span_mV']:.1f} mV")
    print(f"  nominal single device: V50 = {summary['nominal_vth_mV']:.1f} mV, "
          f"span = {summary['nominal_span_mV']:.1f} mV")
    print("Wrote device_threshold_scatter.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--device", type=int, default=None)
    ap.add_argument("--trials", type=int, default=TRIALS)
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        samples = rebuild_devices()
        print(f"[smoke] rebuilt {len(samples)} devices; "
              f"R_W spread {min(s.constants.R_W for s in samples):.1f}-"
              f"{max(s.constants.R_W for s in samples):.1f} ohm")
        out = run_device(0, trials=4)
        print("[smoke] device 0 psw:", out["psw"], f"({out['wall_s']:.0f} s)")
        return
    if args.device is not None:
        print(f"[dev-thr] device {args.device}: {len(SCAN_UA)} pts x {args.trials} trials",
              flush=True)
        out = run_device(args.device, trials=args.trials)
        p = HERE / f"dev_thr_d{args.device}.json"
        p.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[dev-thr] wrote {p.name} (wall {out['wall_s']/60:.1f} min)", flush=True)
        return
    if args.analyze:
        analyze()
        return
    ap.error("choose one of --device / --analyze / --smoke")


if __name__ == "__main__":
    main()
