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


def crossing(V_mV, P, trials):
    V = np.asarray(V_mV, float); P = np.asarray(P, float)
    o = np.argsort(V); V, P = V[o], P[o]
    for i in range(len(V) - 1):
        if (P[i] - 0.5) * (P[i + 1] - 0.5) <= 0 and P[i] != P[i + 1]:
            s = (P[i + 1] - P[i]) / (V[i + 1] - V[i])
            return float(V[i] + (0.5 - P[i]) / s), float(np.sqrt(0.25 / trials) / abs(s))
    return float("nan"), float("nan")


def analyze():
    shards = []
    for f in sorted(glob.glob(str(HERE / "dev_thr_d*.json"))):
        shards.append(json.loads(Path(f).read_text(encoding="utf-8")))
    shards.sort(key=lambda d: d["device"])
    if not shards:
        raise SystemExit("No dev_thr_d*.json found.")

    rows = []
    for d in shards:
        v50, sig = crossing(d["V_mV"], d["psw"], d["trials"])
        P = np.asarray(d["psw"], float)
        rows.append(dict(device=d["device"], vth_mV=v50, sigma_mV=sig,
                         plateau=float(P[-3:].max()), R_W=d["R_W"],
                         Ki_scale=d["scales"].get("Ki_proxy"),
                         Rsot_scale=d["scales"].get("Rsot")))
    vth = np.array([r["vth_mV"] for r in rows], float)
    ok = np.isfinite(vth)
    summary = dict(
        devices=len(rows), resolved=int(ok.sum()), trials_per_point=TRIALS,
        scan_mV=[float(v) for v in np.asarray(SCAN_UA) * 1e-6 * PhysicalConstantsConfig().R_W * 1e3],
        rows=rows,
        vth_mean_mV=float(np.mean(vth[ok])) if ok.any() else None,
        vth_std_mV=float(np.std(vth[ok], ddof=1)) if ok.sum() > 1 else None,
        vth_min_mV=float(np.min(vth[ok])) if ok.any() else None,
        vth_max_mV=float(np.max(vth[ok])) if ok.any() else None,
        nominal_vth_mV=896.1,          # 图2.17 nominal raw crossing
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    if summary["vth_std_mV"]:
        summary["vth_cv"] = summary["vth_std_mV"] / summary["vth_mean_mV"]
    (HERE / "device_threshold_scatter.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print("=" * 66)
    print("Per-device 50% thresholds (0.75 ns, P->AP, Cayley, self-heating ON)")
    print("=" * 66)
    for r in rows:
        v = r["vth_mV"]
        vs = f"{v:7.1f}" if np.isfinite(v) else "    n/a"
        print(f"  device {r['device']:2d}: Vth = {vs} mV  "
              f"(R_W = {r['R_W']:6.1f} ohm, plateau {r['plateau']:.2f})")
    if summary["vth_std_mV"]:
        print(f"\n  resolved {summary['resolved']}/{summary['devices']}: "
              f"mean {summary['vth_mean_mV']:.1f} mV, std {summary['vth_std_mV']:.1f} mV "
              f"(CV {summary['vth_cv']*100:.1f}%), range "
              f"{summary['vth_min_mV']:.0f}-{summary['vth_max_mV']:.0f} mV")
        print(f"  nominal reference: {summary['nominal_vth_mV']:.0f} mV")
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
