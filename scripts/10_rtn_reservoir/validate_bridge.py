"""Stage-2 sLLG↔RTN bridge validation (background-runnable).

Drives the repo's macrospin LLG engine in a low-barrier free-running regime and
checks the three RTN signatures + calibrates tau0 and the bias→V map. Writes a
JSON + human summary next to this script. See this directory's README.

Run:  PYTHONPATH=src python scripts/10_rtn_reservoir/validate_bridge.py
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.rtn import bridge

OUT = os.path.join(os.path.dirname(__file__), "bridge_results.json")
DT = 4e-12          # 4 ps: validated t_step-invariant vs 1 ps (10 ps drifts)
c0 = PhysicalConstantsConfig()
c0 = type(c0)(**{**c0.__dict__, "t_step": DT})


def section(title):
    print(f"\n=== {title} ===", flush=True)


results = {"t_step_s": DT, "tilt_per_field_per_Am": bridge.tilt_per_field(c0)}
t_start = time.time()

# ---- 1. main run at Delta=2, zero bias: dwell + PSD + tau0 + ergodic <mz> ----
section("Delta=2.0 main run (dwell / PSD / tau0)")
DELTA = 2.0
c = bridge.low_barrier_constants(DELTA, c0)
n_main = 2_000_000                       # 8 us
t0 = time.time()
mz = bridge.free_run(c, n_main, h_ex_z=0.0, seed=101)
dw = bridge.dwell_times(mz, DT, thr=0.5)
psd = bridge.psd_lorentzian(mz, DT, n_seg=8)
tau0 = bridge.infer_tau0(dw.tau_mean_ns, DELTA)
print(f"  {n_main} steps ({n_main*DT*1e9:.0f} ns) in {time.time()-t0:.0f}s", flush=True)
print(f"  <mz>={mz.mean():+.3f}  frac|mz|>0.5={np.mean(np.abs(mz)>0.5):.2f}", flush=True)
print(f"  flips={dw.n_flips}  tau_dwell={dw.tau_mean_ns:.1f} ns  CV={dw.cv:.2f}", flush=True)
print(f"  tau0_inferred={tau0:.1f} ns (RTN model assumes 1 ns)", flush=True)
print(f"  PSD corner f_c={psd.f_corner_hz:.3e} Hz -> tau_c={psd.tau_corner_ns:.1f} ns "
      f"(expect ~tau_dwell/2={dw.tau_mean_ns/2:.1f} ns)", flush=True)
results["main_delta2"] = {
    "delta": DELTA, "n_steps": n_main, "mz_mean": float(mz.mean()),
    "frac_in_state": float(np.mean(np.abs(mz) > 0.5)),
    "n_flips": dw.n_flips, "tau_dwell_ns": dw.tau_mean_ns, "dwell_cv": dw.cv,
    "tau0_inferred_ns": tau0, "psd_fc_hz": psd.f_corner_hz, "psd_tau_c_ns": psd.tau_corner_ns,
    "dwell_ns": dw.dwell_ns.tolist(),
}

# ---- 2. tau0 vs Delta (is the attempt time ~Delta-independent?) ----
section("tau0 vs Delta")
scan = []
for d in (1.5, 2.0, 2.5):
    cd = bridge.low_barrier_constants(d, c0)
    n = 1_200_000
    t0 = time.time()
    mzd = bridge.free_run(cd, n, h_ex_z=0.0, seed=202)
    dwd = bridge.dwell_times(mzd, DT, thr=0.5)
    t0i = bridge.infer_tau0(dwd.tau_mean_ns, d)
    print(f"  Delta={d}: flips={dwd.n_flips} tau_dwell={dwd.tau_mean_ns:.1f}ns "
          f"CV={dwd.cv:.2f} tau0={t0i:.1f}ns  ({time.time()-t0:.0f}s)", flush=True)
    scan.append({"delta": d, "n_flips": dwd.n_flips, "tau_dwell_ns": dwd.tau_mean_ns,
                 "cv": dwd.cv, "tau0_inferred_ns": t0i})
results["tau0_vs_delta"] = scan

# ---- 3. <m_z> vs longitudinal bias -> tanh calibration (Delta=2) ----
section("<m_z> vs bias (tanh calibration)")
hz_list = (-3000.0, -1500.0, 0.0, 1500.0, 3000.0)
hz, tilt, mzmean, tanh_pred = bridge.mean_mz_vs_bias(
    c, hz_list, 1_200_000, seeds=(301, 302))
for h, ti, mm, tp in zip(hz, tilt, mzmean, tanh_pred):
    print(f"  hz={h:+7.0f} A/m  tilt={ti:+5.2f}  <mz>={mm:+.3f}  tanh={tp:+.3f}", flush=True)
# linear-region slope check: d<mz>/d(tilt) near 0 should be ~1 for a clean tanh
results["tanh_calib"] = {"hz": hz.tolist(), "tilt": tilt.tolist(),
                         "mz_mean": mzmean.tolist(), "tanh_pred": tanh_pred.tolist()}

results["elapsed_s"] = time.time() - t_start
with open(OUT, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nWrote {OUT}  (total {results['elapsed_s']:.0f}s)", flush=True)
