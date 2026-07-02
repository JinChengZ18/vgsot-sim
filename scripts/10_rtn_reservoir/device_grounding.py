"""Stage-3.5 device grounding: translate the RC operating point into energy / latency / window.

Anchored to the PhysicalConstantsConfig device numbers and the §2.3 write energy
(0.78 pJ). Run: PYTHONPATH=src python scripts/10_rtn_reservoir/device_grounding.py

Caveat: the RTN 'V' is a phenomenological tilt whose mapping to a real terminal
drive is unresolved (Stage 2), so the per-node energy is an order-of-magnitude
estimate (two bounding electrical paths are shown), not a calibrated figure.
"""
from __future__ import annotations

import json
import os

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.rtn import telegraph as tg

OUT = os.path.join(os.path.dirname(__file__), "device_grounding.json")
c = PhysicalConstantsConfig()
R_SOT = c.R_W
R_MTJ = c.RA / c.A1
E_WRITE_PJ = 0.78

res = {"R_SOT_ohm": R_SOT, "R_MTJ_ohm": R_MTJ, "E_write_pJ": E_WRITE_PJ, "points": []}
print(f"R_SOT = {R_SOT:.0f} ohm   R_MTJ = {R_MTJ/1e3:.1f} kohm   E_write(sec 2.3) = {E_WRITE_PJ} pJ\n")

for tau0_label, tau0 in [("model 1 ns", 1.0), ("LLG-measured ~20 ns", 20.0)]:
    print(f"--- attempt time tau0 = {tau0_label} (sets dt) ---")
    for V in (0.1, 0.2, 0.3):
        dt = tau0 * 1e-9
        E_mtj = V ** 2 / R_MTJ * dt                  # bias across the MTJ (read path)
        E_sot = V ** 2 / R_SOT * dt                  # bias across the SOT channel
        tau_mem = float(tg.relaxation_time(V, tau0=tau0, Delta=2.0, Vc0=0.884))
        thr = 1.0 / dt / 1e6                          # MOPS per node
        print(f"  V={V} V, dt={tau0:.0f} ns: E_MTJ={E_mtj*1e15:6.1f} fJ  E_SOT={E_sot*1e12:4.2f} pJ  "
              f"tau_mem(Delta=2)={tau_mem:.0f} ns  throughput={thr:.0f} MOPS/node")
        res["points"].append({"tau0_ns": tau0, "V": V, "E_mtj_fJ": E_mtj * 1e15,
                              "E_sot_pJ": E_sot * 1e12, "tau_mem_ns": tau_mem, "mops_per_node": thr})
    print()

print("usable Delta window for RC (tau0=1 ns units):")
window = []
for D in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0):
    tmax = tg.tau_max(Delta=D)
    slope = D / 0.884
    v_sat = 0.884 * 2.0 / D
    print(f"  Delta={D}: tau_max={tmax:.1f} ns  d<s>/dV|0={slope:.2f}/V  tanh saturates ~|V|={v_sat:.2f} V")
    window.append({"Delta": D, "tau_max_ns": tmax, "slope_per_V": slope, "v_sat_V": v_sat})
res["delta_window"] = window

with open(OUT, "w") as f:
    json.dump(res, f, indent=2)
print(f"\nWrote {OUT}")
