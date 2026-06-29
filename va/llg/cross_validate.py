#!/usr/bin/env python3
"""Cross-engine regression: Python `switching_vector` vs `vgsot_llg.va` (ngspice).

Drives the SAME deterministic SOT setup through both engines and asserts their
m_z(t) trajectories agree. Two engines, one physics. Graceful-skip (exit 0) if
OpenVAF / ngspice are not installed (like eda/testbenches/run_regression.py).

Tool discovery: env OPENVAF / NGSPICE, then PATH.
Run:  python va/llg/cross_validate.py
"""
from __future__ import annotations
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]                 # va/llg -> va -> <vgsot-sim repo root>
sys.path.insert(0, str(REPO / "src"))
VA = HERE / "vgsot_llg.va"
OSDI = HERE / "vgsot_llg.osdi"
TB = HERE / "_xval_tb.spice"
OUT = HERE / "_xval_out.csv"
SUMMARY = HERE / "cross_validate_summary.json"

# --- shared deterministic setup (super-threshold SOT, V_MTJ=0, start near +z) ---
ISOT = -3.0e-3        # A
VMTJ = 0.0            # V  (no STT)
M0 = np.array([0.141, 0.0, 0.99]); M0 = M0 / np.linalg.norm(M0)
DUR_NS = 3.0
RSOT = 776.0
TOL = 0.02            # max |Δm_z| over the window


def find_tool(env_key, names):
    v = os.environ.get(env_key)
    if v and Path(v).exists():
        return v
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None


def python_ref():
    from vgsot_sim.configs import PhysicalConstantsConfig
    from vgsot_sim.dynamic_switching_vector import switching_vector
    from vgsot_sim.initialize import compute_Rp
    cc = PhysicalConstantsConfig()
    n = int(round(DUR_NS * 1e-9 / cc.t_step))
    m = M0.copy()
    Rp = compute_Rp(cc)
    mz = np.empty(n); t = np.empty(n)
    for i in range(n):
        m = switching_vector(m, VMTJ, ISOT, Rp, 0, 1, VNV=0, NON=0,
                             R_SOT_FL_DL=0.83, constants=cc)
        mz[i] = m[2]; t[i] = (i + 1) * cc.t_step
    return t, mz


def write_tb():
    TB.write_text(
        "* cross-engine deterministic cross-validation testbench\n"
        ".model vgsot_llg vgsot_llg\n"
        "N1 sot_p sot_n mtj_p mx my mz hx hy hz vgsot_llg\n"
        f"Vsot sot_p 0 dc {ISOT * RSOT}\n"
        "Vsn sot_n 0 dc 0\n"
        f"Vmtj mtj_p 0 dc {VMTJ}\n"
        "Vhx hx 0 dc 0\nVhy hy 0 dc 0\nVhz hz 0 dc 0\n"
        f".ic v(mx)={M0[0]} v(my)={M0[1]} v(mz)={M0[2]}\n"
        ".control\n"
        f"  tran 1p {DUR_NS}n uic\n"
        f"  wrdata {OUT.name} v(mx) v(my) v(mz)\n"
        "  quit\n.endc\n.end\n"
    )


def main():
    openvaf = find_tool("OPENVAF", ["openvaf", "openvaf-r"])
    ngspice = find_tool("NGSPICE", ["ngspice_con", "ngspice"])
    if not openvaf or not ngspice:
        print("OpenVAF/ngspice not found; cross-engine check skipped (see va/README.md).")
        return 0

    subprocess.run([openvaf, str(VA), "-o", str(OSDI)], check=True, cwd=HERE)
    (HERE / ".spiceinit").write_text(f"osdi {OSDI.name}\n")
    write_tb()
    subprocess.run([ngspice, "-b", TB.name], check=True, cwd=HERE)

    r = np.loadtxt(OUT)
    t_va, mz_va = r[:, 0], r[:, 5]
    t_py, mz_py = python_ref()
    mz_py_i = np.interp(t_va, t_py, mz_py)
    err = np.abs(mz_py_i - mz_va)
    maxerr = float(err.max())
    eq_err = float(abs(mz_py[-1] - mz_va[-1]))
    ok = maxerr < TOL

    SUMMARY.write_text(json.dumps({
        "kind": "cross_engine_llg_va_vs_python",
        "note": "deterministic m_z(t) agreement; two engines, one physics.",
        "I_SOT_A": ISOT, "V_MTJ": VMTJ, "dur_ns": DUR_NS,
        "max_abs_dmz": maxerr, "equilibrium_abs_dmz": eq_err,
        "tol": TOL, "pass": ok, "n_va_points": int(len(t_va)),
        "openvaf": str(openvaf), "ngspice": str(ngspice),
    }, indent=2))
    print(f"cross-engine Python vs vgsot_llg.va: max|dmz|={maxerr:.4f} "
          f"eq|dmz|={eq_err:.4f} (tol {TOL}) -> {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
