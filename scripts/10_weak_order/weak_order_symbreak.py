"""Symmetry-breaking robustness check for the E-weak p_weak~2 result.

The committed weak_order_audit.py measures weak order at an AXIALLY SYMMETRIC
operating point (pure PMA+demag well, isotropic noise, observable m_z is the
symmetry invariant). If the observed p_weak~2.0 is a symmetry artifact (leading
odd-order weak-error terms cancelling by rotational symmetry about z), it will
NOT survive a transverse field that breaks the symmetry. If it is genuine, it
will. Same vectorised stepper as the committed script (verified == scalar API
to 2.5e-16), plus a constant transverse field H_x.

Operating point: Ms_scale=1e-3 (same), H_x = 5e7 A/m ~ 7% of the scaled
H_PMA(2Ki/(u0 Ms tf) ~ 7.4e8 A/m) -> solidly symmetry-broken, still a
perturbation of the well.
"""
import json
import math
import time

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.demag import demag_factors

MS_SCALE = 1e-3
H_X = 5.0e7                      # transverse symmetry-breaking field [A/m]
M0 = np.array([math.sqrt(1 - 0.5 ** 2), 0.0, 0.5])
T_FIN = 2.0e-13
N_REF = 1024
FIT_NSTEPS = (256, 128, 64, 32, 16)
M_TRAJ = 150000

cc = PhysicalConstantsConfig()
cc.Ms = PhysicalConstantsConfig().Ms * MS_SCALE
gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
alpha = cc.alpha
Ki, Ms, u0, tf = cc.Ki, cc.Ms, cc.u0, cc.tf
Nx, Ny, Nz = demag_factors(cc, mode="ellipsoid")
A0 = math.sqrt(2.0 * cc.kb * cc.T * cc.alpha / (cc.u0 * cc.Ms * cc.gamma * cc.v))
H_PMA_mag = 2 * Ki / (u0 * Ms * tf)
print(f"H_PMA={H_PMA_mag:.3e} A/m, H_x={H_X:.3e} ({H_X/H_PMA_mag:.1%} of well)", flush=True)


def run_ensemble(dt, N, M, stepper, seed):
    rng = np.random.default_rng(seed)
    sqdt = math.sqrt(dt)
    m = np.tile(M0, (M, 1))

    def omega_vec(mv, h_th):
        Hx = -Ms * Nx * mv[:, 0] + H_X
        Hy = -Ms * Ny * mv[:, 1]
        Hz = -Ms * Nz * mv[:, 2] + H_PMA_mag * mv[:, 2]
        H = np.stack([Hx, Hy, Hz], axis=1) + h_th
        return gamma_red * H + alpha * gamma_red * np.cross(mv, H)

    def cayley_raw(mv, om):
        s = 0.5 * dt
        c1 = np.cross(om, mv)
        c2 = np.cross(om, c1)
        osq = np.sum(om * om, axis=1)
        fac = (2.0 * s) / (1.0 + s * s * osq)
        return mv + fac[:, None] * (c1 + s * c2)

    for _ in range(N):
        dW = rng.normal(0.0, sqdt, (M, 3))
        h_th = A0 * dW / dt
        if stepper == "published":
            om = omega_vec(m, h_th)
            mn = cayley_raw(m, om)
        else:  # midpoint
            mn = m.copy()
            for _ in range(3):
                mmid = 0.5 * (m + mn)
                mmid = mmid / np.linalg.norm(mmid, axis=1)[:, None]
                om = omega_vec(mmid, h_th)
                mn = cayley_raw(m, om)
        m = mn / np.linalg.norm(mn, axis=1)[:, None]
    mz = m[:, 2]
    return float(mz.mean()), float((mz * mz).mean()), float(mz.std() / math.sqrt(M))


def fit_order(dts, errs):
    x = np.log(np.asarray(dts, float)); y = np.log(np.asarray(errs, float))
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    r2 = 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)
    return float(coef[0]), float(r2)


t0 = time.time()
dt_ref = T_FIN / N_REF
Emz_A, Emz2_A, se_A = run_ensemble(dt_ref, N_REF, M_TRAJ * 3, "published", 91000)
Emz_B, Emz2_B, se_B = run_ensemble(dt_ref, N_REF, M_TRAJ * 3, "midpoint", 92000)
print(f"ref A(self) E[mz]={Emz_A:.6f} (se {se_A:.1e}); ref B(mid) E[mz]={Emz_B:.6f}; |A-B|={abs(Emz_A-Emz_B):.2e}", flush=True)

dts, errA, errB, errA2, errB2 = [], [], [], [], []
for n in FIT_NSTEPS:
    dt = T_FIN / n
    mz, mz2, se = run_ensemble(dt, n, M_TRAJ, "published", 30000 + n)
    dts.append(dt)
    errA.append(abs(mz - Emz_A)); errB.append(abs(mz - Emz_B))
    errA2.append(abs(mz2 - Emz2_A)); errB2.append(abs(mz2 - Emz2_B))
    print(f"  n={n:4d} dt={dt:.3e} E[mz]={mz:.6f}(se {se:.1e}) errA={errA[-1]:.2e} errB={errB[-1]:.2e}", flush=True)

res = {"H_x_over_well": H_X / H_PMA_mag}
for name, e in [("mz_selfref", errA), ("mz_midref", errB),
                ("mz2_selfref", errA2), ("mz2_midref", errB2)]:
    p, r2 = fit_order(dts, e)
    p3, _ = fit_order(dts[:3], e[:3])
    res[name] = dict(p_full=p, r2=r2, p_finest3=p3)
    print(f"  SYMBREAK WEAK fit {name}: p_full={p:.3f} (R2={r2:.3f}) p_finest3={p3:.3f}", flush=True)

pooled = np.median([res[k]["p_full"] for k in ("mz_selfref", "mz_midref", "mz2_selfref", "mz2_midref")])
print(f"\nSYMBREAK HEADLINE pooled p_weak = {pooled:.3f}  (runtime {time.time()-t0:.0f}s)", flush=True)
print("p~2 SURVIVES symmetry breaking" if pooled > 1.7 else
      ("p DEGRADES toward 1 -> original 2.0 was a symmetry artifact" if pooled < 1.3 else
       "INTERMEDIATE -> needs finer look"), flush=True)
with open("result/sec_2_2_3_2/E/weak_order_symbreak.json", "w") as f:
    json.dump(res, f, indent=2)
