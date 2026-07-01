"""Exp E (weak order) — weak stochastic convergence order of the published
explicit-omega Cayley sLLG step.

Thesis section: §2.2.3.2, Experiment E.

WHAT THIS MEASURES
------------------
§2.2.3.2 inherits (by citation [31], implicit-midpoint property) WEAK order 2.0.
The published code uses a LEFT-endpoint explicit-omega Cayley step.  Hypothesis
H1: with no midpoint symmetrisation and no Levy-area term, the WEAK order is ~1,
not 2.

METHOD
------
Fixed short horizon T_fin (pre-equilibrium) from a FIXED initial state m0.
For functionals g in {m_z, m_z^2} estimate
    weak error  e_g(dt) = | E[g(m_dt(T_fin))] - E[g(m_ref(T_fin))] |
where E[.] is a Monte-Carlo average over M INDEPENDENT noise realisations (fresh
per-step Gaussian increments, the engine's native noise model — weak order does
NOT require a common Brownian path).  The fine reference uses dt_ref << every
fitted dt.  Order p_weak = slope of log e vs log dt.

TWO independent reference constructions (red-team cross-check):
  (A) self-Cayley: reference = published explicit-Cayley at the fine grid.
  (B) Stratonovich-Heun midpoint: reference = an INDEPENDENT integrator (true
      implicit-midpoint Cayley, 3 fixed-point iterations) at the fine grid.
A weak order that agrees between the two references is trustworthy.

NOISE INJECTION
---------------
Per step we inject  h_th_ext = A0 * dW / dt,  dW ~ N(0, dt) drawn fresh each
step (independent across steps and trajectories), with the dt-independent FDT
amplitude  A0 = sqrt(2 kB T alpha /(u0 Ms gamma v)).  This reproduces the
engine's internal field A0/sqrt(dt)*N(0,1) exactly (FDT convention) while letting
the harness own the RNG for reproducibility and for using the SAME amplitude
across all dt (variance per unit time conserved).

DECISIVE
--------
p_weak ~ 1 (both references agree, upper CI < ~1.5) -> the cited weak order 2.0
is falsified.

Run:  PYTHONPATH=src python scripts/10_weak_order/weak_order_audit.py
"""
from __future__ import annotations

import json
import math
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.anisotropy import field
from vgsot_sim.dynamic_switching_vector import _omega_from_state, cayley_step

OUTDIR = os.path.join("result", "sec_2_2_3_2", "E")
os.makedirs(OUTDIR, exist_ok=True)

# Pre-equilibrium operating point.  m0 is tilted WELL off the easy axis
# (m_z0 = 0.5) so E[m_z], E[m_z^2] move substantially and SMOOTHLY over the
# horizon; pure multiplicative noise (ENE=0, no drive).  At Ms_scale=1e-3 the
# PMA restoring field is large, so the relaxation is fast: T_fin = 2e-13 keeps
# the run genuinely pre-equilibrium (E[m_z] moves 0.5 -> ~0.94, not yet at the
# ~0.998 steady state).  A horizon long enough to equilibrate would turn the
# weak-error test into a stationary-distribution (Boltzmann) test instead.
M0 = np.array([np.sqrt(1.0 - 0.5 ** 2), 0.0, 0.5])
T_FIN = 2.0e-13
MS_SCALE = 1e-3                 # noise-dominated (same regime as strong-order Exp)
N_REF = 1024                    # reference grid steps over T_FIN (dt_ref << every fitted dt)
FIT_NSTEPS = (256, 128, 64, 32, 16)   # fitted coarse grids (all coarser than ref)
M_TRAJ = 150000                # MC trajectories per (dt, scheme)


def fit_order(dts, errs):
    x = np.log(np.asarray(dts, float))
    y = np.log(np.asarray(errs, float))
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(coef[0]), float(coef[1]), float(r2)


def make_cc(Ms_scale=MS_SCALE):
    cc = PhysicalConstantsConfig()
    cc.Ms = PhysicalConstantsConfig().Ms * Ms_scale
    return cc


def fdt_amplitude(cc, T=None):
    T_use = cc.T if T is None else T
    return math.sqrt(2.0 * cc.kb * T_use * cc.alpha
                     / (cc.u0 * cc.Ms * cc.gamma * cc.v))


def omega_at(m, cc, A0, h_th):
    m = np.asarray(m, float)
    theta = float(np.arccos(np.clip(m[2], -1.0, 1.0)))
    phi = float(np.arctan2(m[1], m[0]))
    H_eff, _ = field(theta, phi, 0.0, n=1, NON=1, ENE=0, VNV=0,
                     constants=cc, demag_mode="ellipsoid", h_th_ext=h_th, T=None)
    H_eff = np.asarray(H_eff, float)
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
    return _omega_from_state(m, H_eff, np.array([-1.0, 0.0, 0.0]),
                             H_DL_SOT=0.0, H_FL_SOT=0.0, H_DL_STT=0.0,
                             H_FL_STT=0.0, sigma_STT=np.array([0.0, 0.0, 1.0]),
                             alpha=cc.alpha, gamma_red=gamma_red)


def step_published(m, cc, dt, dW, A0):
    h = A0 * dW / dt
    om = omega_at(m, cc, A0, h)
    mn = cayley_step(m, om, dt)
    return mn / np.linalg.norm(mn)


def step_midpoint(m, cc, dt, dW, A0, iters=3):
    h = A0 * dW / dt
    mn = m.copy()
    for _ in range(iters):
        mmid = (m + mn) / 2.0
        mmid = mmid / np.linalg.norm(mmid)
        om = omega_at(mmid, cc, A0, h)
        mn = cayley_step(m, om, dt)
    return mn / np.linalg.norm(mn)


def mc_expectations(cc, dt, N, A0, M, stepper, seed):
    """Monte-Carlo E[m_z], E[m_z^2] at the terminal time over M trajectories.

    Vectorised over trajectories: m is (M,3); per-step fresh Gaussian (M,3).
    The published `_omega_from_state`/`field` path is scalar, so we vectorise a
    minimal copy of the pure-noise omega here for speed (identical math)."""
    rng = np.random.default_rng(seed)
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
    alpha = cc.alpha
    # constant (state-independent up to m) effective-field pieces for ENE=0,
    # VNV=0, pure noise:  H_eff = H_PMA(m_z) + H_D(m) + H_TH.
    Ki, Ms, u0, tf = cc.Ki, cc.Ms, cc.u0, cc.tf
    from vgsot_sim.demag import demag_factors
    Nx, Ny, Nz = demag_factors(cc, mode="ellipsoid")
    sqdt = math.sqrt(dt)

    m = np.tile(np.asarray(M0, float), (M, 1))

    def omega_vec(mv, h_th):
        # H_PMA along z
        Hx = -Ms * Nx * mv[:, 0]
        Hy = -Ms * Ny * mv[:, 1]
        Hz = -Ms * Nz * mv[:, 2] + (2 * Ki / (u0 * Ms * tf)) * mv[:, 2]
        H = np.stack([Hx, Hy, Hz], axis=1) + h_th
        mxH = np.cross(mv, H)
        return gamma_red * H + alpha * gamma_red * mxH

    def cayley_vec(mv, om):
        s = 0.5 * dt
        c1 = np.cross(om, mv)
        c2 = np.cross(om, c1)
        osq = np.sum(om * om, axis=1)
        fac = (2.0 * s) / (1.0 + s * s * osq)
        mn = mv + fac[:, None] * (c1 + s * c2)
        return mn / np.linalg.norm(mn, axis=1)[:, None]

    def cayley_mid(mv, h_th):
        mn = mv.copy()
        for _ in range(3):
            mmid = 0.5 * (mv + mn)
            mmid = mmid / np.linalg.norm(mmid, axis=1)[:, None]
            om = omega_vec(mmid, h_th)
            mn = cayley_vec_raw(mv, om)
        return mn / np.linalg.norm(mn, axis=1)[:, None]

    def cayley_vec_raw(mv, om):
        s = 0.5 * dt
        c1 = np.cross(om, mv)
        c2 = np.cross(om, c1)
        osq = np.sum(om * om, axis=1)
        fac = (2.0 * s) / (1.0 + s * s * osq)
        return mv + fac[:, None] * (c1 + s * c2)

    use_mid = (stepper == "midpoint")
    for _ in range(N):
        dW = rng.normal(0.0, sqdt, (M, 3))
        h_th = A0 * dW / dt
        if use_mid:
            m = cayley_mid(m, h_th)
        else:
            om = omega_vec(m, h_th)
            m = cayley_vec(m, om)
    mz = m[:, 2]
    return mz, (mz * mz)


def main():
    t0 = time.time()
    res = {"experiment": "E_weak_order", "section": "2.2.3.2",
           "T_fin_s": T_FIN, "N_ref": N_REF, "fit_nsteps": list(FIT_NSTEPS),
           "M_traj": M_TRAJ, "Ms_scale": MS_SCALE, "m0": M0.tolist()}
    cc = make_cc()
    A0 = fdt_amplitude(cc)
    res["A0"] = A0

    # reference grids
    dt_ref = T_FIN / N_REF
    # (A) self-Cayley reference and (B) midpoint reference, each with large M.
    M_ref = M_TRAJ * 3
    print(f"[REF] M_ref={M_ref}, N_ref={N_REF}, dt_ref={dt_ref:.3e}", flush=True)
    mz_A, mz2_A = mc_expectations(cc, dt_ref, N_REF, A0, M_ref, "published", seed=10001)
    Emz_A, Emz2_A = float(mz_A.mean()), float(mz2_A.mean())
    mz_B, mz2_B = mc_expectations(cc, dt_ref, N_REF, A0, M_ref, "midpoint", seed=20002)
    Emz_B, Emz2_B = float(mz_B.mean()), float(mz2_B.mean())
    # MC std-error on the reference means (for honest CI on weak error floor)
    se_Emz_A = float(mz_A.std() / math.sqrt(M_ref))
    se_Emz_B = float(mz_B.std() / math.sqrt(M_ref))
    print(f"  ref A (self-Cayley): E[mz]={Emz_A:.6f} E[mz^2]={Emz2_A:.6f} (se_mz={se_Emz_A:.2e})", flush=True)
    print(f"  ref B (midpoint)   : E[mz]={Emz_B:.6f} E[mz^2]={Emz2_B:.6f} (se_mz={se_Emz_B:.2e})", flush=True)
    res["reference"] = dict(self_cayley=dict(Emz=Emz_A, Emz2=Emz2_A, se_Emz=se_Emz_A),
                            midpoint=dict(Emz=Emz_B, Emz2=Emz2_B, se_Emz=se_Emz_B),
                            consistency_mz=abs(Emz_A - Emz_B),
                            consistency_mz2=abs(Emz2_A - Emz2_B))
    print(f"  reference consistency |A-B|: mz={abs(Emz_A-Emz_B):.2e} mz2={abs(Emz2_A-Emz2_B):.2e}", flush=True)

    dts = [T_FIN / n for n in FIT_NSTEPS]
    err_mz_A, err_mz2_A, err_mz_B, err_mz2_B = [], [], [], []
    se_mz, se_mz2 = [], []
    for n in FIT_NSTEPS:
        dt = T_FIN / n
        mz, mz2 = mc_expectations(cc, dt, n, A0, M_TRAJ, "published", seed=30000 + n)
        Emz, Emz2 = float(mz.mean()), float(mz2.mean())
        s_mz = float(mz.std() / math.sqrt(M_TRAJ))
        s_mz2 = float(mz2.std() / math.sqrt(M_TRAJ))
        err_mz_A.append(abs(Emz - Emz_A)); err_mz2_A.append(abs(Emz2 - Emz2_A))
        err_mz_B.append(abs(Emz - Emz_B)); err_mz2_B.append(abs(Emz2 - Emz2_B))
        se_mz.append(s_mz); se_mz2.append(s_mz2)
        print(f"  n={n:4d} dt={dt:.3e} E[mz]={Emz:.6f}(se {s_mz:.1e}) "
              f"err_mz(A)={abs(Emz-Emz_A):.2e} err_mz(B)={abs(Emz-Emz_B):.2e}", flush=True)

    fits = {}
    for name, errs in [("mz_selfref", err_mz_A), ("mz2_selfref", err_mz2_A),
                       ("mz_midref", err_mz_B), ("mz2_midref", err_mz2_B)]:
        p, b, r2 = fit_order(dts, errs)
        # asymptotic sub-window slope (finest 3 dt) — least contaminated by the
        # large-error (coarse-dt) nonlinear regime.
        p_fine, _, _ = fit_order(dts[:3], errs[:3])
        fits[name] = dict(p=p, r2=r2, p_finest3=p_fine, errs=errs)
        print(f"  WEAK fit {name}: p_full={p:.3f} (R2={r2:.3f}) p_finest3={p_fine:.3f}", flush=True)
    res["dts"] = dts
    res["se_mz"] = se_mz
    res["se_mz2"] = se_mz2
    res["fits"] = fits

    # headline = m_z weak order, both references
    p_mz_A = fits["mz_selfref"]["p"]
    p_mz_B = fits["mz_midref"]["p"]
    p_mz2_A = fits["mz2_selfref"]["p"]
    p_mz2_B = fits["mz2_midref"]["p"]
    # asymptotic (finest-window) slopes, the least-contaminated estimate
    pf_mz_A = fits["mz_selfref"]["p_finest3"]
    pf_mz_B = fits["mz_midref"]["p_finest3"]
    pf_mz2_A = fits["mz2_selfref"]["p_finest3"]
    pf_mz2_B = fits["mz2_midref"]["p_finest3"]
    res["headline"] = dict(p_mz_selfref=p_mz_A, p_mz_midref=p_mz_B,
                           p_mz2_selfref=p_mz2_A, p_mz2_midref=p_mz2_B,
                           p_finest3_mz_selfref=pf_mz_A, p_finest3_mz_midref=pf_mz_B,
                           p_finest3_mz2_selfref=pf_mz2_A, p_finest3_mz2_midref=pf_mz2_B)
    p_pool = float(np.median([p_mz_A, p_mz_B, p_mz2_A, p_mz2_B]))
    p_pool_fine = float(np.median([pf_mz_A, pf_mz_B, pf_mz2_A, pf_mz2_B]))
    refs_agree = max(abs(p_mz_A - p_mz_B), abs(p_mz2_A - p_mz2_B)) < 0.5
    # Honest verdict.  H1 expected weak ~1 (falsify 2.0).  We report whichever
    # the data show: weak-2.0 is "falsified" only if the asymptotic order is
    # clearly below ~1.5; if it sits near 2 the cited order HOLDS.
    weak2_falsified = (p_pool_fine < 1.5) and refs_agree
    weak2_holds = (p_pool_fine >= 1.6) and refs_agree
    res["weak_order_pooled_median_full"] = p_pool
    res["weak_order_pooled_median_finest3"] = p_pool_fine
    res["references_agree"] = bool(refs_agree)
    res["decisive_weak_order_2p0_falsified"] = bool(weak2_falsified)
    res["decisive_weak_order_2p0_holds"] = bool(weak2_holds)
    res["runtime_s"] = time.time() - t0
    print(f"\nHEADLINE weak order: pooled median full={p_pool:.3f}, "
          f"finest3={p_pool_fine:.3f}; refs agree={refs_agree}", flush=True)
    print(f"DECISIVE: cited weak-order-2.0 falsified={weak2_falsified}, "
          f"holds={weak2_holds}", flush=True)

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    dd = np.array(dts)
    for ax, gA, gB, lab in [(axes[0], err_mz_A, err_mz_B, "m_z"),
                            (axes[1], err_mz2_A, err_mz2_B, "m_z^2")]:
        ax.loglog(dd, gA, "o-", label=f"self-Cayley ref p={fit_order(dts,gA)[0]:.2f}")
        ax.loglog(dd, gB, "s-", label=f"midpoint ref p={fit_order(dts,gB)[0]:.2f}")
        ax.loglog(dd, gA[0] * (dd / dd[0]) ** 1.0, "k--", alpha=0.6, label="slope 1.0")
        ax.loglog(dd, gA[0] * (dd / dd[0]) ** 2.0, "k:", alpha=0.6, label="slope 2.0")
        ax.set_xlabel("dt [s]"); ax.set_ylabel(f"weak error |E[{lab}]_dt - E[{lab}]_ref|")
        ax.set_title(f"Weak order, g={lab}")
        ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    figpath = os.path.join(OUTDIR, "weak_order.png")
    fig.savefig(figpath, dpi=140)
    res["figure"] = figpath

    jpath = os.path.join(OUTDIR, "weak_order_results.json")
    with open(jpath, "w") as f:
        json.dump(res, f, indent=2, default=float)
    print(f"\nwrote {figpath}\nwrote {jpath}", flush=True)
    return res


if __name__ == "__main__":
    main()
