"""Exp E (strong order) — strong stochastic convergence order of the published
explicit-omega Cayley sLLG step, via common-Brownian-path self-coupling.

Thesis section: §2.2.3.2, Experiment E.

WHAT THIS MEASURES
------------------
§2.2.3.2 writes the update as an *implicit-midpoint* Cayley step and inherits,
by citation [31], strong convergence order 1.0.  The PUBLISHED code
`switching_vector()` instead evaluates omega at the LEFT endpoint m_n and applies
ONE closed-form Cayley rotation (no midpoint iteration).  Hypothesis H1: the
multiplicative noise enters non-commutatively (m x . and m x (m x .) do not
commute) and the left-endpoint sampling carries no Levy-area / Milstein term, so
the STRONG order falls BELOW the cited 1.0 (toward ~0.5 in the noise-dominated
limit).

METHOD (Kloeden-Platen common-Brownian-path coupling)
-----------------------------------------------------
* Per seed draw ONE fine 3-vector Wiener path on the finest grid dt_min (the
  REFERENCE grid, level 0).  Coarse increments are built by SUMMING fine
  increments (variance-conserving), so every level shares ONE realisation.
* The realisation is injected through the verified `h_th_ext` hook:
      h_th_ext = A0 * dW_step / dt_step ,
  A0 = sqrt(2 kB T alpha / (u0 Ms gamma v))   (dt-INDEPENDENT FDT amplitude).
  This reproduces the engine's internal noise exactly, because internally
  H_TH = sqrt(2 kB T alpha /(u0 Ms gamma v dt)) * N(0,1) = A0 * dW/dt with
  dW ~ N(0, dt).  (FDT amplitude convention; verified against anisotropy.field.)
* REFERENCE is STRICTLY FINER than every fitted level (fit l = 1..L,
  dt_l = dt_min * 2^l; reference = level 0) -> avoids the p~0 self-coupling
  degeneracy of "same path is both reference and a fitted point".
* CRITICAL red-team fix: the reference uses an INDEPENDENT higher strong-order
  integrator (true implicit-midpoint Cayley, 3 fixed-point iterations on the
  Stratonovich midpoint), NOT the published scheme itself.  A self-reference
  (published-as-its-own-reference) measures convergence to the SCHEME'S OWN
  dt->0 limit and inflates the apparent order (probe: self-ref gives p~1.0,
  independent-ref gives p~0.7-0.9 at the same operating point).
* Strong error e(dt) = E_seed[ | m_dt(T) - m_ref(T) | ] (terminal Euclidean
  distance).  Order p = slope of log e vs log dt; bootstrap CI over seeds.

GATES (mandatory; abort if any fail)
-------------------------------------
GATE-0 : scalar GBM on the SAME order-fit harness gives EM p ~ 0.5 +- 0.1 AND
         Milstein p ~ 1.0 +- 0.1 (proves the fit resolves 0.5 vs 1.0).
GATE-1 : deterministic (NON=0) published sLLG p_det >= 0.95.
GATE-2 : noise-dominance ratio = (stochastic m-perpendicular rotation per step)
         / (deterministic m-perpendicular rotation per step) >= 50, reached by
         shrinking Ms via Ms_scale in {1e-2, 3e-3, 1e-3}.  The metric uses the
         m-PERPENDICULAR component |omega x m| (the part that actually rotates
         m); |omega| alone is dominated by the m-parallel PMA field and is not a
         rotation.

DECISIVE
--------
Published explicit-Cayley strong p, vs the independent midpoint reference, lands
strictly BELOW the cited 1.0 (upper CI < ~0.95) and is separated from the
deterministic order anchor -> the cited strong order 1.0 is falsified.

Run:  PYTHONPATH=src python scripts/11_strong_order/strong_order_brownian.py
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
M0 = np.array([0.0, 0.0, 1.0])          # on the easy axis: PMA is m-parallel
T_HORIZON = 3.0e-13                       # short, pre-equilibrium horizon
DT_MIN_EXP = 11                           # reference grid = 2^11 steps (midpoint
                                          # reference is 3x cost/step; 2^11 keeps
                                          # the full run to ~minutes while leaving
                                          # the reference >=4x finer than every
                                          # fitted level)
LEVELS = (1, 2, 3, 4, 5)                  # fitted coarse levels (all >= 1)
MS_SCALES = (1e-2, 3e-3, 1e-3)


# ---------------------------------------------------------------------------
# order fit + bootstrap CI on the slope
# ---------------------------------------------------------------------------
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


def bootstrap_slope_ci(dts, per_seed_errs, n_boot=1500, seed=12345):
    rng = np.random.default_rng(seed)
    n_seed = per_seed_errs.shape[0]
    slopes = []
    for _ in range(n_boot):
        idx = rng.integers(0, n_seed, n_seed)
        e = per_seed_errs[idx].mean(axis=0)
        p, _, _ = fit_order(dts, e)
        slopes.append(p)
    slopes = np.array(slopes)
    return (float(np.median(slopes)),
            float(np.percentile(slopes, 2.5)),
            float(np.percentile(slopes, 97.5)))


# ===========================================================================
# GATE-0 : scalar GBM on the same order-fit harness
# ===========================================================================
def gbm_gate(n_seed=4000, mu=1.5, sigma=1.0, X0=1.0, T=1.0,
             dt_min_exp=12, levels=LEVELS):
    """dX = mu X dt + sigma X dW.  Strong reference = exact solution on the same
    fine path.  EM -> p~0.5, Milstein -> p~1.0.  Identical coarse-summing
    structure as the sLLG harness."""
    N_min = 2 ** dt_min_exp
    dt_min = T / N_min
    rng = np.random.default_rng(20240601)
    dts = [dt_min * (2 ** l) for l in levels]
    em = np.zeros((n_seed, len(levels)))
    mil = np.zeros((n_seed, len(levels)))
    for s in range(n_seed):
        dW_fine = rng.normal(0.0, math.sqrt(dt_min), N_min)
        W_T = float(np.sum(dW_fine))
        X_exact = X0 * math.exp((mu - 0.5 * sigma ** 2) * T + sigma * W_T)
        for li, l in enumerate(levels):
            step = 2 ** l
            N = N_min // step
            dt = dt_min * step
            dW = dW_fine[: N * step].reshape(N, step).sum(axis=1)
            xe = xm = X0
            for k in range(N):
                dWk = dW[k]
                xe = xe + mu * xe * dt + sigma * xe * dWk
                xm = (xm + mu * xm * dt + sigma * xm * dWk
                      + 0.5 * sigma * sigma * xm * (dWk * dWk - dt))
            em[s, li] = abs(xe - X_exact)
            mil[s, li] = abs(xm - X_exact)
    p_em, _, r2_em = fit_order(dts, em.mean(axis=0))
    p_mil, _, r2_mil = fit_order(dts, mil.mean(axis=0))
    _, em_lo, em_hi = bootstrap_slope_ci(dts, em)
    _, mil_lo, mil_hi = bootstrap_slope_ci(dts, mil)
    return dict(dts=dts, em_mean=em.mean(axis=0).tolist(),
                mil_mean=mil.mean(axis=0).tolist(),
                p_em=p_em, r2_em=r2_em, em_ci=[em_lo, em_hi],
                p_mil=p_mil, r2_mil=r2_mil, mil_ci=[mil_lo, mil_hi])


# ===========================================================================
# sLLG harness
# ===========================================================================
def make_cc(Ms_scale=1.0):
    cc = PhysicalConstantsConfig()
    cc.Ms = PhysicalConstantsConfig().Ms * Ms_scale
    return cc


def fdt_amplitude(cc, T=None):
    T_use = cc.T if T is None else T
    return math.sqrt(2.0 * cc.kb * T_use * cc.alpha
                     / (cc.u0 * cc.Ms * cc.gamma * cc.v))


def omega_at(m, cc, *, V_MTJ, I_SOT, R_MTJ, ESOT, ESTT, R_SOT_FL_DL,
             VNV, NON, h_th_ext=None, ENE=1, T=None):
    m = np.asarray(m, float)
    theta = float(np.arccos(np.clip(m[2], -1.0, 1.0)))
    phi = float(np.arctan2(m[1], m[0]))
    H_eff, _ = field(theta, phi, V_MTJ, n=1, NON=NON, ENE=ENE, VNV=VNV,
                     constants=cc, demag_mode="ellipsoid", h_th_ext=h_th_ext, T=T)
    H_eff = np.asarray(H_eff, float)
    Ms_use = cc.Ms
    I_MTJ = V_MTJ / R_MTJ if R_MTJ else 0.0
    J_STT = I_MTJ / cc.A1
    J_SOT = I_SOT / cc.A2
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
    H_DL_STT = ESTT * cc.h_bar * cc.P * J_STT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_DL_SOT = ESOT * cc.h_bar * cc.theta_SH * J_SOT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_FL_SOT = R_SOT_FL_DL * H_DL_SOT
    return _omega_from_state(m, H_eff, np.array([-1.0, 0.0, 0.0]),
                             H_DL_SOT=H_DL_SOT, H_FL_SOT=H_FL_SOT,
                             H_DL_STT=H_DL_STT, H_FL_STT=0.0,
                             sigma_STT=np.array([0.0, 0.0, 1.0]),
                             alpha=cc.alpha, gamma_red=gamma_red)


def _h_th(dW, dt, A0):
    return A0 * np.asarray(dW, float) / dt


def step_published(m, cc, dt, dW, drive, NON=1):
    """Published explicit-omega left-endpoint Cayley step."""
    h = _h_th(dW, dt, drive["A0"]) if (NON and dW is not None) else None
    om = omega_at(m, cc, h_th_ext=h, NON=NON, **drive["omega_kw"])
    mn = cayley_step(m, om, dt)
    return mn / np.linalg.norm(mn)


def step_midpoint(m, cc, dt, dW, drive, NON=1, iters=3):
    """Independent higher-strong-order reference: true implicit-midpoint Cayley
    (fixed-point on the Stratonovich midpoint m_{1/2} = (m_n + m_{n+1})/2)."""
    h = _h_th(dW, dt, drive["A0"]) if (NON and dW is not None) else None
    mn = m.copy()
    for _ in range(iters):
        mmid = (m + mn) / 2.0
        mmid = mmid / np.linalg.norm(mmid)
        om = omega_at(mmid, cc, h_th_ext=h, NON=NON, **drive["omega_kw"])
        mn = cayley_step(m, om, dt)
    return mn / np.linalg.norm(mn)


def integrate(cc, dt, N, dW, drive, stepper, NON=1):
    m = np.array(M0, float)
    for k in range(N):
        dWk = None if dW is None else dW[k]
        m = stepper(m, cc, dt, dWk, drive, NON=NON)
    return m


def build_drive(cc, *, ENE):
    """Pure-noise operating point: I_SOT=0, ESOT=0, VNV=0.  ENE selects whether
    the constant -50 Oe in-plane field is present (ENE=1, used for the GATE-2
    finite deterministic-torque ratio) or absent (ENE=0, used for the DECISIVE
    pure-multiplicative-noise strong-order measurement so the only torque is the
    thermal field)."""
    return dict(A0=fdt_amplitude(cc),
                omega_kw=dict(V_MTJ=0.0, I_SOT=0.0, R_MTJ=0.0, ESOT=0, ESTT=0,
                              R_SOT_FL_DL=0.83, VNV=0, ENE=ENE, T=None))


# ---------------------------------------------------------------------------
# GATE-2 noise-dominance diagnostic (m-perpendicular rotation per step)
# ---------------------------------------------------------------------------
def noise_dominance(cc, dt, drive, n_probe=600, seed=7):
    rng = np.random.default_rng(seed)
    A0 = drive["A0"]
    m = np.array(M0, float)
    det_a, sto_a = [], []
    for _ in range(n_probe):
        om_det = omega_at(m, cc, h_th_ext=None, NON=0, **drive["omega_kw"])
        det_a.append(np.linalg.norm(np.cross(om_det, m)) * dt)
        dW = rng.normal(0.0, math.sqrt(dt), 3)
        om_sto = omega_at(m, cc, h_th_ext=A0 * dW / dt, NON=1, **drive["omega_kw"])
        sto_a.append(np.linalg.norm(np.cross(om_sto - om_det, m)) * dt)
    det = float(np.mean(det_a))
    sto = float(np.mean(sto_a))
    return (sto / det if det > 0 else float("inf")), det, sto


# ===========================================================================
# strong-order driver (common Brownian path; independent midpoint reference)
# ===========================================================================
def strong_order(cc, drive, *, fit_stepper, ref_stepper, NON=1,
                 T_horizon=T_HORIZON, dt_min_exp=DT_MIN_EXP, levels=LEVELS,
                 n_seed=200, base_seed=11):
    N_min = 2 ** dt_min_exp
    dt_min = T_horizon / N_min
    dts = [dt_min * (2 ** l) for l in levels]
    errs = np.zeros((n_seed, len(levels)))
    mr = np.random.default_rng(base_seed)
    for s in range(n_seed):
        rng = np.random.default_rng(int(mr.integers(0, 2 ** 31 - 1)))
        dW_fine = (rng.normal(0.0, math.sqrt(dt_min), (N_min, 3))
                   if NON else None)
        ref = integrate(cc, dt_min, N_min, dW_fine, drive, ref_stepper, NON=NON)
        for li, l in enumerate(levels):
            step = 2 ** l
            N = N_min // step
            dt = dt_min * step
            dW_l = (dW_fine[: N * step].reshape(N, step, 3).sum(axis=1)
                    if NON else None)
            mf = integrate(cc, dt, N, dW_l, drive, fit_stepper, NON=NON)
            errs[s, li] = float(np.linalg.norm(mf - ref))
    mean_err = errs.mean(axis=0)
    p, b, r2 = fit_order(dts, mean_err)
    pmed, plo, phi = bootstrap_slope_ci(dts, errs)
    return dict(dts=dts, mean_err=mean_err.tolist(), p=p, r2=r2,
                p_ci=[plo, phi], p_med=pmed)


def main():
    t0 = time.time()
    res = {"experiment": "E_strong_order", "section": "2.2.3.2",
           "horizon_s": T_HORIZON, "dt_min_exp": DT_MIN_EXP,
           "levels": list(LEVELS), "Ms_scales": list(MS_SCALES),
           "m0": M0.tolist()}

    # -------- GATE-0 : GBM EM vs Milstein --------
    print("[GATE-0] scalar GBM EM vs Milstein ...", flush=True)
    g0 = gbm_gate()
    res["gate0_gbm"] = g0
    print(f"  EM p={g0['p_em']:.3f} CI={g0['em_ci']} R2={g0['r2_em']:.4f}", flush=True)
    print(f"  Milstein p={g0['p_mil']:.3f} CI={g0['mil_ci']} R2={g0['r2_mil']:.4f}", flush=True)
    gate0 = abs(g0["p_em"] - 0.5) <= 0.1 and abs(g0["p_mil"] - 1.0) <= 0.1
    res["gate0_pass"] = bool(gate0)
    print(f"  GATE-0 {'PASS' if gate0 else 'FAIL'}", flush=True)
    if not gate0:
        raise SystemExit("GATE-0 failed: harness cannot resolve 0.5 vs 1.0.")

    # -------- GATE-1 : deterministic published order --------
    print("[GATE-1] deterministic (NON=0) published-step order ...", flush=True)
    cc1 = make_cc(1.0)
    dr1 = build_drive(cc1, ENE=1)
    g1 = strong_order(cc1, dr1, fit_stepper=step_published,
                      ref_stepper=step_published, NON=0, n_seed=1)
    res["gate1_det"] = dict(dts=g1["dts"], mean_err=g1["mean_err"],
                            p=g1["p"], r2=g1["r2"])
    print(f"  p_det={g1['p']:.3f} R2={g1['r2']:.4f}", flush=True)
    gate1 = g1["p"] >= 0.95
    res["gate1_pass"] = bool(gate1)
    print(f"  GATE-1 {'PASS' if gate1 else 'FAIL'}", flush=True)
    if not gate1:
        raise SystemExit("GATE-1 failed: deterministic order < 0.95.")

    # -------- GATE-2 : noise dominance over Ms_scale grid --------
    print("[GATE-2] noise-dominance ratio (perp rotation) over Ms grid ...", flush=True)
    dt_probe = T_HORIZON / (2 ** DT_MIN_EXP) * (2 ** max(LEVELS))
    gate2 = {}
    passing = []
    for sc in MS_SCALES:
        cc = make_cc(sc)
        dr = build_drive(cc, ENE=1)   # ENE=1 -> finite det perp torque
        ratio, det, sto = noise_dominance(cc, dt_probe, dr)
        gate2[f"{sc:g}"] = dict(ratio=ratio, det_angle=det, sto_angle=sto)
        print(f"  Ms_scale={sc:g}: ratio={ratio:.1f} (det={det:.2e} sto={sto:.2e})", flush=True)
        if ratio >= 50:
            passing.append(sc)
    res["gate2"] = gate2
    gate2_pass = len(passing) > 0
    res["gate2_pass"] = bool(gate2_pass)
    res["gate2_passing_scales"] = passing
    print(f"  GATE-2 {'PASS' if gate2_pass else 'FAIL'} (scales>=50: {passing})", flush=True)
    if not gate2_pass:
        raise SystemExit("GATE-2 failed: no Ms_scale reaches noise/det >= 50.")

    # -------- DECISIVE : published strong order vs independent midpoint ref --------
    print("[DECISIVE] published explicit-Cayley strong order "
          "(independent midpoint reference, pure noise) ...", flush=True)
    n_seed = 400
    stoch = {}
    for sc in MS_SCALES:
        cc = make_cc(sc)
        dr = build_drive(cc, ENE=0)   # pure multiplicative noise
        r = strong_order(cc, dr, fit_stepper=step_published,
                         ref_stepper=step_midpoint, NON=1,
                         n_seed=n_seed, base_seed=int(1000 + round(1e6 * sc)))
        ratio = gate2[f"{sc:g}"]["ratio"]
        # asymptotic (finest-3) sub-window slope: least contaminated by the
        # large-rotation saturation that steepens the coarse-dt points.
        p_fine, _, _ = fit_order(r["dts"][:3], r["mean_err"][:3])
        stoch[f"{sc:g}"] = dict(dts=r["dts"], mean_err=r["mean_err"],
                                p=r["p"], r2=r["r2"], p_ci=r["p_ci"],
                                p_med=r["p_med"], p_finest3=p_fine, ratio=ratio)
        print(f"  Ms_scale={sc:g} (ratio={ratio:.0f}): strong p_full={r['p']:.3f} "
              f"CI={r['p_ci']} R2={r['r2']:.4f}  p_finest3={p_fine:.3f}", flush=True)
    res["stochastic"] = stoch

    # control: published-as-its-own-reference (the contaminated self-coupling)
    print("[CONTROL] self-reference (published vs published) — "
          "demonstrates the degeneracy the independent ref avoids ...", flush=True)
    cc = make_cc(MS_SCALES[-1])
    dr = build_drive(cc, ENE=0)
    self_ref = strong_order(cc, dr, fit_stepper=step_published,
                            ref_stepper=step_published, NON=1,
                            n_seed=250, base_seed=55)
    res["self_reference_control"] = dict(
        Ms_scale=f"{MS_SCALES[-1]:g}", p=self_ref["p"], p_ci=self_ref["p_ci"],
        r2=self_ref["r2"], mean_err=self_ref["mean_err"], dts=self_ref["dts"])
    print(f"  self-ref p={self_ref['p']:.3f} CI={self_ref['p_ci']} "
          f"(vs independent-ref) — inflation diagnostic", flush=True)

    # headline = CLEANEST scale = highest R2 (least saturated).  Stronger noise
    # (smaller Ms) drives the coarse-dt points into the large-rotation
    # saturation regime, which steepens and biases the slope; the cleanest fit
    # is the honest asymptotic measurement.
    head = max(stoch, key=lambda k: stoch[k]["r2"])
    p_strong = stoch[head]["p"]
    p_ci = stoch[head]["p_ci"]
    p_strong_fine = stoch[head]["p_finest3"]
    res["headline"] = dict(Ms_scale=head, p_full=p_strong, p_ci=p_ci,
                           p_finest3=p_strong_fine,
                           p_det=g1["p"], p_milstein_anchor=g0["p_mil"],
                           p_self_reference=self_ref["p"])

    # honest verdict.  H1 expected strong ~0.5 (falsify 1.0).  Report what the
    # data show: the cited order 1.0 is "falsified" only if the clean asymptotic
    # order is clearly below ~0.8; if it approaches ~1.0 the cited order HOLDS.
    strong1_falsified = p_strong_fine < 0.8 and p_ci[1] < 0.8
    strong1_holds = p_strong_fine >= 0.85
    res["decisive_strong_order_1p0_falsified"] = bool(strong1_falsified)
    res["decisive_strong_order_1p0_holds"] = bool(strong1_holds)
    res["runtime_s"] = time.time() - t0
    print(f"\nHEADLINE (cleanest Ms={head}): strong p_full={p_strong:.3f} "
          f"CI={p_ci}, p_finest3={p_strong_fine:.3f}; "
          f"p_det={g1['p']:.3f}, Milstein-anchor={g0['p_mil']:.3f}, "
          f"self-ref={self_ref['p']:.3f}", flush=True)
    print(f"DECISIVE: cited strong-order-1.0 falsified={strong1_falsified}, "
          f"holds={strong1_holds}", flush=True)

    # -------- figure --------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    d0 = np.array(g0["dts"])
    ax.loglog(d0, g0["em_mean"], "o-", label=f"GBM EM p={g0['p_em']:.2f}")
    ax.loglog(d0, g0["mil_mean"], "s-", label=f"GBM Milstein p={g0['p_mil']:.2f}")
    ax.loglog(d0, g0["em_mean"][0] * (d0 / d0[0]) ** 0.5, "k--", alpha=0.5, label="slope 0.5")
    ax.loglog(d0, g0["mil_mean"][0] * (d0 / d0[0]) ** 1.0, "k:", alpha=0.5, label="slope 1.0")
    ax.set_xlabel("dt"); ax.set_ylabel("strong error |X_dt(T)-X_exact|")
    ax.set_title("GATE-0: harness resolves 0.5 vs 1.0")
    ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)

    ax = axes[1]
    for sc in MS_SCALES:
        k = f"{sc:g}"
        d = np.array(stoch[k]["dts"])
        ax.loglog(d, stoch[k]["mean_err"], "o-",
                  label=f"Ms={sc:g} (r={stoch[k]['ratio']:.0f}) p={stoch[k]['p']:.2f}")
    dd = np.array(stoch[head]["dts"])
    e0 = stoch[head]["mean_err"][0]
    ax.loglog(dd, e0 * (dd / dd[0]) ** 0.5, "k--", alpha=0.6, label="slope 0.5")
    ax.loglog(dd, e0 * (dd / dd[0]) ** 1.0, "k:", alpha=0.6, label="slope 1.0")
    ax.set_xlabel("dt [s]"); ax.set_ylabel("strong error |m_dt(T)-m_ref(T)|")
    ax.set_title("DECISIVE: published Cayley strong order (indep. midpoint ref)")
    ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    figpath = os.path.join(OUTDIR, "strong_order.png")
    fig.savefig(figpath, dpi=140)
    res["figure"] = figpath

    jpath = os.path.join(OUTDIR, "strong_order_results.json")
    with open(jpath, "w") as f:
        json.dump(res, f, indent=2, default=float)
    print(f"\nwrote {figpath}\nwrote {jpath}", flush=True)
    return res


if __name__ == "__main__":
    main()
