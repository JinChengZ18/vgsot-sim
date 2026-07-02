#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Experiment B — Boltzmann effective-temperature dt-sweep  (thesis §2.2.3.2, claim SC5)
=====================================================================================

Tests the chapter's strongest theoretical claim:

    SC5 : the explicit-omega Cayley step samples the *correct* Boltzmann
          distribution  p(m_z) ∝ exp(-Δ (1 - m_z²))  with NO Itô↔Stratonovich
          drift correction, i.e. T_eff/T → 1 as Δt → 0.

The hypothesis under attack (H1) is that the *published* kernel evaluates ω at
the LEFT end point m_n, which breaks the time-symmetry that the implicit-midpoint
convergence theorem relies on, so the stationary distribution is *over-heated*
(T_eff/T > 1, growing with Δt).  A genuine implicit-midpoint Cayley step is
included as the control that should restore T_eff/T → 1.

Physical set-up (zero-drive equilibrium, pure uniaxial well)
------------------------------------------------------------
* I_SOT = 0, V_MTJ = 0, ESOT = 0, ESTT = 0  → no torque drive.
* ENE = 0                                    → the −50 Oe in-plane h_ex bias is OFF.
* VNV = 0                                    → VCMA off.
* NON = 1                                    → Brown-1963 / FDT thermal field ON.

With ENE = VNV = 0 the remaining energy is PMA + demag.  Because the ellipsoid
demag has N_x = N_y, the in-plane part is isotropic, so the total energy is
*exactly* uniaxial,

    E(m) = -K_eff · v · m_z²  + const ,
    K_eff = Ki/tf - 0.5·u0·Ms²·(N_z - N_x) ,

and the equilibrium marginal is analytic,

    p(m_z) ∝ exp(-Δ (1 - m_z²)) ,   Δ = K_eff · v / (k_B T)  ≈ 48.515 .

(The FDT amplitude in anisotropy.field() is exactly the Brown-1963 standard
form whose Stratonovich continuous limit is exp(-E/kT); see plan §2 ref.)

Red-team must-fixes implemented
-------------------------------
1. **No harmonic estimator.**  T_eff is read from (a) an EXACT bounded-basin
   moment inversion <1-m_z²> = f(Δ) [single-basin theory, primary], and (b) a
   WEIGHTED-LS slope fit of ln p(m_z) vs (1-m_z²) over a high-m_z window
   [cross-check].  The harmonic Δ_eff = 1/<1-m_z²> is computed ONLY to print
   its analytic construction bias (0.989 at Δ=48.5) and is never used for the
   verdict.
2. **Fixed PHYSICAL sample time** T_phys = n_steps·Δt across the dt sweep
   (n_steps = round(T_phys/Δt)), NOT fixed step count.  The FDT amplitude
   auto-rescales ∝ 1/√Δt so the Brownian increment scales ∝ √Δt (Euler-Maruyama).
3. **Stationarity proven first.**  τ_int of (1-m_z²) and a trajectory-length
   plateau of T_eff/T are reported; the burn-in is discarded.
4. **Δt→0 Richardson extrapolation** of T_eff/T with block-bootstrap CIs.
   SC5 is declared FALSE iff the Δt→0 intercept differs from 1 by ≥3σ OR the
   slope d(T_eff/T)/d(Δt) > 0 at ≥3σ.

Run modes
---------
    PYTHONPATH=src python scripts/11_boltzmann_teff/teff_dtsweep.py            # PILOT
    PYTHONPATH=src python scripts/11_boltzmann_teff/teff_dtsweep.py --full     # FULL

Outputs (under result/sec_2_2_3_2/B/):
    teff_dtsweep_pilot.json / teff_dtsweep_full.json
    teff_dtsweep_pilot.png  / teff_dtsweep_full.png
    teff_stationarity_pilot.png / teff_stationarity_full.png
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.optimize import brentq

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.anisotropy import field
from vgsot_sim.demag import demag_factors
from vgsot_sim.dynamic_switching_vector import _omega_from_state, cayley_step
from vgsot_sim.dynamic_switching import switching


# ─────────────────────────────────────────────────────────────────────────
# Analytic constants
# ─────────────────────────────────────────────────────────────────────────
def analytic_delta(cc: PhysicalConstantsConfig, demag_mode: str = "ellipsoid") -> float:
    """Δ = K_eff·v/(k_B T) for the pure uniaxial well (PMA + isotropic-in-plane demag)."""
    Nx, Ny, Nz = demag_factors(cc, mode=demag_mode)
    Ku = cc.Ki / cc.tf                                  # PMA uniaxial density
    K_shape = 0.5 * cc.u0 * cc.Ms ** 2 * (Nz - Nx)      # demag out-of-plane penalty
    K_eff = Ku - K_shape
    return K_eff * cc.v / (cc.kb * cc.T)


# ─────────────────────────────────────────────────────────────────────────
# T_eff estimators (NEVER the harmonic 1/<1-m_z^2>)
# ─────────────────────────────────────────────────────────────────────────
def _mean_u_of_delta(D: float) -> float:
    """Exact single-basin <u>, u = 1-m_z², under p(m_z) ∝ exp(-D u), m_z∈[0,1]."""
    num, _ = quad(lambda mz: (1.0 - mz * mz) * np.exp(-D * (1.0 - mz * mz)), 0.0, 1.0)
    den, _ = quad(lambda mz: np.exp(-D * (1.0 - mz * mz)), 0.0, 1.0)
    return num / den


def delta_from_moment(u_samples: np.ndarray) -> float:
    """EXACT bounded-basin moment inversion: solve <u>_theory(Δ) = mean(u).

    Construction-bias-free (recovers Δ exactly on analytic samples).
    """
    mbar = float(np.mean(u_samples))
    # bracket generously; mean_u is monotone decreasing in D
    return float(brentq(lambda D: _mean_u_of_delta(D) - mbar, 1.0, 5000.0, xtol=1e-6))


def delta_from_wls_slope(mz: np.ndarray, mz_window: float = 0.85, n_bins: int = 25):
    """WEIGHTED-LS slope of ln p(m_z) vs u=(1-m_z²) over the high-m_z window.

    Slope of ln p vs u is exactly -Δ analytically → unbiased by construction
    (small residual histogram-binning bias only).  Returns (Δ_fit, r2, n_in_window).
    """
    mz_pos = np.abs(mz)                       # fold both basins onto m_z∈[0,1]
    sel = mz_pos >= mz_window
    if sel.sum() < 4 * n_bins:
        return np.nan, np.nan, int(sel.sum())
    u = 1.0 - mz_pos[sel] ** 2
    edges = np.linspace(u.min(), u.max(), n_bins + 1)
    cnt, _ = np.histogram(u, bins=edges)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]
    good = cnt > 0
    if good.sum() < 3:
        return np.nan, np.nan, int(sel.sum())
    x = ctr[good]
    y = np.log(cnt[good] / width)
    w = cnt[good].astype(float)               # Poisson weights ∝ counts
    Sw = w.sum(); Swx = (w * x).sum(); Swy = (w * y).sum()
    Swxx = (w * x * x).sum(); Swxy = (w * x * y).sum()
    denom = Sw * Swxx - Swx ** 2
    slope = (Sw * Swxy - Swx * Swy) / denom
    intercept = (Swy - slope * Swx) / Sw
    yhat = slope * x + intercept
    ss_res = (w * (y - yhat) ** 2).sum()
    ybar = Swy / Sw
    ss_tot = (w * (y - ybar) ** 2).sum()
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(-slope), float(r2), int(sel.sum())


# ─────────────────────────────────────────────────────────────────────────
# omega path (verified harness; explicit-omega + true-midpoint share this)
# ─────────────────────────────────────────────────────────────────────────
def omega_zero_drive(m, cc, h_th):
    """ω for zero-drive pure-uniaxial well with externally injected thermal field h_th."""
    m = np.asarray(m, float)
    theta = float(np.arccos(np.clip(m[2], -1.0, 1.0)))
    phi = float(np.arctan2(m[1], m[0]))
    H_eff, _ = field(theta, phi, 0.0, n=1, NON=1, ENE=0, VNV=0,
                     constants=cc, demag_mode="ellipsoid", h_th_ext=h_th)
    H_eff = np.asarray(H_eff, float)
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
    return _omega_from_state(m, H_eff, np.array([-1.0, 0.0, 0.0]),
                             H_DL_SOT=0.0, H_FL_SOT=0.0, H_DL_STT=0.0, H_FL_STT=0.0,
                             sigma_STT=np.array([0.0, 0.0, 1.0]),
                             alpha=cc.alpha, gamma_red=gamma_red)


def fdt_amplitude(cc, dt):
    """Brown-1963 FDT thermal-field amplitude per unit N(0,1) at step size dt."""
    return np.sqrt(2.0 * cc.kb * cc.T * cc.alpha /
                   (cc.u0 * cc.Ms * cc.gamma * cc.v * dt))


def _make_eq_sampler(Delta, n_grid=400001):
    """Return a closure that draws m0 from the ANALYTIC Boltzmann distribution.

    Starting every trajectory at the *target* analytic equilibrium removes the
    cold-pole transient (a single trajectory from the pole needs ≫ τ_int to
    relax). If the integrator's true stationary law is hotter than analytic,
    eq-started trajectories drift outward and T_eff/T climbs toward the
    integrator value from below — so the length-plateau is the convergence
    diagnostic, NOT a fixed burn-in.
    """
    mzg = np.linspace(0.0, 1.0, n_grid)
    w = np.exp(-Delta * (1.0 - mzg ** 2))
    w /= w.sum()
    cdf = np.cumsum(w)

    def sample(rng):
        mz = mzg[np.searchsorted(cdf, rng.random())]
        sign = 1.0 if rng.random() < 0.5 else -1.0    # either basin (symmetric)
        phi = rng.uniform(0.0, 2.0 * np.pi)
        st = np.sqrt(max(1.0 - mz ** 2, 0.0))
        return np.array([st * np.cos(phi), st * np.sin(phi), sign * mz])

    return sample


# ─────────────────────────────────────────────────────────────────────────
# Single trajectory for one integrator (shared injected Brownian path)
# ─────────────────────────────────────────────────────────────────────────
def run_trajectory(cc, dt, n_steps, m0, rng, integrator, n_mid=3):
    """Return m_z time series (length n_steps) for the chosen integrator.

    The SAME thermal-field realisation h_th is used at every ω-evaluation
    within a step (so the true-midpoint scheme reuses the step's noise rather
    than redrawing it — the physically-correct Stratonovich midpoint).
    """
    amp = fdt_amplitude(cc, dt)
    m = np.asarray(m0, float).copy()
    m /= np.linalg.norm(m)
    mz = np.empty(n_steps, dtype=float)

    if integrator == "euler_spherical":
        theta = float(np.arccos(np.clip(m[2], -1.0, 1.0)))
        phi = float(np.arctan2(m[1], m[0]))
        for i in range(n_steps):
            h_th = amp * rng.normal(0.0, 1.0, 3)
            mz_i, phi, theta = switching(
                V_MTJ=0.0, I_SOT=0.0, R_MTJ=1.0, theta=theta, phi=phi,
                ESTT=0, ESOT=0, VNV=0, NON=1, constants=cc, h_th_ext=h_th)
            mz[i] = mz_i
        return mz

    if integrator == "cayley_explicit":
        for i in range(n_steps):
            h_th = amp * rng.normal(0.0, 1.0, 3)
            om = omega_zero_drive(m, cc, h_th)
            m = cayley_step(m, om, dt)
            m /= np.linalg.norm(m)
            mz[i] = m[2]
        return mz

    if integrator == "cayley_true_midpoint":
        for i in range(n_steps):
            h_th = amp * rng.normal(0.0, 1.0, 3)
            m_next = m.copy()                       # initial guess = m_n
            for _ in range(n_mid):
                m_mid = 0.5 * (m + m_next)
                nrm = np.linalg.norm(m_mid)
                if nrm > 0:
                    m_mid = m_mid / nrm
                om = omega_zero_drive(m_mid, cc, h_th)
                m_next = cayley_step(m, om, dt)      # exact orthogonal (no renorm)
            m = m_next / np.linalg.norm(m_next)
            mz[i] = m[2]
        return mz

    raise ValueError(f"unknown integrator {integrator!r}")


# ─────────────────────────────────────────────────────────────────────────
# Stationarity diagnostics
# ─────────────────────────────────────────────────────────────────────────
def integrated_autocorr_time(x, c=5.0, max_lag=None):
    """τ_int (in steps) of series x via the automatic-windowing (Sokal) method."""
    x = np.asarray(x, float)
    n = len(x)
    x = x - x.mean()
    if max_lag is None:
        max_lag = n // 2
    # autocovariance via FFT
    f = np.fft.rfft(x, n=2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n].real
    acf /= acf[0]
    tau = 1.0
    for w in range(1, max_lag):
        tau = 1.0 + 2.0 * np.sum(acf[1:w + 1])
        if w >= c * tau:
            break
    return float(tau), acf


# ─────────────────────────────────────────────────────────────────────────
# Block bootstrap over trajectories
# ─────────────────────────────────────────────────────────────────────────
def block_bootstrap_teff(per_traj_u, per_traj_mz, Delta, n_boot, rng,
                         mz_window=0.85, n_bins=25):
    """Bootstrap T_eff/T over trajectories (each trajectory is one block).

    Returns dicts for the moment-inversion estimator and the WLS-slope estimator,
    each with point estimate (from the full pooled sample) and bootstrap σ / CI.
    """
    n_traj = len(per_traj_u)

    # point estimates on the full pooled sample
    u_all = np.concatenate(per_traj_u)
    mz_all = np.concatenate(per_traj_mz)
    D_mom = delta_from_moment(u_all)
    D_slo, r2_slo, n_win = delta_from_wls_slope(mz_all, mz_window, n_bins)
    teff_mom_pt = Delta / D_mom
    teff_slo_pt = Delta / D_slo if np.isfinite(D_slo) else np.nan

    boot_mom = np.empty(n_boot)
    boot_slo = np.empty(n_boot)
    idx_all = np.arange(n_traj)
    for b in range(n_boot):
        pick = rng.choice(idx_all, size=n_traj, replace=True)
        u_b = np.concatenate([per_traj_u[k] for k in pick])
        mz_b = np.concatenate([per_traj_mz[k] for k in pick])
        boot_mom[b] = Delta / delta_from_moment(u_b)
        Ds, _, _ = delta_from_wls_slope(mz_b, mz_window, n_bins)
        boot_slo[b] = Delta / Ds if np.isfinite(Ds) else np.nan

    def summarise(point, boot):
        boot = boot[np.isfinite(boot)]
        return dict(point=float(point),
                    sigma=float(np.std(boot, ddof=1)),
                    ci_lo=float(np.percentile(boot, 2.5)),
                    ci_hi=float(np.percentile(boot, 97.5)),
                    n_boot=int(len(boot)))

    return (dict(estimator="moment_inversion", **summarise(teff_mom_pt, boot_mom)),
            dict(estimator="wls_slope", r2=float(r2_slo), n_window=int(n_win),
                 **summarise(teff_slo_pt, boot_slo)))


# ─────────────────────────────────────────────────────────────────────────
# Main sweep
# ─────────────────────────────────────────────────────────────────────────
def run_sweep(cfg, checkpoint_path=None):
    # Zero the −50 Oe in-plane exchange bias so the well is PURELY uniaxial and
    # the analytic Boltzmann p(m_z) ∝ exp(-Δ(1-m_z²)) is exact.  The Cayley path
    # uses ENE=0 anyway, but the legacy spherical-Euler `switching()` hardcodes
    # ENE=1, so we must zero h_ex in the config for a fair common-well anchor.
    cc = dataclasses.replace(PhysicalConstantsConfig(),
                             h_ex_x=0.0, h_ex_y=0.0, h_ex_z=0.0)
    Delta = analytic_delta(cc)

    # harmonic-estimator construction bias (reported, NEVER used for verdict)
    mz_grid = np.linspace(0.0, 1.0, 2_000_001)
    w = np.exp(-Delta * (1.0 - mz_grid ** 2)); w /= w.sum()
    mean_u_analytic = float(np.sum((1.0 - mz_grid ** 2) * w))
    harmonic_bias = (1.0 / mean_u_analytic) / Delta   # 0.989 → harmonic over-reports T_eff/T

    integrators = cfg["integrators"]
    integ_idx = {name: i for i, name in enumerate(integrators)}
    dt_list = cfg["dt_list_ps"]
    T_phys = cfg["T_phys_s"]
    burn_frac = cfg["burn_frac"]
    n_traj = cfg["n_traj"]
    n_boot = cfg["n_boot"]
    base_seed = cfg["base_seed"]
    mz_window = cfg["mz_window"]
    n_bins = cfg["n_bins"]

    # Each trajectory starts at an INDEPENDENT analytic-equilibrium draw — no
    # cold-pole transient. A small fixed burn-in absorbs the discretization-
    # specific relaxation onset; the length-plateau confirms stationarity.
    eq_sample = _make_eq_sampler(Delta)

    results = {"meta": dict(Delta=Delta, mean_u_analytic=mean_u_analytic,
                            harmonic_construction_bias=harmonic_bias,
                            T_phys_s=T_phys, n_traj=n_traj, n_boot=n_boot,
                            burn_frac=burn_frac, dt_list_ps=dt_list,
                            integrators=integrators, mz_window=mz_window,
                            n_bins=n_bins, init="analytic_equilibrium"),
               "sweep": {}}

    stationarity = {}     # per integrator at smallest dt: tau_int, plateau curve

    t_start = time.time()
    for integ in integrators:
        results["sweep"][integ] = []
        for dt_ps in dt_list:
            dt = dt_ps * 1e-12
            cc_dt = dataclasses.replace(cc, t_step=dt)
            n_steps = int(round(T_phys / dt))
            n_burn = int(round(burn_frac * n_steps))

            per_traj_u = []
            per_traj_mz = []
            plateau_keep = []         # full mz of the first few trajs (plateau only)
            n_plateau = min(8, n_traj)
            tau_acc = []
            for k in range(n_traj):
                rng = np.random.default_rng(base_seed + 100003 * k +
                                            7919 * integ_idx[integ] +
                                            int(round(dt_ps * 1000)))
                m0 = eq_sample(rng)
                mz_ts = run_trajectory(cc_dt, dt, n_steps, m0, rng, integ)
                mz_eq = mz_ts[n_burn:]
                per_traj_u.append(1.0 - mz_eq ** 2)
                per_traj_mz.append(mz_eq)
                if dt_ps == min(dt_list) and k < n_plateau:
                    plateau_keep.append(mz_ts)    # retain only a few, only at min dt
                tau, _ = integrated_autocorr_time(1.0 - mz_eq ** 2)
                tau_acc.append(tau)

            # POOLED length-plateau (over the retained trajectories) at the
            # smallest dt: T_eff/T computed on a growing window from each start.
            if dt_ps == min(dt_list):
                fracs = np.linspace(0.1, 1.0, 10)
                plateau = []
                for fr in fracs:
                    cut = max(int(fr * n_steps), n_bins * 4)
                    uu = np.concatenate([1.0 - mz[:cut] ** 2
                                         for mz in plateau_keep])
                    plateau.append(Delta / delta_from_moment(uu))
                stationarity[integ] = dict(
                    plateau=dict(fracs=fracs.tolist(),
                                 t_ns=[fr * T_phys * 1e9 for fr in fracs],
                                 teff=plateau),
                    tau_int_ps=float(np.mean(tau_acc) * dt_ps))

            mom, slo = block_bootstrap_teff(per_traj_u, per_traj_mz, Delta,
                                            n_boot, np.random.default_rng(base_seed + 1),
                                            mz_window, n_bins)
            entry = dict(dt_ps=dt_ps, n_steps=n_steps, n_burn=n_burn,
                         tau_int_steps=float(np.mean(tau_acc)),
                         tau_int_ps=float(np.mean(tau_acc) * dt_ps),
                         n_eff_per_traj=float((n_steps - n_burn) / max(np.mean(tau_acc), 1.0)),
                         moment=mom, slope=slo)
            results["sweep"][integ].append(entry)
            print(f"[{integ:>22}] dt={dt_ps:5.3f}ps n_steps={n_steps:>7} "
                  f"tau_int={entry['tau_int_ps']:.3f}ps  "
                  f"T_eff/T(mom)={mom['point']:.4f}±{mom['sigma']:.4f}  "
                  f"T_eff/T(slope)={slo['point']:.4f}±{slo['sigma']:.4f}",
                  flush=True)
            if checkpoint_path is not None:
                # crash/reboot resilience: persist every finished cell so a
                # killed run loses at most one (integrator, dt) cell
                with open(checkpoint_path, "w", encoding="utf-8") as f:
                    json.dump(results, f, indent=2)

    results["stationarity"] = stationarity

    # ── Δt→0 Richardson (linear-in-dt) extrapolation per integrator ──
    extrap = {}
    for integ in integrators:
        rows = results["sweep"][integ]
        x = np.array([r["dt_ps"] for r in rows])
        for est_key in ("moment", "slope"):
            y = np.array([r[est_key]["point"] for r in rows])
            s = np.array([max(r[est_key]["sigma"], 1e-9) for r in rows])
            ok = np.isfinite(y)
            if ok.sum() >= 2:
                # weighted linear fit y = a + b*dt
                W = 1.0 / s[ok] ** 2
                X = np.vstack([np.ones(ok.sum()), x[ok]]).T
                WX = X * W[:, None]
                cov = np.linalg.inv(X.T @ WX)
                beta = cov @ (WX.T @ y[ok])
                a, b = beta
                sa = np.sqrt(cov[0, 0]); sb = np.sqrt(cov[1, 1])
                extrap.setdefault(integ, {})[est_key] = dict(
                    intercept=float(a), intercept_sigma=float(sa),
                    slope=float(b), slope_sigma=float(sb),
                    intercept_dev_sigma=float((a - 1.0) / sa),
                    slope_pos_sigma=float(b / sb))
    results["extrapolation"] = extrap

    # ── Verdict (red-team decisive criterion) applied to cayley_explicit ──
    published = "cayley_explicit"
    pub = extrap.get(published, {}).get("moment", {})
    sc5_false = False
    reasons = []
    if pub:
        if abs(pub["intercept_dev_sigma"]) >= 3.0:
            sc5_false = True
            reasons.append(f"dt->0 intercept {pub['intercept']:.3f} deviates from 1 "
                           f"by {pub['intercept_dev_sigma']:.1f} sigma")
        if pub["slope_pos_sigma"] >= 3.0:
            sc5_false = True
            reasons.append(f"d(T_eff/T)/d(dt) = {pub['slope']:.3f}/ps > 0 "
                           f"at {pub['slope_pos_sigma']:.1f} sigma")
    results["verdict"] = dict(sc5_false=bool(sc5_false), reasons=reasons,
                              published_integrator=published)

    results["meta"]["wall_time_s"] = time.time() - t_start
    return results


# ─────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────
def plot_results(results, out_png, tag):
    integrators = results["meta"]["integrators"]
    Delta = results["meta"]["Delta"]
    colors = {"cayley_explicit": "#c0392b", "euler_spherical": "#2980b9",
              "cayley_true_midpoint": "#27ae60"}
    labels = {"cayley_explicit": "Cayley explicit-ω (published)",
              "euler_spherical": "Euler spherical (anchor)",
              "cayley_true_midpoint": "Cayley true midpoint (control)"}

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: ln p(m_z) vs (1-m_z²) at the smallest dt, plus analytic line.
    # Reconstruct a representative histogram from the slope-window stats is not
    # stored; instead we draw the analytic reference and annotate fitted Δ.
    u = np.linspace(0.0, 0.25, 200)
    axL.plot(u, -Delta * u, "k--", lw=2, label=f"analytic slope -Δ={Delta:.2f}")
    for integ in integrators:
        row0 = results["sweep"][integ][0]   # smallest dt
        D_slo = Delta / row0["slope"]["point"] if np.isfinite(row0["slope"]["point"]) else np.nan
        if np.isfinite(D_slo):
            axL.plot(u, -D_slo * u, color=colors[integ], lw=1.6,
                     label=f"{labels[integ]}: Δ_fit={D_slo:.2f}")
    axL.set_xlabel(r"$1 - m_z^2$")
    axL.set_ylabel(r"$\ln p(m_z)$ (slope = $-\Delta$)")
    axL.set_title(f"Boltzmann slope ({tag}, dt={results['meta']['dt_list_ps'][0]} ps)")
    axL.legend(fontsize=8)
    axL.grid(alpha=0.3)

    # Right: construction-bias-corrected T_eff/T vs dt, three integrators.
    for integ in integrators:
        rows = results["sweep"][integ]
        x = [r["dt_ps"] for r in rows]
        y = [r["moment"]["point"] for r in rows]
        e = [r["moment"]["sigma"] for r in rows]
        axR.errorbar(x, y, yerr=e, marker="o", capsize=3, color=colors[integ],
                     label=labels[integ])
        ex = results["extrapolation"].get(integ, {}).get("moment")
        if ex:
            xx = np.linspace(0, max(x) * 1.05, 50)
            axR.plot(xx, ex["intercept"] + ex["slope"] * xx, ":",
                     color=colors[integ], alpha=0.7)
            axR.scatter([0], [ex["intercept"]], marker="*", s=120,
                        color=colors[integ], zorder=5,
                        edgecolor="k", linewidth=0.5)
    axR.axhline(1.0, color="k", ls="-", lw=1, alpha=0.6)
    axR.set_xlabel(r"$\Delta t$ (ps)")
    axR.set_ylabel(r"$T_{\rm eff}/T$  (moment-inversion estimator)")
    axR.set_title("Effective temperature vs step size (★ = Δt→0 extrapolation)")
    axR.legend(fontsize=8)
    axR.grid(alpha=0.3)
    axR.set_xlim(left=0)

    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def plot_stationarity(results, out_png, tag):
    stat = results.get("stationarity", {})
    if not stat:
        return
    colors = {"cayley_explicit": "#c0392b", "euler_spherical": "#2980b9",
              "cayley_true_midpoint": "#27ae60"}
    fig, ax = plt.subplots(figsize=(7, 5))
    for integ, d in stat.items():
        pl = d["plateau"]
        xx = pl.get("t_ns", pl["fracs"])
        ax.plot(xx, pl["teff"], marker="o", color=colors.get(integ, "k"),
                label=f"{integ} (τ_int={d['tau_int_ps']:.1f} ps)")
    ax.axhline(1.0, color="k", lw=1, alpha=0.5)
    ax.set_xlabel("trajectory length used (ns)")
    ax.set_ylabel(r"$T_{\rm eff}/T$ (moment inversion)")
    ax.set_title(f"Stationarity plateau ({tag})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────
# Configs
# ─────────────────────────────────────────────────────────────────────────
# τ_int ≈ 0.28 ns in this deep well (Δ≈48.5, α=0.05) → trajectories must span
# many τ_int and the eq-start length-plateau must be flat before T_eff/T is
# trusted. The PILOT is deliberately under-length (≈ 8 ns ≈ 28 τ_int): it
# validates the harness and shows the climbing trend, but its T_eff/T is a
# *lower bound* on the stationary value and is reported as PRELIMINARY only.
PILOT_CFG = dict(
    integrators=["cayley_explicit", "euler_spherical", "cayley_true_midpoint"],
    dt_list_ps=[0.05, 0.1, 0.2],
    T_phys_s=3.0e-9,        # 3 ns physical window (fixed across dt)
    burn_frac=0.10,         # eq-start: small burn-in only
    n_traj=12,
    n_boot=300,
    base_seed=20260630,
    mz_window=0.85,
    n_bins=25,
)

# FULL run: long enough (≈ 60 ns ≈ 210 τ_int) for the eq-start plateau to
# settle, 4 dt points for a clean Δt→0 Richardson fit, 48 trajectories for
# block-bootstrap σ(T_eff/T) ≲ 0.03. ≈ 9 h single-core (≈ overnight); the
# expensive dt=0.025 ps point is dropped (it alone doubles wall time and the
# linear-in-dt fit is well constrained by 0.05–0.4 ps). Re-add it / raise
# n_traj if more machine time is available.
FULL_CFG = dict(
    integrators=["cayley_explicit", "euler_spherical", "cayley_true_midpoint"],
    dt_list_ps=[0.05, 0.1, 0.2, 0.4],
    T_phys_s=60.0e-9,       # 60 ns physical window (fixed across dt)
    burn_frac=0.25,
    n_traj=48,
    n_boot=2000,
    base_seed=20260630,
    mz_window=0.85,
    n_bins=30,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true", help="run the full overnight sweep")
    ap.add_argument("--integrators", type=str, default=None,
                    help="comma-separated subset of integrators to run "
                         "(e.g. after a partial/killed run: "
                         "euler_spherical,cayley_true_midpoint)")
    ap.add_argument("--tag", type=str, default=None,
                    help="output-filename tag override (avoids clobbering a "
                         "previous partial run's outputs)")
    args = ap.parse_args()

    cfg = dict(FULL_CFG if args.full else PILOT_CFG)
    tag = args.tag or ("full" if args.full else "pilot")
    if args.integrators:
        cfg["integrators"] = args.integrators.split(",")

    out_dir = os.path.join("result", "sec_2_2_3_2", "B")
    os.makedirs(out_dir, exist_ok=True)

    print(f"=== Experiment B — Boltzmann T_eff ({tag}) ===", flush=True)
    print(f"config: {cfg}", flush=True)

    results = run_sweep(cfg, checkpoint_path=os.path.join(
        out_dir, f"teff_dtsweep_{tag}_checkpoint.json"))

    json_path = os.path.join(out_dir, f"teff_dtsweep_{tag}.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    png_path = os.path.join(out_dir, f"teff_dtsweep_{tag}.png")
    plot_results(results, png_path, tag)
    stat_png = os.path.join(out_dir, f"teff_stationarity_{tag}.png")
    plot_stationarity(results, stat_png, tag)

    print("\n=== SUMMARY ===", flush=True)
    print(f"Delta (analytic)          = {results['meta']['Delta']:.4f}")
    print(f"harmonic construction bias= {results['meta']['harmonic_construction_bias']:.4f} "
          f"(why we do NOT use 1/<1-m_z^2>)")
    for integ in cfg["integrators"]:
        ex = results["extrapolation"].get(integ, {}).get("moment", {})
        if ex:
            print(f"  {integ:>22}: dt->0 T_eff/T = {ex['intercept']:.4f} "
                  f"± {ex['intercept_sigma']:.4f} "
                  f"(dev {ex['intercept_dev_sigma']:+.1f}σ); "
                  f"slope = {ex['slope']:+.4f}/ps ({ex['slope_pos_sigma']:+.1f}σ)")
    print(f"VERDICT SC5_false = {results['verdict']['sc5_false']}")
    for r in results["verdict"]["reasons"]:
        print(f"   - {r}")
    print(f"wall time = {results['meta']['wall_time_s']:.1f} s")
    print(f"\nwrote: {json_path}\n       {png_path}\n       {stat_png}")


if __name__ == "__main__":
    main()
