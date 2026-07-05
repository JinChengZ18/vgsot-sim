"""Discriminate the mechanism of the euler_spherical T_eff = T/2 defect.

Three variants of the spherical-Euler stepper, identical injected noise
stream structure, deep uniaxial well (h_ex zeroed, ENE effect nulled by
config), analytic-equilibrium starts:

  (a) baseline           — shipped switching()               -> expect T_eff/T ~ 0.50
  (b) + Wong-Zakai drift — theta += D_theta*cot(theta)*dt    -> expect ~ 1.00 (mechanism = missing
                            noise-induced drift of the (theta,phi) chart)
  (c) amp * sqrt(2)      — doubled injected variance          -> moment ~ 1.0 BUT slope ~ 1.14
                            (an amplitude error CANNOT repair both estimators)

The two estimators (bounded-basin moment inversion; WLS histogram slope of
ln p vs 1-m_z^2) respond differently under (b) vs (c), so the pair
(moment, slope) identifies the mechanism uniquely:
  missing-drift theory: (a) (0.50, ~0.57)  (b) (1.0, 1.0)  (c) (1.0, ~1.14)
  amplitude theory:     (a) (0.50, 0.50)   (b) n/a          (c) (1.0, 1.0)
Measured full-run (a) = (0.488-0.507, 0.559-0.575) already favours
missing-drift; (b)+(c) close the case.

Run:  PYTHONPATH=src python scripts/11_boltzmann_teff/euler_drift_discrimination.py
"""
from __future__ import annotations

import json
import math
import os
import time

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.demag import demag_factors
from vgsot_sim.dynamic_switching import switching

OUT = os.path.join("result", "sec_2_2_3_2", "B")
os.makedirs(OUT, exist_ok=True)

DT_PS = 0.2
T_PHYS = 12.0e-9
BURN = 0.25
N_TRAJ = 16
T_K = 300.0
SEED0 = 20260703

cc = PhysicalConstantsConfig()
cc.h_ex_x = cc.h_ex_y = cc.h_ex_z = 0.0          # pure uniaxial well
dt = DT_PS * 1e-12
cc.t_step = dt

Nx, Ny, Nz = demag_factors(cc, mode="ellipsoid")
Keff = cc.Ki / cc.tf - 0.5 * cc.u0 * cc.Ms ** 2 * (Nz - Nx)
v = cc.tf * cc.A1
Delta = Keff * v / (cc.kb * T_K)
amp0 = math.sqrt(2.0 * cc.kb * T_K * cc.alpha
                 / (cc.u0 * cc.Ms * cc.gamma * v * dt))     # Brown/FDT per-component sigma
gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
# theta-noise diffusion const: Var(dtheta) = gamma_red^2 (1+alpha^2) sigma_H^2 dt^2 = 2 D dt
D_theta = 0.5 * gamma_red ** 2 * (1.0 + cc.alpha ** 2) * amp0 ** 2 * dt
print(f"Delta={Delta:.3f}  amp0={amp0:.3f} A/m  D_theta={D_theta:.4e} rad^2/s", flush=True)


def mean_u_of(delta):
    zn = quad(lambda m: math.exp(-delta * (1 - m * m)), 0, 1)[0]
    zu = quad(lambda m: (1 - m * m) * math.exp(-delta * (1 - m * m)), 0, 1)[0]
    return zu / zn


MEAN_U_TRUE = mean_u_of(Delta)


def teff_moment(u_obs):
    f = lambda d: mean_u_of(d) - u_obs
    d_eff = brentq(f, 1.0, 5000.0)
    return Delta / d_eff


def teff_slope(mz_eq, mz_window=0.85, n_bins=30):
    m = np.abs(mz_eq)
    m = m[m >= mz_window]
    hist, edges = np.histogram(m, bins=n_bins, density=True)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    ok = hist > 0
    u = 1.0 - ctr[ok] ** 2
    y = np.log(hist[ok])
    w = hist[ok]
    A = np.vstack([u, np.ones_like(u)]).T
    W = np.diag(w)
    beta = np.linalg.solve(A.T @ W @ A, A.T @ W @ y)
    return Delta / (-beta[0])


def eq_sample(rng):
    mzg = np.linspace(0.0, 1.0, 200001)
    w = np.exp(-Delta * (1.0 - mzg ** 2)); w /= w.sum()
    cdf = np.cumsum(w)
    mz = mzg[np.searchsorted(cdf, rng.random())]
    sign = 1.0 if rng.random() < 0.5 else -1.0
    return math.acos(max(min(sign * mz, 1.0), -1.0)), rng.uniform(0, 2 * math.pi)


def run_variant(name, add_drift, amp_scale):
    n_steps = int(round(T_PHYS / dt))
    n_burn = int(round(BURN * n_steps))
    amp = amp0 * amp_scale
    us, mzs = [], []
    for k in range(N_TRAJ):
        rng = np.random.default_rng(SEED0 + 977 * k)
        theta, phi = eq_sample(rng)
        mz_ts = np.empty(n_steps)
        for i in range(n_steps):
            h_th = amp * rng.normal(0.0, 1.0, 3)
            mz_i, phi, theta = switching(
                V_MTJ=0.0, I_SOT=0.0, R_MTJ=1.0, theta=theta, phi=phi,
                ESTT=0, ESOT=0, VNV=0, NON=1, constants=cc, h_th_ext=h_th)
            if add_drift:
                st = math.sin(theta)
                if abs(st) > 1e-8:
                    theta += D_theta * (math.cos(theta) / st) * dt
                    mz_i = math.cos(theta)
            mz_ts[i] = mz_i
        eq = mz_ts[n_burn:]
        us.append(1.0 - eq ** 2)
        mzs.append(eq)
    u_all = np.concatenate(us)
    mz_all = np.concatenate(mzs)
    tm = teff_moment(float(u_all.mean()))
    ts = teff_slope(mz_all)
    print(f"[{name:>22}] T_eff/T moment={tm:.4f}  slope={ts:.4f}", flush=True)
    return dict(moment=tm, slope=ts)


t0 = time.time()
res = {"Delta": Delta, "dt_ps": DT_PS, "T_phys_ns": T_PHYS * 1e9, "n_traj": N_TRAJ,
       "variants": {}}
res["variants"]["baseline"] = run_variant("baseline euler", False, 1.0)
res["variants"]["wong_zakai_drift"] = run_variant("euler + D*cot(theta)", True, 1.0)
res["variants"]["amp_sqrt2"] = run_variant("euler, amp*sqrt(2)", False, math.sqrt(2.0))
res["runtime_s"] = time.time() - t0

b, d, a = (res["variants"][k] for k in ("baseline", "wong_zakai_drift", "amp_sqrt2"))
drift_wins = (abs(d["moment"] - 1) < 0.12 and abs(d["slope"] - 1) < 0.15
              and a["slope"] > 1.05)
res["verdict"] = ("missing_wong_zakai_drift" if drift_wins else "inconclusive_see_numbers")
print(f"\nVERDICT: {res['verdict']}  (runtime {res['runtime_s']:.0f}s)", flush=True)
print("  theory grid — drift-fix should give (~1,~1); amp*sqrt2 should give (~1,>1.1)", flush=True)
with open(os.path.join(OUT, "euler_drift_discrimination.json"), "w") as f:
    json.dump(res, f, indent=2)
