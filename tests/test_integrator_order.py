"""Deterministic order-of-accuracy regression for §2.2.3.2 (Experiment A).

Guards the chapter's "core symptom" finding: the published Cayley stepper
``switching_vector`` evaluates omega ONLY at the left endpoint m_n (explicit-
omega), so it is GLOBALLY FIRST-ORDER (p~1) under a time-varying omega(m(t)) --
NOT the second-order implicit midpoint the text claims. A true iterated-midpoint
built from the SAME ``_omega_from_state`` / ``cayley_step`` primitives recovers
p~2, and the legacy spherical-Euler stepper is an independent first-order anchor.

Red-team hardening baked into this test:
  * the order-2 slope is fitted on the SMALL-dt ASYMPTOTIC subset only (the
    coarsest dt levels, where omega*dt ~ O(1), are dropped) -- otherwise the
    pre-asymptotic curvature biases the published slope upward;
  * a TRUE constant-omega control (alpha=0, K_i off, demag off, I_SOT=0, only a
    uniform field so omega is constant) proves the published step is EXACTLY 2nd
    order when omega is frozen, isolating the order loss to the TIME-VARIATION of
    omega(m(t)) rather than damping or the geometric Cayley map itself.

Acceptance gates (per the experiment SPEC):
    p_published in [0.85, 1.15]   and   p_midpoint in [1.85, 2.15].

Only public API is imported; no edits to src/. Kept fast (sub-second) by using a
modest dt ladder that is still safely inside the asymptotic regime.
"""
from __future__ import annotations

import numpy as np
import pytest

from vgsot_sim.anisotropy import field
from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.dynamic_switching import switching
from vgsot_sim.dynamic_switching_vector import _omega_from_state, cayley_step

# Deterministic super-threshold SOT pulse (NON=0). Same workhorse config as the
# committed script scripts/02_integrator/order_of_accuracy.py.
PARAMS = dict(
    V_MTJ=0.0, I_SOT=-1.5e-3, R_MTJ=5000.0, ESOT=1, ESTT=0,
    R_SOT_FL_DL=0.83, VNV=0, NON=0, ENE=1,
)
T_HORIZON = 0.2e-9
# Asymptotic ladder used by the test (smaller than the script's, for speed):
N_STEPS_LEVELS = [64, 128, 256, 512, 1024]
N_REF = 32768                 # >=32x finer than the finest fit level (1024)
N_DROP_COARSE = 1             # drop the coarsest dt for the asymptotic fit


def _normalize(v):
    return v / np.linalg.norm(v)


def _initial_state():
    return _normalize(np.array([np.sin(0.3), 0.0, np.cos(0.3)]))


def _omega_at(m, cc):
    m = np.asarray(m, float)
    theta = float(np.arccos(np.clip(m[2], -1.0, 1.0)))
    phi = float(np.arctan2(m[1], m[0]))
    H_eff, _ = field(theta, phi, PARAMS["V_MTJ"], n=1, NON=PARAMS["NON"],
                     ENE=PARAMS["ENE"], VNV=PARAMS["VNV"], constants=cc,
                     demag_mode="ellipsoid")
    H_eff = np.asarray(H_eff, float)
    Ms_use = cc.Ms
    I_MTJ = PARAMS["V_MTJ"] / PARAMS["R_MTJ"] if PARAMS["R_MTJ"] else 0.0
    J_STT = I_MTJ / cc.A1
    J_SOT = PARAMS["I_SOT"] / cc.A2
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
    H_DL_STT = PARAMS["ESTT"] * cc.h_bar * cc.P * J_STT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_DL_SOT = PARAMS["ESOT"] * cc.h_bar * cc.theta_SH * J_SOT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_FL_SOT = PARAMS["R_SOT_FL_DL"] * H_DL_SOT
    return _omega_from_state(
        m, H_eff, np.array([-1.0, 0.0, 0.0]),
        H_DL_SOT=H_DL_SOT, H_FL_SOT=H_FL_SOT, H_DL_STT=H_DL_STT, H_FL_STT=0.0,
        sigma_STT=np.array([0.0, 0.0, 1.0]),
        alpha=cc.alpha, gamma_red=gamma_red,
    )


def _run_published(m0, dt, n):
    cc = PhysicalConstantsConfig()
    cc.t_step = dt
    m = m0.copy()
    for _ in range(n):
        m = _normalize(cayley_step(m, _omega_at(m, cc), dt))
    return m


def _run_midpoint(m0, dt, n, niter=3):
    cc = PhysicalConstantsConfig()
    cc.t_step = dt
    m = m0.copy()
    for _ in range(n):
        mn = m.copy()
        for _ in range(niter):
            om = _omega_at(_normalize((m + mn) / 2.0), cc)
            mn = cayley_step(m, om, dt)
        m = _normalize(mn)
    return m


def _run_euler_spherical(m0, dt, n):
    cc = PhysicalConstantsConfig()
    cc.t_step = dt
    theta = float(np.arccos(np.clip(m0[2], -1.0, 1.0)))
    phi = float(np.arctan2(m0[1], m0[0]))
    for _ in range(n):
        _, phi, theta = switching(
            PARAMS["V_MTJ"], PARAMS["I_SOT"], PARAMS["R_MTJ"], theta, phi,
            ESTT=0, ESOT=1, VNV=0, NON=0,
            R_SOT_FL_DL=PARAMS["R_SOT_FL_DL"], constants=cc)
    return np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi), np.cos(theta)])


def _fit_slope(dts, errs):
    x, y = np.log(np.asarray(dts)), np.log(np.asarray(errs))
    p, b = np.polyfit(x, y, 1)
    yhat = p * x + b
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(p), float(r2)


def _curve(runner, m0, m_ref):
    dts, errs = [], []
    for n in N_STEPS_LEVELS:
        dt = T_HORIZON / n
        errs.append(float(np.linalg.norm(runner(m0, dt, n) - m_ref)))
        dts.append(dt)
    return np.array(dts), np.array(errs)


def _asym_slope(dts, errs):
    sl = slice(N_DROP_COARSE, None)
    return _fit_slope(dts[sl], errs[sl])


@pytest.fixture(scope="module")
def reference():
    m0 = _initial_state()
    return m0, _run_midpoint(m0, T_HORIZON / N_REF, N_REF)


def test_published_is_first_order(reference):
    """Published explicit-omega Cayley step: GLOBAL first order (p ~ 1)."""
    m0, m_ref = reference
    dts, errs = _curve(_run_published, m0, m_ref)
    p, r2 = _asym_slope(dts, errs)
    assert 0.85 <= p <= 1.15, f"p_published={p:.4f} (R2={r2:.4f}) not ~1"
    assert r2 > 0.99, f"published asymptotic fit R2={r2:.4f} too low"


def test_midpoint_is_second_order(reference):
    """True iterated-midpoint (same primitives): SECOND order (p ~ 2)."""
    m0, m_ref = reference
    dts, errs = _curve(_run_midpoint, m0, m_ref)
    p, r2 = _asym_slope(dts, errs)
    assert 1.85 <= p <= 2.15, f"p_midpoint={p:.4f} (R2={r2:.4f}) not ~2"
    assert r2 > 0.99, f"midpoint asymptotic fit R2={r2:.4f} too low"


def test_euler_spherical_anchor_is_first_order(reference):
    """Legacy spherical-Euler stepper: independent first-order anchor (p ~ 1)."""
    m0, m_ref = reference
    dts, errs = _curve(_run_euler_spherical, m0, m_ref)
    p, r2 = _asym_slope(dts, errs)
    assert 0.85 <= p <= 1.15, f"p_euler_spherical={p:.4f} (R2={r2:.4f}) not ~1"


def test_constant_omega_control_is_second_order():
    """Constant-omega masking control: published step is EXACTLY 2nd order.

    alpha=0, K_i off & demag off (omega built directly from a uniform field),
    I_SOT=0 -> omega is constant. The exact solution is a fixed-axis rotation
    (Rodrigues). p~2 here proves the published step's order loss in the full
    problem comes from the TIME-VARIATION of omega(m(t)), not from damping or
    the Cayley map.
    """
    cc = PhysicalConstantsConfig()
    cc.alpha = 0.0
    H = np.array([3000.0, 1000.0, 2000.0])
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)

    def omega_const(m):
        return _omega_from_state(
            m, H, np.array([-1.0, 0.0, 0.0]),
            H_DL_SOT=0.0, H_FL_SOT=0.0, H_DL_STT=0.0, H_FL_STT=0.0,
            sigma_STT=np.array([0.0, 0.0, 1.0]),
            alpha=cc.alpha, gamma_red=gamma_red)

    m0 = _normalize(np.array([np.sin(0.5), 0.2, np.cos(0.5)]))
    om0 = omega_const(m0)
    axis = om0 / np.linalg.norm(om0)
    ang = np.linalg.norm(om0) * T_HORIZON
    m_ref = (m0 * np.cos(ang) + np.cross(axis, m0) * np.sin(ang)
             + axis * np.dot(axis, m0) * (1.0 - np.cos(ang)))

    dts, errs = [], []
    for n in N_STEPS_LEVELS:
        dt = T_HORIZON / n
        m = m0.copy()
        for _ in range(n):
            m = _normalize(cayley_step(m, omega_const(m), dt))
        dts.append(dt)
        errs.append(float(np.linalg.norm(m - m_ref)))
    p, r2 = _fit_slope(dts, errs)
    assert 1.9 <= p <= 2.1, f"constant-omega p={p:.4f} (R2={r2:.4f}) not ~2"


def test_published_and_midpoint_separated(reference):
    """The two slopes must be cleanly separated (no overlap of CI bands)."""
    m0, m_ref = reference
    p_pub, _ = _asym_slope(*_curve(_run_published, m0, m_ref))
    p_mid, _ = _asym_slope(*_curve(_run_midpoint, m0, m_ref))
    assert p_mid - p_pub > 0.7, (
        f"published (p={p_pub:.3f}) and midpoint (p={p_mid:.3f}) not separated")
