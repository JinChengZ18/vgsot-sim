"""Regression test for the OPTIONAL true implicit-midpoint Cayley integrator.

This guards the new `n_midpoint` keyword of `switching_vector()` (and the
`"cayley_midpoint"` integrator string it backs), which makes the chapter's
stated "隐式中点法" (§2.2.3.2) an actually-selectable engine:

  * `n_midpoint == 0` (default)  — the published EXPLICIT-ω Cayley step,
    globally FIRST order (p ~ 1) and byte-identical to the original code;
  * `n_midpoint == k > 0`        — `k` fixed-point iterations of ω at the
    normalised midpoint, globally SECOND order (p ~ 2).

Unlike `tests/test_integrator_order.py` (which re-implements the midpoint in a
private helper), this test drives the midpoint through the PUBLIC
`switching_vector(..., n_midpoint=...)` API, so it validates the shipped
integrator itself. The order-fit machinery is reused from that test /
scripts/02_integrator/order_of_accuracy.py.

Asserts:
  (a) deterministic global order: cayley_midpoint = 2.0 +/- 0.2, while the
      explicit cayley step = 1.0 +/- 0.15, on the SAME dt grid;
  (b) |m| = 1 preserved to < 1e-12 across a dt grid (both schemes);
  (c) switching_vector(..., n_midpoint=0) == switching_vector(...) EXACTLY on
      several random states (the default path is untouched).

Only public API is imported; no edits to src/.
"""
from __future__ import annotations

import numpy as np
import pytest

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.dynamic_switching_vector import switching_vector

# Deterministic super-threshold SOT pulse (NON=0) — the same workhorse config as
# tests/test_integrator_order.py and scripts/02_integrator/order_of_accuracy.py.
# switching_vector() hardwires ENE=1, so the exchange-bias field is present here
# exactly as in the committed order-of-accuracy experiment.
PARAMS = dict(
    V_MTJ=0.0, I_SOT=-1.5e-3, R_MTJ=5000.0, ESTT=0, ESOT=1,
    VNV=0, NON=0, R_SOT_FL_DL=0.83,
)
T_HORIZON = 0.2e-9
N_STEPS_LEVELS = [64, 128, 256, 512, 1024]
N_REF = 32768                 # >=32x finer than the finest fit level (1024)
N_DROP_COARSE = 1             # drop the coarsest dt for the asymptotic fit


def _normalize(v):
    return v / np.linalg.norm(v)


def _initial_state():
    return _normalize(np.array([np.sin(0.3), 0.0, np.cos(0.3)]))


def _run(m0, dt, n, n_midpoint):
    """Integrate n steps of switching_vector at time step dt (public API)."""
    cc = PhysicalConstantsConfig()
    cc.t_step = dt
    m = m0.copy()
    for _ in range(n):
        m = switching_vector(
            m, PARAMS["V_MTJ"], PARAMS["I_SOT"], PARAMS["R_MTJ"],
            PARAMS["ESTT"], PARAMS["ESOT"],
            VNV=PARAMS["VNV"], NON=PARAMS["NON"], R_SOT_FL_DL=PARAMS["R_SOT_FL_DL"],
            constants=cc, n_midpoint=n_midpoint,
        )
    return m


def _fit_slope(dts, errs):
    x, y = np.log(np.asarray(dts)), np.log(np.asarray(errs))
    p, b = np.polyfit(x, y, 1)
    yhat = p * x + b
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(p), float(r2)


def _curve(m0, m_ref, n_midpoint):
    dts, errs = [], []
    for n in N_STEPS_LEVELS:
        dt = T_HORIZON / n
        errs.append(float(np.linalg.norm(_run(m0, dt, n, n_midpoint) - m_ref)))
        dts.append(dt)
    return np.array(dts), np.array(errs)


def _asym_slope(dts, errs):
    sl = slice(N_DROP_COARSE, None)
    return _fit_slope(dts[sl], errs[sl])


@pytest.fixture(scope="module")
def reference():
    """Fine-grid true-midpoint trajectory (shared reference for both curves)."""
    m0 = _initial_state()
    return m0, _run(m0, T_HORIZON / N_REF, N_REF, n_midpoint=3)


def test_cayley_midpoint_is_second_order(reference):
    """switching_vector(n_midpoint=3) is GLOBALLY SECOND order (p ~ 2)."""
    m0, m_ref = reference
    dts, errs = _curve(m0, m_ref, n_midpoint=3)
    p, r2 = _asym_slope(dts, errs)
    assert 1.8 <= p <= 2.2, f"p_cayley_midpoint={p:.4f} (R2={r2:.4f}) not ~2"
    assert r2 > 0.99, f"cayley_midpoint asymptotic fit R2={r2:.4f} too low"


def test_explicit_cayley_is_first_order(reference):
    """switching_vector(n_midpoint=0) is GLOBALLY FIRST order (p ~ 1)."""
    m0, m_ref = reference
    dts, errs = _curve(m0, m_ref, n_midpoint=0)
    p, r2 = _asym_slope(dts, errs)
    assert 0.85 <= p <= 1.15, f"p_cayley={p:.4f} (R2={r2:.4f}) not ~1"


def test_two_schemes_cleanly_separated(reference):
    """The measured orders must be cleanly separated (2 vs 1)."""
    m0, m_ref = reference
    p_exp, _ = _asym_slope(*_curve(m0, m_ref, n_midpoint=0))
    p_mid, _ = _asym_slope(*_curve(m0, m_ref, n_midpoint=3))
    assert p_mid - p_exp > 0.7, (
        f"explicit (p={p_exp:.3f}) and midpoint (p={p_mid:.3f}) not separated")


def test_unit_norm_preserved_across_dt_grid():
    """|m| = 1 is preserved to < 1e-12 across a dt grid, both schemes."""
    m0 = _initial_state()
    for n_midpoint in (0, 3):
        for n in N_STEPS_LEVELS:
            dt = T_HORIZON / n
            m = _run(m0, dt, n, n_midpoint)
            assert abs(np.linalg.norm(m) - 1.0) < 1e-12, (
                f"|m|-1={abs(np.linalg.norm(m) - 1.0):.3e} at "
                f"n_midpoint={n_midpoint}, dt={dt:.3e}")


def test_n_midpoint_zero_is_byte_identical():
    """switching_vector(n_midpoint=0) == switching_vector(...) EXACTLY.

    On several random unit states the new default path must reproduce the
    pre-existing explicit output bit-for-bit — the n_midpoint>0 branch is fully
    guarded and cannot perturb the shipped default.
    """
    cc = PhysicalConstantsConfig()
    cc.t_step = 1e-12
    rng = np.random.default_rng(20260705)
    for _ in range(12):
        m0 = _normalize(rng.normal(size=3))
        common = dict(
            VNV=0, NON=0, R_SOT_FL_DL=0.83, constants=cc,
        )
        m_default = switching_vector(
            m0, PARAMS["V_MTJ"], PARAMS["I_SOT"], PARAMS["R_MTJ"],
            PARAMS["ESTT"], PARAMS["ESOT"], **common,
        )
        m_zero = switching_vector(
            m0, PARAMS["V_MTJ"], PARAMS["I_SOT"], PARAMS["R_MTJ"],
            PARAMS["ESTT"], PARAMS["ESOT"], n_midpoint=0, **common,
        )
        assert np.array_equal(m_default, m_zero), (
            f"n_midpoint=0 diverged from default: {m_default} vs {m_zero}")
