"""Fast unit tests for the sLLG<->RTN bridge (vgsot_sim.rtn.bridge).

These exercise the barrier/bias algebra and the analysis helpers on synthetic
traces. The expensive LLG validation (long free-running runs) lives in
scripts/10_rtn_reservoir/validate_bridge.py, not here.
"""
from __future__ import annotations

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.rtn import bridge


def test_ki_delta_roundtrip():
    c = PhysicalConstantsConfig()
    for d in (1.5, 2.0, 5.0, 20.0):
        ki = bridge.ki_for_delta(d, c)
        assert np.isclose(bridge.delta_of_ki(ki, c), d, rtol=1e-9)


def test_low_barrier_is_perpendicular_not_in_plane():
    """ki_for_delta keeps K_u^eff positive (perpendicular well), unlike a PMA-only Ki."""
    c = PhysicalConstantsConfig()
    from vgsot_sim.material_temperature import k_u_eff_of_T
    c_low = bridge.low_barrier_constants(2.0, c)
    assert k_u_eff_of_T(c_low.T, c_low) > 0.0            # still perpendicular
    assert c_low.h_ex_x == 0.0 and c_low.h_ex_y == 0.0    # in-plane bias zeroed
    # a naive PMA-only Ki for Delta=2 (ignoring demag) would go in-plane:
    import dataclasses
    ki_pma = 2.0 * c.kb * c.T / c.A1                      # Ki*A1/(kBT)=Delta
    c_bad = dataclasses.replace(c, Ki=ki_pma)
    assert k_u_eff_of_T(c_bad.T, c_bad) < 0.0             # demag flips it in-plane


def test_tilt_per_field_value():
    c = PhysicalConstantsConfig()
    assert np.isclose(bridge.tilt_per_field(c), c.u0 * c.Ms * c.v / (c.kb * c.T))


def test_free_run_shape_and_range():
    c = bridge.low_barrier_constants(2.0, PhysicalConstantsConfig())
    mz = bridge.free_run(c, 200, seed=0)                  # tiny: just shape/range
    assert mz.shape == (200,)
    assert np.all(np.abs(mz) <= 1.0 + 1e-9)


def test_dwell_times_on_square_wave():
    """A clean ±1 square wave has dwell = half-period, CV→0."""
    dt = 1e-9
    period = 100                                          # samples
    s = np.tile(np.r_[np.ones(period // 2), -np.ones(period // 2)], 20)
    dw = bridge.dwell_times(s, dt, thr=0.5)
    assert dw.n_flips > 10
    assert np.isclose(dw.tau_mean_ns, (period // 2) * dt * 1e9, rtol=0.05)
    assert dw.cv < 0.05                                   # perfectly regular


def test_dwell_times_too_few_flips_returns_nan():
    dw = bridge.dwell_times(np.ones(1000), 1e-9, thr=0.5)
    assert dw.n_flips == 0 and np.isnan(dw.tau_mean_ns)


def test_psd_lorentzian_runs_on_synthetic():
    rng = np.random.default_rng(0)
    x = np.cumsum(rng.standard_normal(8192))             # red-ish noise
    x = np.tanh(0.01 * x)
    psd = bridge.psd_lorentzian(x, 1e-9, n_seg=4)
    assert psd.P.ndim == 1 and np.all(psd.P >= 0.0)
