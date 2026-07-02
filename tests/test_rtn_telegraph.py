"""Tests for the RTN reservoir primitive (vgsot_sim.rtn.telegraph)."""
from __future__ import annotations

import numpy as np
import pytest

from vgsot_sim.rtn import (
    TelegraphArray,
    TelegraphParams,
    neel_brown_rate,
    relaxation_time,
    simulate_trace,
    stationary_mean,
    tau_max,
    up_down_rates,
)
from vgsot_sim.analysis.nb_fit import NBFitResult, psw_nb


def test_stationary_mean_equals_tanh_and_rate_ratio():
    V = np.linspace(-0.3, 0.3, 25)
    Delta, Vc0 = 5.15, 0.884
    s_tanh = stationary_mean(V, Delta=Delta, Vc0=Vc0)
    assert np.allclose(s_tanh, np.tanh(Delta * V / Vc0))
    r_up, r_dn = up_down_rates(V, Delta=Delta, Vc0=Vc0)
    s_rates = (r_up - r_dn) / (r_up + r_dn)
    assert np.allclose(s_tanh, s_rates, atol=1e-12)


def test_tau_peaks_at_zero_bias():
    V = np.linspace(-0.3, 0.3, 51)
    tau = relaxation_time(V, tau0=1.0, Delta=5.15, Vc0=0.884)
    i0 = np.argmin(np.abs(V))
    assert i0 == np.argmax(tau)                                   # peak at V=0
    assert np.isclose(tau[i0], tau_max(tau0=1.0, Delta=5.15))     # tau0*exp(Delta)/2
    # monotone decay away from zero
    assert np.all(np.diff(tau[i0:]) < 0) and np.all(np.diff(tau[:i0 + 1]) > 0)


def test_relaxation_time_is_inverse_total_rate():
    V = np.linspace(-0.2, 0.2, 11)
    r_up, r_dn = up_down_rates(V, tau0=1.0, Delta=5.15, Vc0=0.884)
    assert np.allclose(relaxation_time(V, tau0=1.0, Delta=5.15, Vc0=0.884),
                       1.0 / (r_up + r_dn))


def test_exact_propagator_recovers_stationary_mean():
    """Ensemble mean after many steps must match tanh(Delta V/Vc0)."""
    p = TelegraphParams(tau0=1.0, Delta=5.15, Vc0=0.884)
    n = 20000
    for V in (0.0, 0.1, -0.15):
        arr = TelegraphArray(n, p, seed=123)
        for _ in range(40):                       # 40 * 20 ns >> tau -> stationary
            arr.step(V, dt=20.0)
        m = arr.state.mean()
        assert abs(m - stationary_mean(V, Delta=p.Delta, Vc0=p.Vc0)) < 0.02, (V, m)


def test_states_are_binary_and_seed_reproducible():
    p = TelegraphParams()
    a = TelegraphArray(500, p, seed=7)
    b = TelegraphArray(500, p, seed=7)
    for _ in range(10):
        sa = a.step(0.05, dt=10.0)
        sb = b.step(0.05, dt=10.0)
        assert set(np.unique(sa)).issubset({-1.0, 1.0})
        assert np.array_equal(sa, sb)             # same seed -> identical stream


def test_from_nb_fit_copies_delta_vc0():
    fit = NBFitResult(a=0.9, b=0.1, Delta=3.8, Vc0=0.86, tau_ret_ns=44.7)
    p = TelegraphParams.from_nb_fit(fit, tau0=1.0)
    assert p.Delta == 3.8 and p.Vc0 == 0.86 and p.tau0 == 1.0


def test_simulate_trace_shape_and_values():
    V_t = np.concatenate([np.full(50, 0.2), np.full(50, -0.2)])
    out = simulate_trace(V_t, dt=5.0, n=8, seed=1)
    assert out.shape == (100, 8)
    assert set(np.unique(out)).issubset({-1.0, 1.0})


def test_rate_clipped_at_attempt_frequency():
    """Rate never exceeds the attempt frequency 1/tau0 (fixes the |V|>Vc0 blow-up)."""
    tau0, Delta, Vc0 = 1.0, 5.15, 0.884
    V = np.linspace(-1.5, 1.5, 121)
    r = neel_brown_rate(V, tau0=tau0, Delta=Delta, Vc0=Vc0)
    assert np.all(r <= 1.0 / tau0 + 1e-12)                       # capped at attempt freq
    # inside the physical domain |V|<Vc0 the clip is inactive -> matches the raw law
    inside = np.abs(V) < Vc0
    raw = (1.0 / tau0) * np.exp(-Delta * (1.0 - V[inside] / Vc0))
    assert np.allclose(r[inside], raw)
    # exactly the floored exponent everywhere
    expected = (1.0 / tau0) * np.exp(-np.maximum(Delta * (1.0 - V / Vc0), 0.0))
    assert np.allclose(r, expected)


def test_rate_matches_psw_nb_instantaneous():
    """neel_brown_rate equals the instantaneous rate implied by nb_fit.psw_nb."""
    tau0, Delta, Vc0 = 1.0, 5.15, 0.884
    V = np.linspace(-1.2, 1.2, 49)
    t = 1e-4                                                     # any t: relation is exact
    implied = -np.log1p(-psw_nb(V, t, Delta, Vc0, tau0=tau0)) / t
    assert np.allclose(implied, neel_brown_rate(V, tau0=tau0, Delta=Delta, Vc0=Vc0))


def test_step_warns_outside_domain():
    """step() flags |V|>Vc0 as outside the two-state model's physical domain."""
    arr = TelegraphArray(16, TelegraphParams(Vc0=0.884), seed=0)
    with pytest.warns(RuntimeWarning):
        arr.step(1.2, dt=1.0)
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("error")                                # no warning inside domain
        arr.step(0.3, dt=1.0)
