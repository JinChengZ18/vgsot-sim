"""Tests for the RTN reservoir primitive (vgsot_sim.rtn.telegraph)."""
from __future__ import annotations

import numpy as np

from vgsot_sim.rtn import (
    TelegraphArray,
    TelegraphParams,
    relaxation_time,
    simulate_trace,
    stationary_mean,
    tau_max,
    up_down_rates,
)
from vgsot_sim.analysis.nb_fit import NBFitResult


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
