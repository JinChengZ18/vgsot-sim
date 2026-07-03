"""Fast unit tests for the RTN reservoir layer (vgsot_sim.rtn.reservoir).

Small/short configs so the suite stays quick; the full benchmark lives in
scripts/10_rtn_reservoir/benchmark_reservoir.py.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

import pytest

from vgsot_sim.rtn.reservoir import (
    DelayReservoir,
    Reservoir,
    ReservoirConfig,
    information_processing_capacity,
    memory_capacity,
    narma10,
    parity_task,
    ridge_fit,
    ridge_predict,
    ring_reservoir,
)


def test_ridge_recovers_linear_map():
    rng = np.random.default_rng(0)
    X = rng.standard_normal((400, 5))
    W_true = rng.standard_normal((5, 1))
    Y = X @ W_true + 0.7
    W = ridge_fit(X, Y, alpha=1e-10)
    pred = ridge_predict(X, W)
    assert np.allclose(pred, Y, atol=1e-4)


def test_reservoir_run_shape_and_range():
    r = Reservoir(ReservoirConfig(n_nodes=20), seed=0)
    u = np.random.default_rng(0).uniform(-1, 1, size=300)
    X = r.run(u, washout=50)
    assert X.shape == (250, 20)
    assert np.all(np.abs(X) <= 1.0 + 1e-9)             # node states are bounded


def test_heterogeneous_beats_broadcast_baseline():
    """W_in + heterogeneity must give materially more memory than identical broadcast (A2)."""
    cfg = ReservoirConfig(n_nodes=40)
    het = Reservoir(cfg, seed=2)
    mc_het = memory_capacity(het, n_samples=900, max_delay=20, washout=100, seed=2)

    base = Reservoir(replace(cfg, delta_range=(2.0, 2.0), bias_spread=0.0), seed=2)
    base.W_in[:] = 1.0                                  # identical nodes, same input
    mc_base = memory_capacity(base, n_samples=900, max_delay=20, washout=100, seed=2)

    assert mc_base.mc_total < 1.5
    assert mc_het.mc_total > 3.0
    assert mc_het.mc_total > 3.0 * mc_base.mc_total


def test_memory_capacity_decays_with_delay():
    r = Reservoir(ReservoirConfig(n_nodes=40), seed=3)
    mc = memory_capacity(r, n_samples=900, max_delay=15, washout=100, seed=3)
    assert mc.mc_k[0] > mc.mc_k[-1]                     # immediate recall > long delay
    assert mc.mc_k[0] > 0.5                             # recovers u[t-1] well


def test_narma10_series_properties():
    u = np.random.default_rng(0).uniform(0.0, 0.5, size=200)
    y = narma10(u)
    assert y.shape == u.shape
    assert np.all(y[:10] == 0.0)                        # warmup
    assert np.all(np.isfinite(y))


def test_stochastic_mode_runs():
    r = Reservoir(ReservoirConfig(n_nodes=16), seed=4)
    u = np.random.default_rng(0).uniform(-1, 1, size=120)
    X = r.run(u, mode="stochastic", washout=20, n_replicas=2, seed=4)
    assert X.shape == (100, 16)
    assert np.all(np.abs(X) <= 1.0 + 1e-9)


def test_parity_task_solvable_nonlinearly():
    """The reservoir's nonlinear mixing solves parity-2 (a linear readout alone cannot)."""
    r = Reservoir(ReservoirConfig(n_nodes=60), seed=2)
    acc = parity_task(r, n_bits=2, n_samples=1500, washout=100, seed=2)
    assert acc > 0.9


def test_ipc_has_linear_and_nonlinear_capacity():
    r = Reservoir(ReservoirConfig(n_nodes=60), seed=3)
    ipc = information_processing_capacity(r, max_delay=12, pair_delay=6,
                                          n_samples=2000, seed=3)
    assert ipc.deg1 > 1.0 and ipc.deg2 > 0.2          # both memory and nonlinearity
    assert np.isclose(ipc.total, ipc.deg1 + ipc.deg2)


def test_delay_reservoir_runs_and_has_memory():
    dr = DelayReservoir(n_virtual=30, fb_gain=0.3, seed=7)
    u = np.random.default_rng(0).uniform(-1, 1, size=200)
    X = dr.run(u, washout=50)
    assert X.shape == (150, 30)
    assert np.all(np.abs(X) <= 1.0 + 1e-9)
    mc = memory_capacity(dr, max_delay=20, n_samples=900, seed=7)
    assert mc.mc_total > 1.5                            # feedback gives nontrivial memory


def test_ring_reservoir_beats_filter_bank_mc():
    """The tuned simple-cycle delay line must clearly exceed the filter-bank MC."""
    ring = ring_reservoir(40, seed=3)
    mc_ring = memory_capacity(ring, max_delay=50, n_samples=1500, washout=150, seed=3)
    fb = Reservoir(ReservoirConfig(n_nodes=40), seed=3)
    mc_fb = memory_capacity(fb, max_delay=50, n_samples=1500, washout=150, seed=3)
    assert mc_ring.mc_total > 1.5 * mc_fb.mc_total


def test_coupling_validation_and_stochastic_guard():
    with pytest.raises(ValueError):
        Reservoir(ReservoirConfig(n_nodes=8, coupling_radius=0.5,
                                  coupling_topology="hexagon"), seed=0)
    with pytest.raises(ValueError):
        Reservoir(ReservoirConfig(n_nodes=8, input_mode="everywhere"), seed=0)
    r = ring_reservoir(8, seed=0)
    u = np.random.default_rng(0).uniform(-1, 1, size=50)
    with pytest.raises(NotImplementedError):
        r.run(u, mode="stochastic", seed=0)


def test_random_coupling_runs_and_is_bounded():
    r = Reservoir(ReservoirConfig(n_nodes=16, coupling_radius=0.6,
                                  coupling_topology="random"), seed=1)
    u = np.random.default_rng(1).uniform(-1, 1, size=200)
    X = r.run(u, washout=40)
    assert X.shape == (160, 16)
    assert np.all(np.abs(X) <= 1.0 + 1e-9)


def test_mackey_glass_series_and_task():
    """Smoke test only: the MG readout Gram is ill-conditioned (cond 1e14-1e20),
    so exact NRMSE values drift across BLAS invocations — assert capability, not
    a benchmark figure (see scripts/10_rtn_reservoir/README.md)."""
    from vgsot_sim.rtn.reservoir import mackey_glass, mackey_glass_task
    x = mackey_glass(1500, seed=0)
    assert x.shape == (1500,)
    assert 0.2 < x.min() and x.max() < 1.6          # bounded chaotic attractor
    assert x.std() > 0.1                             # not collapsed to fixed point
    r = Reservoir(ReservoirConfig(n_nodes=40), seed=4)
    t = mackey_glass_task(r, horizon=10, n_samples=1200, washout=150,
                          alpha=1e-3, seed=4)
    assert np.isfinite(t.nrmse) and t.nrmse < 1.2    # short-horizon is learnable


def test_kernel_quality_and_esp():
    from vgsot_sim.rtn.reservoir import esp_convergence, kernel_quality
    fb = Reservoir(ReservoirConfig(n_nodes=24), seed=3)
    kr = kernel_quality(fb, n_streams=16, t_len=40, washout=80, seed=1)
    gr = kernel_quality(fb, n_streams=16, t_len=40, washout=80,
                        generalization=True, seed=1)
    assert 1 <= kr <= 16 and 1 <= gr <= 16
    assert gr <= kr                                   # shared tails collapse rank
    assert esp_convergence(fb, t_len=300, seed=2) < 1e-8   # leaky nodes forget init
