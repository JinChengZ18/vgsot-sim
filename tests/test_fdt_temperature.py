"""V4: the FDT thermal-field amplitude must track the lattice temperature T.

Brown 1963 / FDT: |H_th| ∝ sqrt(T). When self-heating is on, the noise variance
should rise with the instantaneous T(t). This pins field()'s new `T` argument.
"""
from __future__ import annotations

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.anisotropy import field


_KW = dict(theta=0.5, phi=0.3, V_MTJ=0.0, n=1, ENE=1, VNV=0)


def test_fdt_amplitude_scales_with_sqrt_T():
    cc = PhysicalConstantsConfig()
    # Deterministic part (NON=0): identical regardless of T / noise stream.
    H_det, _ = field(NON=0, constants=cc, **_KW)
    # Same seeded noise draw at two temperatures -> only the amplitude differs.
    H1, _ = field(NON=1, constants=cc, T=300.0, rng=np.random.default_rng(0), **_KW)
    H2, _ = field(NON=1, constants=cc, T=1200.0, rng=np.random.default_rng(0), **_KW)
    noise1 = np.asarray(H1) - np.asarray(H_det)
    noise2 = np.asarray(H2) - np.asarray(H_det)
    # |H_th| ∝ sqrt(T): ratio = sqrt(1200/300) = 2 on every component.
    assert np.allclose(noise2, 2.0 * noise1, rtol=1e-9, atol=0.0)
    assert np.linalg.norm(noise1) > 0.0           # noise actually present


def test_fdt_default_T_is_constants_T():
    cc = PhysicalConstantsConfig()
    H_none, _ = field(NON=1, constants=cc, T=None, rng=np.random.default_rng(1), **_KW)
    H_ct, _ = field(NON=1, constants=cc, T=cc.T, rng=np.random.default_rng(1), **_KW)
    assert np.allclose(np.asarray(H_none), np.asarray(H_ct), rtol=0, atol=0)


def test_h_th_ext_injects_external_field():
    """h_th_ext is added directly as H_TH (the harness-owned-noise injection path
    that matches vgsot_llg.va's hx/hy/hz nodes)."""
    cc = PhysicalConstantsConfig()
    v = np.array([1.2e4, -3.0e3, 5.0e3])
    H_inj, _ = field(NON=1, constants=cc, h_th_ext=v, **_KW)
    H_det, _ = field(NON=0, constants=cc, **_KW)
    assert np.allclose(np.asarray(H_inj) - np.asarray(H_det), v, rtol=0, atol=1e-6)


def test_switching_vector_h_th_ext_is_deterministic():
    """With an injected field (no internal draw) the step is RNG-free/reproducible."""
    from vgsot_sim.dynamic_switching_vector import switching_vector
    from vgsot_sim.initialize import compute_Rp
    cc = PhysicalConstantsConfig()
    m0 = np.array([0.1, 0.0, 0.995]); m0 = m0 / np.linalg.norm(m0)
    v = np.array([1.0e4, 0.0, 0.0])
    Rp = compute_Rp(cc)
    kw = dict(VNV=0, NON=1, R_SOT_FL_DL=0.83, constants=cc, h_th_ext=v)
    a = switching_vector(m0, 0.0, -1.0e-3, Rp, 0, 1, **kw)
    b = switching_vector(m0, 0.0, -1.0e-3, Rp, 0, 1, **kw)
    assert np.array_equal(a, b)
