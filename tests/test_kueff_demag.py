"""V3: k_u_eff_of_T default demag must be the exact oblate-ellipsoid tensor
(consistent with the field that drives the LLG), not the thin-disk approximation.
"""
from __future__ import annotations

from math import pi

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.material_temperature import k_u_eff_of_T
from vgsot_sim.demag import demag_factors


def test_kueff_default_uses_ellipsoid_demag():
    cc = PhysicalConstantsConfig()
    Nx, _, Nz = demag_factors(cc, mode="ellipsoid")
    explicit = k_u_eff_of_T(350.0, cc, Nz_minus_Nx=(Nz - Nx))
    default = k_u_eff_of_T(350.0, cc)
    assert np.isclose(default, explicit, rtol=1e-12)


def test_kueff_default_differs_from_thin_disk():
    cc = PhysicalConstantsConfig()
    Nx_td = pi * cc.tf / (4 * cc.D_elec)          # old thin-disk approximation
    thin = k_u_eff_of_T(350.0, cc, Nz_minus_Nx=(1 - 3 * Nx_td))
    default = k_u_eff_of_T(350.0, cc)
    assert not np.isclose(default, thin, rtol=1e-9)
