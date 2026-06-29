"""V2: R_P must be the honest RA/A1, and the BDR predictor must genuinely
depend on the barrier (the old compute_Rp dressed RA/A1 in BDR algebra that
cancelled to be barrier-independent).
"""
from __future__ import annotations

import copy

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.initialize import compute_Rp, resistance_area_bdr


def test_compute_Rp_is_RA_over_A1():
    cc = PhysicalConstantsConfig()
    assert abs(compute_Rp(cc) - cc.RA / cc.A1) < 1e-9
    assert 4.0e3 < compute_Rp(cc) < 6.0e3        # ~5 kOhm at the default RA


def test_compute_Rp_independent_of_barrier():
    """R_P = RA/A1 does not depend on phi_bar/tox (RA is the handle)."""
    cc = PhysicalConstantsConfig()
    c2 = copy.deepcopy(cc); c2.phi_bar = 0.7; c2.tox = 2.0e-9
    assert abs(compute_Rp(cc) - compute_Rp(c2)) < 1e-9


def test_bdr_predictor_depends_on_barrier():
    """resistance_area_bdr MUST change with phi_bar and tox (genuine BDR)."""
    cc = PhysicalConstantsConfig()
    ra0 = resistance_area_bdr(cc)
    c_phi = copy.deepcopy(cc); c_phi.phi_bar = 0.5
    c_tox = copy.deepcopy(cc); c_tox.tox = 1.8e-9
    assert resistance_area_bdr(c_phi) > ra0 * 1.5      # higher barrier -> higher RA
    assert resistance_area_bdr(c_tox) > ra0 * 1.5      # thicker barrier -> higher RA


def test_bdr_predictor_matches_measured_RA_at_mgo_effective_mass():
    """With an MgO-like effective mass (~0.3 m_e), the Simmons predictor lands
    within a factor ~2 of the measured RA at the physical phi_bar=0.4 eV."""
    cc = PhysicalConstantsConfig()
    ra = resistance_area_bdr(cc, m_eff_ratio=0.3)
    assert 0.5 * cc.RA < ra < 2.0 * cc.RA
