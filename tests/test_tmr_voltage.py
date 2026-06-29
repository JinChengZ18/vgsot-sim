"""Pin the TMR(V) PDK / Lorentzian forms and the R(m_z) conductance endpoints."""
from __future__ import annotations

import copy

import pytest

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.initialize import compute_Rp
from vgsot_sim.tmr import tmr, tmr_eff


def test_pdk_tmr_eff_is_TMR_at_zero_bias_and_decays():
    cc = PhysicalConstantsConfig()          # default tmr_model = "pdk"
    assert cc.tmr_model == "pdk"
    assert abs(tmr_eff(0.0, cc) - cc.TMR) < 1e-3        # (1/c - 1)*TMR/k_tmr ~= TMR (PDK norm rounding)
    assert tmr_eff(0.3, cc) < tmr_eff(0.0, cc)          # decays with |V|
    assert abs(tmr_eff(0.3, cc) - tmr_eff(-0.3, cc)) < 1e-12   # even in V


def test_lorentzian_half_at_Vh():
    cc = copy.deepcopy(PhysicalConstantsConfig()); cc.tmr_model = "lorentzian"
    assert abs(tmr_eff(0.0, cc) - cc.TMR) < 1e-12
    assert abs(tmr_eff(cc.Vh, cc) - cc.TMR / 2.0) < 1e-12


def test_R_endpoints_parallel_and_antiparallel():
    cc = PhysicalConstantsConfig()
    Rp = compute_Rp(cc)
    T0 = tmr_eff(0.0, cc)
    assert abs(tmr(0.0, +1.0, cc) - Rp) < 1e-6                  # m_z=+1 -> R_P
    assert abs(tmr(0.0, -1.0, cc) - Rp * (1.0 + T0)) < 1e-6     # m_z=-1 -> R_AP


def test_unknown_tmr_model_raises():
    cc = copy.deepcopy(PhysicalConstantsConfig()); cc.tmr_model = "nope"
    with pytest.raises(ValueError):
        tmr_eff(0.1, cc)
