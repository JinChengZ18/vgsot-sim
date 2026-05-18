"""
Smoke tests for the opt-in toggles introduced in the 2026-05-16 round:

  * `R_series`            — additive parasitic resistance in `tmr()`
  * `rng=...`             — end-to-end Generator plumbing
  * `rng_mode`            — `ser_sot_no_vcma_thermal` legacy/generator switch
  * `integrator=...`      — euler_spherical / cayley
  * `SerResult.psw`       — Psw = 1 − SER alias

The tests are fast (< 10 s) so they can run in CI on every push.
"""
from __future__ import annotations

import numpy as np
import pytest

from vgsot_sim.configs import PhysicalConstantsConfig, SerSotNoVcmaThermalConfig
from vgsot_sim.tmr import tmr
from vgsot_sim.initialize import compute_Rp
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation
from vgsot_sim.ser_cases import ser_sot_no_vcma_thermal


# ─────────────────────────────────────────────────────────────────────────
# 1. R_series toggle
# ─────────────────────────────────────────────────────────────────────────
def test_r_series_default_zero_does_not_change_R():
    """R_series=0 must give identical R to a fresh config."""
    cc = PhysicalConstantsConfig()
    assert cc.R_series == 0.0
    assert abs(tmr(0.0, +1, cc) - compute_Rp(cc)) < 1e-9


def test_r_series_adds_constant_offset():
    """Setting R_series=R0 must add exactly R0 to every R(m_z, V) value."""
    cc = PhysicalConstantsConfig()
    R_P_bare  = tmr(0.0, +1, cc)
    R_AP_bare = tmr(0.0, -1, cc)
    cc.R_series = 350.0
    assert abs(tmr(0.0, +1, cc) - R_P_bare  - 350.0) < 1e-6
    assert abs(tmr(0.0, -1, cc) - R_AP_bare - 350.0) < 1e-6


def test_r_series_include_flag_bypasses_offset():
    """`include_series=False` recovers the bare tunnel resistance."""
    cc = PhysicalConstantsConfig()
    R_P_bare = tmr(0.0, +1, cc)
    cc.R_series = 500.0
    assert abs(tmr(0.0, +1, cc, include_series=False) - R_P_bare) < 1e-9


# ─────────────────────────────────────────────────────────────────────────
# 2. rng= end-to-end (byte-reproducibility)
# ─────────────────────────────────────────────────────────────────────────
def _single_short_run(rng=None):
    cc = PhysicalConstantsConfig()
    sim_end = int(round(0.5e-9 / cc.t_step))
    return run_piecewise_direct_excitation(
        sim_start_step=1, sim_mid1_step=sim_end,
        sim_mid2_step=sim_end, sim_end_step=sim_end,
        pap=1, v_mtj_stage1=0, v_mtj_stage2=0, v_mtj_stage3=0,
        i_sot_stage1=-1000e-6, i_sot_stage2=0, i_sot_stage3=0,
        estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1,
        estt_stage3=0, esot_stage3=1,
        vnv=0, non=1, r_sot_fl_dl=0.83, constants=cc, show_progress=False,
        rng=rng,
    )


def test_rng_generator_byte_reproducible():
    """Same Generator(seed) → identical mz trajectory."""
    rng_a = np.random.default_rng(seed=42)
    rng_b = np.random.default_rng(seed=42)
    r_a = _single_short_run(rng=rng_a)
    r_b = _single_short_run(rng=rng_b)
    assert np.max(np.abs(r_a.mz - r_b.mz)) == 0.0


def test_rng_legacy_seed_still_works():
    """`rng=None` + np.random.seed(...) reproduces the legacy stream byte-for-byte."""
    np.random.seed(42); r_a = _single_short_run()
    np.random.seed(42); r_b = _single_short_run()
    assert np.max(np.abs(r_a.mz - r_b.mz)) == 0.0


# ─────────────────────────────────────────────────────────────────────────
# 3. Integrator toggle (euler_spherical vs cayley)
# ─────────────────────────────────────────────────────────────────────────
def test_integrator_toggle_runs():
    """Both integrators run end-to-end without raising."""
    cc = PhysicalConstantsConfig()
    sim_end = int(round(0.5e-9 / cc.t_step))
    common = dict(
        sim_start_step=1, sim_mid1_step=sim_end, sim_mid2_step=sim_end,
        sim_end_step=sim_end, pap=1,
        v_mtj_stage1=0, v_mtj_stage2=0, v_mtj_stage3=0,
        i_sot_stage1=-1000e-6, i_sot_stage2=0, i_sot_stage3=0,
        estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1,
        estt_stage3=0, esot_stage3=1,
        vnv=0, non=0, r_sot_fl_dl=0.83, constants=cc, show_progress=False,
    )
    np.random.seed(7); r_eul = run_piecewise_direct_excitation(integrator="euler_spherical", **common)
    np.random.seed(7); r_cay = run_piecewise_direct_excitation(integrator="cayley",          **common)
    # Norm should be preserved for both
    m_eul = np.array([np.sin(r_eul.theta) * np.cos(r_eul.phi),
                       np.sin(r_eul.theta) * np.sin(r_eul.phi),
                       r_eul.mz])
    m_cay = np.array([np.sin(r_cay.theta) * np.cos(r_cay.phi),
                       np.sin(r_cay.theta) * np.sin(r_cay.phi),
                       r_cay.mz])
    norms_eul = np.linalg.norm(m_eul, axis=0)
    norms_cay = np.linalg.norm(m_cay, axis=0)
    assert np.max(np.abs(norms_eul - 1.0)) < 1e-6
    assert np.max(np.abs(norms_cay - 1.0)) < 1e-6
    # Final mz should agree to within ~2e-3 at this short pulse (no noise)
    # — they DO diverge at long times near threshold; this test just verifies
    # both ran sanely on the same input.
    assert np.isfinite(r_eul.mz[-1]) and np.isfinite(r_cay.mz[-1])


# ─────────────────────────────────────────────────────────────────────────
# 4. SerResult.psw alias
# ─────────────────────────────────────────────────────────────────────────
def test_serresult_psw_alias_is_one_minus_ser():
    cc = PhysicalConstantsConfig()
    cfg = SerSotNoVcmaThermalConfig(
        i_sot_list=(-1500e-6,), trials=4,
        sim_start_step=1, sim_mid1_step=750, sim_end_step=4000,
        pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
        target_mz=1.0, failure_tol=0.2, constants=cc,
    )
    res = ser_sot_no_vcma_thermal(cfg, show_progress=False, seed=2026)
    np.testing.assert_allclose(res.psw, 1.0 - res.ser)


# ─────────────────────────────────────────────────────────────────────────
# 5. ser_cases rng_mode toggle does not raise
# ─────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("mode", ["legacy", "generator"])
def test_ser_rng_mode_runs(mode):
    cc = PhysicalConstantsConfig()
    cfg = SerSotNoVcmaThermalConfig(
        i_sot_list=(-1500e-6,), trials=3,
        sim_start_step=1, sim_mid1_step=750, sim_end_step=4000,
        pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
        target_mz=1.0, failure_tol=0.2, constants=cc,
    )
    res = ser_sot_no_vcma_thermal(cfg, show_progress=False, seed=2026, rng_mode=mode)
    assert 0.0 <= res.ser[0] <= 1.0


def test_ser_rng_mode_bad_value_raises():
    cc = PhysicalConstantsConfig()
    cfg = SerSotNoVcmaThermalConfig(
        i_sot_list=(-1500e-6,), trials=1,
        sim_start_step=1, sim_mid1_step=750, sim_end_step=4000,
        pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
        target_mz=1.0, failure_tol=0.2, constants=cc,
    )
    with pytest.raises(ValueError):
        ser_sot_no_vcma_thermal(cfg, show_progress=False, seed=2026, rng_mode="nope")


# ─────────────────────────────────────────────────────────────────────────
# 6. run_piecewise_terminal_voltage: integrator toggle works end-to-end
# ─────────────────────────────────────────────────────────────────────────
def test_terminal_voltage_integrator_runs():
    """`run_piecewise_terminal_voltage(integrator="cayley")` runs end-to-end
    and preserves |m|=1; the spherical-Euler branch likewise runs and
    yields a finite final m_z."""
    from vgsot_sim.configs import TerminalVoltageControlConfig
    from vgsot_sim.time_series_cases import run_piecewise_terminal_voltage

    cc = PhysicalConstantsConfig()
    sim_end = int(round(0.5e-9 / cc.t_step))
    cfg = TerminalVoltageControlConfig(
        sim_start_step=1, sim_mid1_step=sim_end // 2,
        sim_mid2_step=sim_end, sim_end_step=sim_end,
        pap=1,
        v_stage1=(0.5, -0.5, 0.0),
        v_stage2=(0.0, 0.0, 0.0),
        v_stage3=(0.0, 0.0, 0.0),
        estt_stage1=0, esot_stage1=1,
        estt_stage2=0, esot_stage2=1,
        estt_stage3=0, esot_stage3=1,
        vnv=0, non=1, r_sot_fl_dl=0.83,
        constants=cc,
    )
    np.random.seed(11)
    r_eul = run_piecewise_terminal_voltage(cfg, show_progress=False,
                                            integrator="euler_spherical")
    np.random.seed(11)
    r_cay = run_piecewise_terminal_voltage(cfg, show_progress=False,
                                            integrator="cayley")
    # Both should produce finite trajectories with bounded m_z
    assert np.all(np.isfinite(r_eul.mz)) and np.all(np.isfinite(r_cay.mz))
    # Cayley norm preservation
    m_cay = np.array([np.sin(r_cay.theta) * np.cos(r_cay.phi),
                       np.sin(r_cay.theta) * np.sin(r_cay.phi),
                       r_cay.mz])
    norms = np.linalg.norm(m_cay, axis=0)
    assert np.max(np.abs(norms - 1.0)) < 1e-6


# ─────────────────────────────────────────────────────────────────────────
# 7. Cayley path is byte-reproducible when fed an explicit rng
# ─────────────────────────────────────────────────────────────────────────
def test_cayley_rng_byte_reproducible():
    """Same Generator(seed) through the Cayley path → identical m_z."""
    cc = PhysicalConstantsConfig()
    sim_end = int(round(0.5e-9 / cc.t_step))
    common = dict(
        sim_start_step=1, sim_mid1_step=sim_end, sim_mid2_step=sim_end,
        sim_end_step=sim_end, pap=1,
        v_mtj_stage1=0, v_mtj_stage2=0, v_mtj_stage3=0,
        i_sot_stage1=-1000e-6, i_sot_stage2=0, i_sot_stage3=0,
        estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1,
        estt_stage3=0, esot_stage3=1,
        vnv=0, non=1, r_sot_fl_dl=0.83, constants=cc, show_progress=False,
        integrator="cayley",
    )
    rng_a = np.random.default_rng(seed=99)
    rng_b = np.random.default_rng(seed=99)
    r_a = run_piecewise_direct_excitation(rng=rng_a, **common)
    r_b = run_piecewise_direct_excitation(rng=rng_b, **common)
    assert np.max(np.abs(r_a.mz - r_b.mz)) == 0.0
