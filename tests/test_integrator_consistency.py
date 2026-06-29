"""Spherical-Euler vs Cayley integrator RHS-consistency.

Guards the fix for the FL-SOT dphi/dt sign anomaly in the spherical-Euler
stepper (the cos(theta)*cos(phi) term in dynamic_switching.py). Both
integrators must reduce to the SAME dm/dt as dt -> 0; before the fix the
FL-SOT contribution to dphi/dt disagreed by up to ~25% in the azimuthal
direction. We probe dm/dt = (m1 - m)/dt at a tiny step and require the two
steppers to agree to the O(dt) Euler-vs-Cayley truncation floor.
"""
from __future__ import annotations

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.dynamic_switching import switching
from vgsot_sim.dynamic_switching_vector import switching_vector


def _dmdt_both(theta, phi, *, V_MTJ, I_SOT, R_MTJ, ESTT, ESOT, R_SOT_FL_DL, cc):
    m = np.array([np.sin(theta) * np.cos(phi),
                  np.sin(theta) * np.sin(phi),
                  np.cos(theta)])
    mz1, phi1, theta1 = switching(V_MTJ, I_SOT, R_MTJ, theta, phi,
                                  ESTT=ESTT, ESOT=ESOT, VNV=0, NON=0,
                                  R_SOT_FL_DL=R_SOT_FL_DL, constants=cc)
    m1_eul = np.array([np.sin(theta1) * np.cos(phi1),
                       np.sin(theta1) * np.sin(phi1),
                       np.cos(theta1)])
    m1_cay = switching_vector(m, V_MTJ, I_SOT, R_MTJ, ESTT=ESTT, ESOT=ESOT,
                              VNV=0, NON=0, R_SOT_FL_DL=R_SOT_FL_DL, constants=cc)
    return (m1_eul - m) / cc.t_step, (m1_cay - m) / cc.t_step


def test_spherical_eq_cayley_sot_with_fieldlike():
    """FL-SOT active (R_SOT_FL_DL=0.83): the two steppers' dm/dt must agree."""
    cc = PhysicalConstantsConfig()
    cc.t_step = 1e-15
    for theta, phi in [(0.6, 0.8), (1.2, -2.0), (0.9, 2.5), (2.3, 1.1)]:
        dm_eul, dm_cay = _dmdt_both(theta, phi, V_MTJ=0.0, I_SOT=-1.0e-3,
                                    R_MTJ=5000.0, ESTT=0, ESOT=1,
                                    R_SOT_FL_DL=0.83, cc=cc)
        rel = np.linalg.norm(dm_eul - dm_cay) / np.linalg.norm(dm_cay)
        assert rel < 1e-4, f"(theta={theta},phi={phi}) rel mismatch {rel:.2e}"


def test_spherical_eq_cayley_no_fieldlike():
    """DL-SOT only (R_SOT_FL_DL=0): agreement holds (control case)."""
    cc = PhysicalConstantsConfig()
    cc.t_step = 1e-15
    for theta, phi in [(0.7, 0.5), (1.5, -1.0)]:
        dm_eul, dm_cay = _dmdt_both(theta, phi, V_MTJ=0.0, I_SOT=-1.0e-3,
                                    R_MTJ=5000.0, ESTT=0, ESOT=1,
                                    R_SOT_FL_DL=0.0, cc=cc)
        rel = np.linalg.norm(dm_eul - dm_cay) / np.linalg.norm(dm_cay)
        assert rel < 1e-4, f"(theta={theta},phi={phi}) rel mismatch {rel:.2e}"
