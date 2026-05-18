"""Temperature dependence of CoFeB free-layer magnetic parameters.

Implements the closed-form scaling laws derived in §2.2.2.2 of the thesis:

- M_s(T): Bloch T^{3/2} law, normalised to room-temperature M_s(RT).
- η(T):  same Bloch scaling (Julliere mean-field), used to dress STT/SOT
         polarisation efficiencies.
- K_i(T): modified Callen–Callen scaling, with exponent fit to FMR data
         (default 2.18 for CoFeB/MgO interfaces).
- K_U_eff(T): combined interface anisotropy minus thin-film demag
              correction, mediated by ellipsoid demag factors (see
              `anisotropy.ellipsoid_demag`).

All scaling factors return a pure multiplicative correction. Callers
multiply their RT base value:

    Ms_T  = constants.Ms * ms_factor(T, constants)
    Ki_T  = constants.Ki * ki_factor(T, constants)
"""
from __future__ import annotations

from .configs import PhysicalConstantsConfig


def _bloch_factor(T: float, constants: PhysicalConstantsConfig) -> float:
    r"""Bloch T^{3/2} normalised scaling.

    Returns (1 - (T/T_C)^{3/2}) / (1 - (T_RT/T_C)^{3/2}) so that the factor
    equals exactly 1 at T = T_RT and drops to 0 at T = T_C.
    """
    if constants.T_C <= 0:
        raise ValueError(f"T_C must be positive, got {constants.T_C}")
    if T < 0:
        raise ValueError(f"T must be non-negative, got {T}")
    numer = 1.0 - (T / constants.T_C) ** 1.5
    denom = 1.0 - (constants.T_RT / constants.T_C) ** 1.5
    return numer / denom


def ms_factor(T: float, constants: PhysicalConstantsConfig) -> float:
    """Bloch-law scaling of saturation magnetisation."""
    return _bloch_factor(T, constants)


def eta_factor(T: float, constants: PhysicalConstantsConfig) -> float:
    """Bloch-law scaling of spin polarisation (Julliere mean field)."""
    return _bloch_factor(T, constants)


def ki_factor(T: float, constants: PhysicalConstantsConfig) -> float:
    """Modified Callen–Callen scaling of interface anisotropy K_i."""
    return _bloch_factor(T, constants) ** constants.cc_exponent


def ms_of_T(T: float, constants: PhysicalConstantsConfig) -> float:
    """M_s(T) in A/m using the RT calibration constants.Ms."""
    return constants.Ms * ms_factor(T, constants)


def ki_of_T(T: float, constants: PhysicalConstantsConfig) -> float:
    """K_i(T) in J/m^2 using the RT calibration constants.Ki."""
    return constants.Ki * ki_factor(T, constants)


def eta_of_T(T: float, constants: PhysicalConstantsConfig, eta_RT: float | None = None) -> float:
    """η(T), with optional explicit RT calibration. Defaults to constants.P."""
    if eta_RT is None:
        eta_RT = constants.P
    return eta_RT * eta_factor(T, constants)


def k_u_eff_of_T(T: float, constants: PhysicalConstantsConfig,
                 Nz_minus_Nx: float | None = None) -> float:
    """Effective uniaxial anisotropy energy density K_U^eff(T).

        K_U^eff(T) = K_i(T)/t_FL - 0.5 * mu_0 * M_s(T)^2 * (N_z - N_x)

    The first term is the interface contribution per unit free-layer
    thickness; the second is the shape-anisotropy penalty. Returns J/m^3.

    If `Nz_minus_Nx` is None, it is computed via the thin-disk approximation
    N_x = π t_f / (4 D_elec); call sites that have already evaluated the
    exact oblate-ellipsoid factors should pass them in to avoid double work.
    """
    from math import pi
    if Nz_minus_Nx is None:
        Nx = pi * constants.tf / (4 * constants.D_elec)
        Nz_minus_Nx = 1 - 3 * Nx
    Ki_T = ki_of_T(T, constants)
    Ms_T = ms_of_T(T, constants)
    return Ki_T / constants.tf - 0.5 * constants.u0 * Ms_T ** 2 * Nz_minus_Nx
