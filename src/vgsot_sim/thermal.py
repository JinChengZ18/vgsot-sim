"""Self-heating dynamics for the SOT/STT/VGSOT device.

Implements the lumped one-dimensional thermal network derived in
§2.2.2.1 of the thesis:

    C_v · t_MTJ · dT/dt = Q_MTJ + Q_SOT - (λ_MgO / t_MgO)(T - T_0)

Q_MTJ is the Joule power per unit MTJ footprint area from the tunnel
current:
    Q_MTJ = V_MTJ² / R_MTJ · (4 / (π D_phys²))

Q_SOT is the heavy-metal channel dissipation rescaled to the MTJ
footprint via the height ratio t_MTJ / L_SOT:
    Q_SOT = V_SOT² / (R_SOT · A_SOT) · (t_MTJ / L_SOT)

Closed-form transient for constant Q_total:
    T(t) = T_0 + ΔT_eq · (1 − exp(−t/τ_th))
    ΔT_eq = Q_total · t_MgO / λ_MgO       (steady-state offset)
    τ_th  = C_v · t_MTJ · t_MgO / λ_MgO   (thermal time constant)
"""
from __future__ import annotations

from math import exp, pi
import numpy as np

from .configs import PhysicalConstantsConfig


def thermal_time_constant(constants: PhysicalConstantsConfig) -> float:
    """τ_th = C_v · t_MTJ · t_MgO / λ_MgO (seconds)."""
    return constants.Cv * constants.t_mtj * constants.t_MgO / constants.lambda_MgO


def q_mtj(V_MTJ: float, R_MTJ: float, constants: PhysicalConstantsConfig) -> float:
    """Joule heat per MTJ-pillar footprint from the tunnel current (W/m²).

    Q_MTJ = V_MTJ² / R_MTJ × (4 / (π D_elec²))

    Uses D_elec for the dissipation area, consistent with the assumption
    that the edge-damaged rim does not carry tunnel current (and hence
    does not dissipate). The lumped LHS still uses t_MTJ for the heat
    capacity, so this returns power per footprint area directly.
    """
    if R_MTJ <= 0:
        return 0.0
    return V_MTJ * V_MTJ / R_MTJ * (4.0 / (pi * constants.D_elec ** 2))


def q_sot(V_SOT: float, constants: PhysicalConstantsConfig,
          R_SOT: float | None = None) -> float:
    """Joule heat per MTJ-pillar footprint from the SOT channel current (W/m²).

    From thesis §2.2.2.1, the volumetric channel dissipation is
        q_vol = V_SOT² / (R_SOT · A_cross · L_SOT)     [W/m³]
    where A_cross = w · d is the channel **cross-section** (not the
    top-view footprint). Integrating over the MTJ-pillar thickness gives
    the per-footprint contribution
        Q_SOT = q_vol · t_MTJ
              = V_SOT² · t_MTJ / (R_SOT · w · d · L_SOT)
              = V_SOT² · t_MTJ / (R_SOT · A2 · L_SOT)
    where `A2 = w·d` is already available from PhysicalConstantsConfig.

    Pass `R_SOT` to override the geometric ρ L/(wd) value when matching a
    measured channel resistance.
    """
    if R_SOT is None:
        R_SOT = constants.R_W
    if R_SOT <= 0:
        return 0.0
    return V_SOT * V_SOT * constants.t_mtj / (R_SOT * constants.A2 * constants.l)


def steady_state_temperature(Q_total: float, T_0: float,
                             constants: PhysicalConstantsConfig) -> float:
    """T_eq = T_0 + Q_total · t_MgO / λ_MgO."""
    dT_eq = Q_total * constants.t_MgO / constants.lambda_MgO
    return T_0 + dT_eq


def transient_temperature(t_seconds, Q_total: float, T_0: float,
                          constants: PhysicalConstantsConfig):
    """Closed-form T(t) under constant Q_total. Accepts scalar or array t."""
    tau = thermal_time_constant(constants)
    dT_eq = Q_total * constants.t_MgO / constants.lambda_MgO
    t_seconds = np.asarray(t_seconds, dtype=float)
    return T_0 + dT_eq * (1.0 - np.exp(-t_seconds / tau))


def self_heating_step(T: float, V_MTJ: float, V_SOT: float, R_MTJ: float,
                      constants: PhysicalConstantsConfig,
                      dt: float | None = None, T_0: float = 300.0,
                      R_SOT: float | None = None) -> float:
    """One forward-Euler step of the lumped thermal equation.

        dT/dt = (Q_total - (λ_MgO/t_MgO)(T - T_0)) / (C_v · t_MTJ)

    Returns the updated temperature. Use this inside the magnetisation
    time-stepping loop when feeding T(t) into `material_temperature.*`.
    """
    if dt is None:
        dt = constants.t_step
    Q = q_mtj(V_MTJ, R_MTJ, constants) + q_sot(V_SOT, constants, R_SOT=R_SOT)
    cooling = (constants.lambda_MgO / constants.t_MgO) * (T - T_0)
    dT_dt = (Q - cooling) / (constants.Cv * constants.t_mtj)
    return T + dT_dt * dt
