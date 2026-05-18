"""Temperature-sweep and self-heating-transient cases.

These are the simulation primitives behind the chapter-2 figures:

- `sweep_material_temperature()`  → §2.2.2.2 fig 2.13 panels (a)–(c): T-dependent
  M_s, K_i, η on the same temperature axis.
- `tmr_voltage_sweep()`           → §2.2.2.3 fig 2.13(d): bias-dependent TMR
  for both Lorentzian and PDK forms.
- `thermal_transient()`           → §2.2.2.1 fig 2.12: lumped RC transient T(t)
  under a constant operating voltage. Multiple call sites generate the
  STT-only / SOT-only / STT+SOT panel comparison.

The chapter-04 helper scripts now consume these instead of re-deriving
the math inline.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .configs import PhysicalConstantsConfig
from . import material_temperature as mt
from . import tmr as tmr_mod
from . import thermal


@dataclass
class MaterialTemperatureSweep:
    """Per-temperature snapshot of magnetic / spin-transport parameters."""
    T_K: np.ndarray
    M_s: np.ndarray             # A/m
    K_i: np.ndarray             # J/m^2
    eta: np.ndarray             # unitless (spin polarisation)
    K_U_eff: np.ndarray         # J/m^3 (interface + demag-corrected)


def sweep_material_temperature(
    constants: PhysicalConstantsConfig | None = None,
    T_min: float = 250.0,
    T_max: float = 500.0,
    n_points: int = 251,
    eta_RT: float | None = None,
) -> MaterialTemperatureSweep:
    """Evaluate M_s(T), K_i(T), η(T), K_U^eff(T) on a temperature grid.

    Defaults to the chapter-2 temperature window 250–500 K.
    """
    if constants is None:
        constants = PhysicalConstantsConfig()
    T_K = np.linspace(T_min, T_max, n_points)
    M_s = np.array([mt.ms_of_T(float(T), constants) for T in T_K])
    K_i = np.array([mt.ki_of_T(float(T), constants) for T in T_K])
    eta = np.array([mt.eta_of_T(float(T), constants, eta_RT=eta_RT) for T in T_K])
    K_U = np.array([mt.k_u_eff_of_T(float(T), constants) for T in T_K])
    return MaterialTemperatureSweep(T_K=T_K, M_s=M_s, K_i=K_i, eta=eta, K_U_eff=K_U)


@dataclass
class TMRVoltageSweep:
    V: np.ndarray
    TMR_pdk: np.ndarray
    TMR_lor: np.ndarray
    R_at_AP_pdk: np.ndarray
    R_at_AP_lor: np.ndarray
    R_P: float


def tmr_voltage_sweep(
    constants: PhysicalConstantsConfig | None = None,
    V_min: float = 0.0,
    V_max: float = 1.5,
    n_points: int = 301,
) -> TMRVoltageSweep:
    """Evaluate TMR(V) under both Lorentzian and PDK models."""
    if constants is None:
        constants = PhysicalConstantsConfig()
    V = np.linspace(V_min, V_max, n_points)
    saved = constants.tmr_model
    try:
        constants.tmr_model = "pdk"
        tmr_pdk = np.array([tmr_mod.tmr_eff(float(v), constants) for v in V])
        R_AP_pdk = np.array([tmr_mod.tmr(float(v), -1.0, constants) for v in V])
        constants.tmr_model = "lorentzian"
        tmr_lor = np.array([tmr_mod.tmr_eff(float(v), constants) for v in V])
        R_AP_lor = np.array([tmr_mod.tmr(float(v), -1.0, constants) for v in V])
    finally:
        constants.tmr_model = saved
    from .initialize import compute_Rp
    return TMRVoltageSweep(
        V=V, TMR_pdk=tmr_pdk, TMR_lor=tmr_lor,
        R_at_AP_pdk=R_AP_pdk, R_at_AP_lor=R_AP_lor,
        R_P=compute_Rp(constants),
    )


@dataclass
class ThermalTransientResult:
    """T(t) under one operating-bias configuration."""
    t_seconds: np.ndarray
    T_K: np.ndarray
    T_eq_K: float
    tau_th_s: float
    Q_total_per_area: float
    label: str


def thermal_transient(
    constants: PhysicalConstantsConfig | None = None,
    V_MTJ: float = 0.0,
    V_SOT: float = 0.0,
    R_MTJ: float | None = None,
    R_SOT: float | None = None,
    T_0: float = 300.0,
    t_seconds=None,
    label: str = "",
) -> ThermalTransientResult:
    """Closed-form T(t) under constant V_MTJ / V_SOT.

    `R_MTJ` defaults to the BDR-derived parallel-state resistance; pass
    a measured value to match a specific operating point. `R_SOT`
    defaults to the geometric ρ L / (w d).
    """
    if constants is None:
        constants = PhysicalConstantsConfig()
    if R_MTJ is None:
        from .initialize import compute_Rp
        R_MTJ = compute_Rp(constants)
    Q = thermal.q_mtj(V_MTJ, R_MTJ, constants) + thermal.q_sot(V_SOT, constants, R_SOT=R_SOT)
    T_eq = thermal.steady_state_temperature(Q, T_0, constants)
    tau = thermal.thermal_time_constant(constants)
    if t_seconds is None:
        t_seconds = np.linspace(0.0, 100e-12, 600)
    else:
        t_seconds = np.asarray(t_seconds, dtype=float)
    T_curve = thermal.transient_temperature(t_seconds, Q, T_0, constants)
    return ThermalTransientResult(
        t_seconds=t_seconds, T_K=T_curve,
        T_eq_K=T_eq, tau_th_s=tau,
        Q_total_per_area=Q, label=label,
    )
