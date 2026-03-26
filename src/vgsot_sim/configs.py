from __future__ import annotations
from dataclasses import dataclass, field
from math import pi
from typing import Optional, Sequence, Tuple




@dataclass
class PhysicalConstantsConfig:
    """
    Physical constants and device parameters used in the simulation.

    Design principles:
    - All simulation parameters are centralized here for reproducibility
    - Parameters are grouped by physical meaning
    - Base parameters are user-configurable; derived quantities should be computed elsewhere

    Units are explicitly specified for each parameter.
    """

    # =========================
    # 1. Elementary constants
    # =========================
    u0: float = 12.56637e-7     # Vacuum permeability μ0 (H/m)
    e: float = 1.6e-19          # Elementary charge (C)
    h_bar: float = 1.054e-34    # Reduced Planck constant ħ (J·s)
    uB: float = 9.274e-24       # Bohr magneton μB (J/T)
    kb: float = 1.38e-23        # Boltzmann constant kB (J/K)
    m: float = 9.11e-31         # Electron mass (kg)

    # =========================
    # 2. Geometry / electronic parameters
    # =========================
    # R_SOT = 776Ω
    l: float = 240e-9            # Length of SOT channel (m)
    w: float = 200e-9            # Width of device (m)
    d: float = 4.3e-9            # Thickness of SOT channel (m)
    rho: float = 278e-8          # Resistivity of SOT material (Ω·m), e.g. β-IrMn

    # =========================
    # 3. Magnetic / material parameters
    # =========================
    Ki: float = 0.32e-3         # Interfacial anisotropy energy density (J/m^2)
    Ms: float = 0.625e6         # Saturation magnetization (A/m)
    beta: float = 60e-15        # VCMA coefficient (J/V·m)
    tf: float = 1.1e-9          # Free layer thickness (m)
    tox: float = 1.4e-9         # Tunnel barrier thickness (m)
    D: float = 80e-9            # MTJ diameter (m)
    T: float = 300              # Temperature (K)
    alpha: float = 0.05         # Gilbert damping constant (unitless)
    t_step: float = 1e-12       # Simulation time step (s)

    # =========================
    # 4. Switching / transport parameters
    # =========================
    Vh: float = 0.5             # TMR half-bias voltage (V), reduces TMR by 50%
    P: float = 0.58             # Spin polarization (unitless, used in STT)
    theta_SH: float = 0.25      # Spin Hall angle (unitless, SOT efficiency)
    phi_bar: float = 0.4        # Barrier height for tunneling (eV)
    TMR: float = 1.19           # Tunnel magnetoresistance ratio (unitless)
    RA: float = 36e-12          # Resistance-area product (Ω·m^2)

    # =========================
    # Notes
    # =========================
    # - Derived quantities such as:
    #     A1 (MTJ area), A2 (SOT cross-section),
    #     v (volume), R_W (wire resistance), Heff (effective field)
    #   should NOT be stored here as fixed values.
    #
    # - They must be computed dynamically from the base parameters above.
    #
    # - This ensures:
    #     ✅ consistency when parameters change
    #     ✅ no hidden coupling
    #     ✅ easier debugging and physical interpretation
    
    @property
    def gamma(self) -> float:
        return 2 * self.u0 * self.uB / self.h_bar

    @property
    def A2(self) -> float:
        return self.d * self.w

    @property
    def R_W(self) -> float:
        return self.rho * self.l / (self.w * self.d)

    @property
    def A1(self) -> float:
        return pi * (self.D ** 2) / 4

    @property
    def v(self) -> float:
        return self.tf * self.A1

    @property
    def Heff(self) -> float:
        return (2 * self.Ki) / (self.tf * self.Ms * self.u0)


@dataclass
class TerminalVoltageControlConfig:
    sim_start_step: int = 1
    sim_mid1_step: int = 2000
    sim_mid2_step: Optional[int] = None
    sim_end_step: int = 5000
    pap: int = 1

    v_stage1: Tuple[float, float, float] = (1.0, 0.0, 0.1)
    v_stage2: Tuple[float, float, float] = (-1.0, 0.0, 0.0)
    v_stage3: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    estt_stage1: int = 0
    esot_stage1: int = 1
    estt_stage2: int = 0
    esot_stage2: int = 1
    estt_stage3: int = 1
    esot_stage3: int = 1

    vnv: int = 1
    non: int = 1
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 5e-10
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)


@dataclass
class SotOnlyConstantCurrentConfig:
    sim_start_step: int = 1
    sim_mid1_step: int = 2000
    sim_end_step: int = 5000
    pap: int = 1

    i_sot_stage1: float = -400e-6
    i_sot_stage2: float = 0.0
    v_mtj_stage1: float = 0.0
    v_mtj_stage2: float = 0.0

    vnv: int = 1
    non: int = 1
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 5e-10
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)


@dataclass
class SotSwitchingNoVcmaConfig:
    i_sot_list: Sequence[float] = (-800e-6, -700e-6, -750e-6, -600e-6, -500e-6, -400e-6, -300e-6, -200e-6)
    sim_start_step: int = 1
    sim_mid1_step: int = 5000
    sim_end_step: int = 10000
    pap: int = 1

    non: int = 1
    v_mtj: float = 0.0
    i_sot_relax: float = 0.0
    vnv: int = 0
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 1e-9
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)


@dataclass
class SerSotNoVcmaThermalConfig:
    i_sot_list: Sequence[float] = (-800e-6, -700e-6, -750e-6, -600e-6, -500e-6, -400e-6, -300e-6, -200e-6)
    trials: int = 200
    sim_start_step: int = 1
    sim_mid1_step: int = 5000
    sim_end_step: int = 10000
    pap: int = 1

    non: int = 1
    v_mtj: float = 0.0
    vnv: int = 0
    r_sot_fl_dl: float = 0.0
    target_mz: float = 1.0
    failure_tol: float = 1e-1
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)


@dataclass
class VcmaAssistedSwitchingIsotSweepConfig:
    v_mtj: float = 1.2
    i_sot_list: Sequence[float] = (-60e-6, -40e-6, -20e-6, -10e-6)
    sim_start_step: int = 1
    sim_end_step: int = 10000
    pap: int = 1

    non: int = 1
    vnv: int = 1
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 1e-9
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)


@dataclass
class VcmaAssistedSwitchingVmtjSweepConfig:
    i_sot: Optional[float] = -20e-6
    v_mtj_list: Sequence[float] = (1.3, 1.2, 1.1, 1.0, 0.9)
    sim_start_step: int = 1
    sim_end_step: int = 10000
    pap: int = 1

    non: int = 1
    vnv: int = 1
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 1e-9
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)


@dataclass
class OptimizedVgsotSwitchingConfig:
    v_mtj_1: float = 1.5
    v_mtj_2: float = -1.0
    i_sot: Optional[float] = -20e-6
    t_pairs_s: Sequence[Tuple[float, float]] = (
        (25e-9, 0.0),
        (1.4e-9, 1.6e-9),
        (1.8e-9, 1.2e-9),
        (2.2e-9, 0.8e-9),
    )
    sim_total_time_s: float = 10e-9
    pap: int = 1

    non: int = 1
    vnv: int = 1
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 1e-9
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)


@dataclass
class SerOptimizedVgsotConfig:
    iterations_num: int = 300
    v_mtj_1: float = 1.5
    v_mtj_2: float = -1.0
    i_sot: Optional[float] = -20e-6
    t1_list_s: Sequence[float] = (1.3e-9, 1.4e-9, 1.5e-9, 1.6e-9, 1.7e-9, 1.8e-9, 1.9e-9)
    total_pulse_s: float = 3e-9
    sim_total_time_s: float = 10e-9
    pap: int = 1

    non: int = 1
    vnv: int = 1
    r_sot_fl_dl: float = 0.83
    target_final_mz: float = 1.0
    failure_tol: float = 1e-1
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)
