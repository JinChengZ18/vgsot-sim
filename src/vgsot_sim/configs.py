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
    # Values follow CODATA 2018 recommendations (kept to 8 significant figures
    # where measurement allows; μ0 and kB are exact in the 2019 SI redefinition
    # but retain the 4π convention for μ0 for legacy compatibility).
    u0: float = 1.25663706e-6   # Vacuum permeability μ0 (H/m)
    e: float = 1.60217663e-19   # Elementary charge (C)
    h_bar: float = 1.05457182e-34  # Reduced Planck constant ħ (J·s)
    uB: float = 9.27401008e-24  # Bohr magneton μB (J/T)
    kb: float = 1.380649e-23    # Boltzmann constant kB (J/K), SI-exact
    m: float = 9.10938371e-31   # Electron mass (kg)

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
    Ki: float = 0.32e-3         # Interfacial anisotropy energy density at RT (J/m^2)
    Ms: float = 0.625e6         # Saturation magnetization at RT (A/m)
    beta: float = 60e-15        # VCMA coefficient (J/V·m)
    tf: float = 1.1e-9          # Free layer thickness (m)
    tox: float = 1.4e-9         # Tunnel barrier thickness (m)
    D: float = 80e-9            # MTJ physical diameter D_phys (m)
    D_elec: float = 65e-9       # MTJ electrical effective diameter (m) — used by BDR and ellipsoid demag
    t_mtj: float = 20e-9        # Total MTJ pillar thickness (m), for RC thermal time constant
    T: float = 300              # Operating temperature (K)
    T_RT: float = 300           # Reference (room) temperature (K) for M_s/K_i/eta scaling
    T_C: float = 1100           # Curie temperature (K) for Bloch / Callen–Callen scaling
    cc_exponent: float = 2.18   # modified Callen–Callen exponent for K_i(T)
    H_k_eff_RT: float = 5.0e4   # Effective anisotropy field at RT (A/m); H_k from FMR ≈ 500 Oe = 39.8 kA/m baseline; here scaled to give K_i(RT) ≈ 0.32 mJ/m² via the calibration formula
    alpha: float = 0.05         # Gilbert damping constant (unitless)
    t_step: float = 1e-12       # Simulation time step (s)
    # MgO thermal properties (used by thermal.self_heating, optional)
    Cv: float = 2.5e6           # Volumetric heat capacity (J/(m^3·K))
    lambda_MgO: float = 4.0     # MgO thermal conductivity (W/(m·K)), thin-film value
    t_MgO: float = 1.4e-9       # MgO thickness used as thermal path (m); defaults to tox

    # =========================
    # 4. Switching / transport parameters
    # =========================
    Vh: float = 0.5             # TMR half-bias voltage (V), reduces TMR by 50%
    P: float = 0.58             # Spin polarization (unitless, used in STT)
    # Effective spin Hall angle. The original value 0.25 came from chapter
    # §2.2.2.5 as a textbook β-W value, but the simplified spherical-Euler
    # LLG + lumped thermal model does not capture all the SOT-channel
    # dissipation paths (Néel-Edelstein, interfacial spin memory loss,
    # parasitic series resistance, etc.). When the deterministic threshold
    # V_th(0.75 ns) is matched to the experimental Device A P→AP detailed-
    # P_sw value of 894 mV (§2.3.3 Sigmoid fit), a calibrated effective
    # θ_SH ≈ 0.066 is required (Cayley integrator; 0.04 was the pre-FL-SOT-fix
    # value against the spherical-Euler stepper). We use this calibrated value as the default
    # so figures generated from vgsot-sim quantitatively reproduce the
    # detailed P_sw curve at the 0.75 ns operating point; the literature
    # β-W value 0.25 is recovered via `theta_SH=0.25` for material studies.
    theta_SH: float = 0.066     # Spin Hall angle (effective; recalibrated to V_th @ 0.75 ns
                                # with the FL-SOT-corrected Cayley integrator — was 0.04 against
                                # the pre-fix spherical-Euler stepper; see docs/maintenance/version_notes.md 2026-06)
    phi_bar: float = 0.4        # Barrier height for tunneling (eV)
    TMR: float = 1.0            # Tunnel magnetoresistance ratio (unitless, ≈100% matching the §2.3.3 hysteresis amplitude R_AP/R_P ≈ 2)
    RA: float = 16.6e-12        # Resistance-area product (Ω·m²) — calibrated to R_P ≈ 5 kΩ at D_elec = 65 nm
                                # (raw experimental hysteresis loops show low-R state at ≈4.9 kΩ, high-R state at ≈10.0 kΩ)

    # ─── Series-resistance toggle ────────────────────────────────────────
    # Additive parasitic resistance applied on top of the intrinsic MTJ
    # tunnel resistance returned by tmr(): R_obs(m_z, V) = R_MTJ(m_z, V) + R_series.
    # Default 0 keeps the legacy "pure-MTJ" behaviour byte-for-byte; set to
    # the contact + lead value to reproduce wafer-level R(V) measurements.
    # See docs/technical_details.md §2.6.
    R_series: float = 0.0       # Ω

    # =========================
    # 5. Bias / exchange field
    # =========================
    # Sign / direction convention (CRITICAL — see docs/technical_details.md
    # §1.4 and §2.7):
    #   The closed-form LLG in `dynamic_switching.py` is expanded with the
    #   implicit spin-Hall polarisation σ_SH = -x̂. For SOT switching to be
    #   deterministic at a given pulse polarity, the in-plane bias field
    #   `H_ex` MUST be PERPENDICULAR to σ_SH (i.e. along ±y or ±z); a field
    #   collinear with σ_SH (along ±x) leaves the equator-equilibrium m-state
    #   symmetric about the easy axis and the final m_z lands in either
    #   basin with ~50% probability, turning the SER curve from a clean
    #   sigmoid into a flat ~0.5 plateau (back-hopping-like random-bit
    #   regime). The chapter §2.3.3 phrase "200 Oe along the current
    #   direction" maps to the simulator's `h_ex_y` (with the appropriate
    #   sign that gives the desired switching chirality) because in the
    #   code's convention current is along ±y.
    #
    # The 50-Oe magnitude (vs the experimental 200 Oe) is a simulator-
    # tuned value: smaller |H_ex| keeps the in-plane equilibrium close to
    # σ_SH, which the simplified SOT/FL torque expansion handles cleanly.
    # The full experimental 200 Oe also gives correct chirality but
    # exhibits a higher residual SER tail (~0.4) due to mechanisms (Néel-
    # type DMI, interface roughness, additional FL components) not yet in
    # this model.
    h_ex_x: float = 0.0                        # A/m
    h_ex_y: float = -50.0 * 1000.0 / (4 * pi)  # A/m — −50 Oe perpendicular to σ_SH
    h_ex_z: float = 0.0                        # A/m

    # =========================
    # 6. TMR bias-dependence model
    # =========================
    # tmr_model: "lorentzian" or "pdk".
    #   - "lorentzian": TMR(V) = TMR0 / (1 + (V/Vh)^2), single-parameter form
    #     (Zhang et al. Verilog-A), kept for legacy compatibility.
    #   - "pdk": TMR(V) = TMR0/k_TMR * [1/(a V^2 + b|V| + c) - 1], three-
    #     parameter quadratic-rational form from Hikstor PDK extraction
    #     (Doc Section 2.2.2.3, parameters in Table 2.2.2.5).
    tmr_model: str = "pdk"
    k_tmr: float = 1.2346       # PDK normalisation (dimensionless)
    a_tmr: float = 0.1729       # PDK V_MTJ^2 coefficient (1/V^2)
    b_tmr: float = 0.1315       # PDK |V_MTJ| coefficient (1/V)
    c_tmr: float = 0.4475       # PDK constant offset (dimensionless)

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
        """Electrical / magnetic effective MTJ area, using D_elec.

        Resistance area product RA, BDR tunnel resistance and the magnetic
        active volume all refer to the electrical-effective area; the
        edge-etch damaged rim contributes neither carriers nor coherent
        magnetic switching. Use `A1_phys` for the geometric area when a
        thermal-pillar surface or area-based heat source is needed.
        """
        return pi * (self.D_elec ** 2) / 4

    @property
    def A1_phys(self) -> float:
        """Geometric MTJ area, using D_phys (for thermal-pillar heat source)."""
        return pi * (self.D ** 2) / 4

    @property
    def v(self) -> float:
        """Magnetic active volume = t_f * A1(D_elec)."""
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


# ─────────────────────────────────────────────────────────────────────────
# Default protocol for all SOT-only cases is the chapter §2.3.3 same-batch
# detailed-P_sw measurement on Device A:
#   - write pulse `t_w = 0.75 ns`  → sim_mid1_step = 750  (t_step = 1 ps)
#   - relaxation tail 3.25 ns      → sim_end_step  = 4000
#   - V_MTJ = 0  (no VCMA in this protocol)
#   - I_SOT swept across the calibrated threshold |I_th| ≈ 1.1 mA
#     (theta_SH = 0.066 calibration; see configs.PhysicalConstantsConfig
#     and docs/technical_details.md §2.3)
# Override fields explicitly when running a different protocol.
# ─────────────────────────────────────────────────────────────────────────


@dataclass
class SotOnlyConstantCurrentConfig:
    sim_start_step: int = 1
    sim_mid1_step: int = 750       # 0.75 ns write pulse (matches §2.3.3 detailed-P_sw)
    sim_end_step: int = 4000       # + 3.25 ns relaxation tail
    pap: int = 1                   # 1 ≡ AP start; success state is `target_mz = +1` (P)

    i_sot_stage1: float = -1500e-6 # super-threshold (~1.4× I_th) for a clean deterministic switch
    i_sot_stage2: float = 0.0
    v_mtj_stage1: float = 0.0      # §2.3.3 protocol drives SOT only; VCMA disabled
    v_mtj_stage2: float = 0.0

    vnv: int = 0                   # VCMA term disabled
    non: int = 1                   # thermal noise enabled (matches experiment)
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 5e-10
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)


@dataclass
class SotSwitchingNoVcmaConfig:
    # Bracket the calibrated I_th ≈ 1.1 mA so the sigmoid is sampled either side.
    i_sot_list: Sequence[float] = (
        -1500e-6, -1300e-6, -1200e-6, -1100e-6, -1000e-6, -900e-6, -700e-6,
    )
    sim_start_step: int = 1
    sim_mid1_step: int = 750
    sim_end_step: int = 4000
    pap: int = 1

    non: int = 1
    v_mtj: float = 0.0
    i_sot_relax: float = 0.0
    vnv: int = 0
    r_sot_fl_dl: float = 0.83
    tick_spacing_s: float = 5e-10
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)


@dataclass
class SerSotNoVcmaThermalConfig:
    # Bracket the calibrated I_th ≈ 1.1 mA so the P_sw sigmoid is sampled either side.
    i_sot_list: Sequence[float] = (
        -1500e-6, -1300e-6, -1200e-6, -1100e-6, -1000e-6, -900e-6, -700e-6,
    )
    trials: int = 100              # interactive-friendly; chapter Chapter02_local_10 uses 80
    sim_start_step: int = 1
    sim_mid1_step: int = 750
    sim_end_step: int = 4000
    pap: int = 1

    non: int = 1
    v_mtj: float = 0.0
    vnv: int = 0
    r_sot_fl_dl: float = 0.83      # unified with the other SOT-only cases (was 0.0)
    target_mz: float = 1.0         # pap=1 starts AP; success = end at P (m_z = +1)
    failure_tol: float = 0.2       # matches the chapter-figure tolerance
    constants: PhysicalConstantsConfig = field(default_factory=PhysicalConstantsConfig)
