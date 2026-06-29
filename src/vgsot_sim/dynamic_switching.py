from math import cos, sin

from .anisotropy import field
from .configs import PhysicalConstantsConfig


def switching(V_MTJ, I_SOT, R_MTJ, theta, phi, ESTT, ESOT, VNV=1, NON=0,
              R_SOT_FL_DL=0.83, R_STT_FL_DL=0,
              constants: PhysicalConstantsConfig | None = None,
              Ki_T: float | None = None, Ms_T: float | None = None,
              T: float | None = None, h_th_ext=None,
              demag_mode: str = "ellipsoid",
              rng=None):
    """One-step explicit-Euler LLG update in spherical coordinates.

    Computes the new (m_z, phi, theta) after a single time step `constants.t_step`,
    given the current angles and the electrical excitation. The torque sources
    are (in order of accumulation): the effective field built by `anisotropy.field()`
    (PMA + VCMA + demag + exchange-bias + optional thermal noise), the
    spin-transfer torque (STT) driven by V_MTJ/R_MTJ through the MTJ, and the
    spin-orbit torque (SOT) driven by I_SOT through the heavy-metal channel.

    Sign / direction conventions inherited from the upstream Verilog-A model
    (Zhang et al., IEEE Access 2020):

    - STT spin polarisation `m_p = +z_hat` (perpendicular reference layer).
    - SOT spin polarisation `sigma_SH = -x_hat` (in-plane, anti-parallel to
      the assumed sign of the spin-Hall accumulation; verified against the
      DL part of dtheta/dt by hand derivation).
    - DL torque uses the "physicist" sign:
          tau_DL = -gamma * H_DL * m x (m x sigma)
      i.e. positive `H_DL_SOT` (with `theta_SH > 0`, `J_SOT > 0`) drives the
      magnetisation toward sigma.
    - FL torque uses the "positive Gilbert" sign:
          tau_FL = +gamma * H_FL * m x sigma

    FL-SOT dphi/dt sign — FIXED (2026-06):
    The closed-form expansion below matches the standard Landau-Lifshitz
    derivation for the H_eff, STT, and DL-SOT contributions to both
    dtheta/dt and dphi/dt, and for the FL-SOT contribution to dtheta/dt.
    The FL-SOT contribution to dphi/dt previously carried the wrong sign on
    the `cos(theta)*cos(phi)` term (an O(1) coupling weighted by
    `R_SOT_FL_DL ~ 0.83`); it has been corrected to `+cos(theta)*cos(phi)`,
    so this spherical-Euler stepper now shares the SAME right-hand side as
    the Cartesian/Cayley vector stepper (`dynamic_switching_vector`). The
    agreement is regression-tested in `tests/test_integrator_consistency.py`
    (dm/dt match to <1e-4 as dt -> 0). The default integrator is now Cayley.

    Parameters
    ----------
    V_MTJ : float
        Voltage across the MTJ (V). Drives STT and VCMA.
    I_SOT : float
        Current through the SOT channel (A). Drives SOT.
    R_MTJ : float
        Current MTJ resistance (Ω). Combined with V_MTJ to produce I_MTJ.
    theta, phi : float
        Current spherical angles (rad).
    ESTT, ESOT : int
        0/1 flags to enable or disable STT and SOT respectively.
    VNV : int, default 1
        0/1 flag to enable VCMA effect inside `anisotropy.field()`.
    NON : int, default 0
        0/1 flag to enable thermal noise inside `anisotropy.field()`.
    R_SOT_FL_DL, R_STT_FL_DL : float
        Field-like / damping-like ratio for SOT and STT respectively.
    constants : PhysicalConstantsConfig
        Device and physical-constant container.

    Returns
    -------
    (m_z, phi_new, theta_new)
    """
    if constants is None:
        raise ValueError("constants must be provided.")

    n = 1
    ENE = 1
    H_EFF, _ = field(theta, phi, V_MTJ, n, NON, ENE, VNV, constants,
                     demag_mode=demag_mode, Ki_T=Ki_T, Ms_T=Ms_T,
                     T=T, h_th_ext=h_th_ext, rng=rng)

    # SOT/STT effective fields scale as 1/M_s; use the T-corrected value
    # when provided so that self-heating reduces the magnetisation that the
    # spin current must rotate.
    Ms_use = constants.Ms if Ms_T is None else Ms_T

    I_MTJ = V_MTJ / R_MTJ
    J_STT = I_MTJ / constants.A1
    J_SOT = I_SOT / constants.A2

    gamma_red = constants.gamma / (1 + constants.alpha**2)
    H_DL_STT = ESTT * constants.h_bar * constants.P * J_STT / (2 * constants.e * constants.u0 * Ms_use * constants.tf)
    H_FL_STT = R_STT_FL_DL * H_DL_STT
    H_DL_SOT = ESOT * constants.h_bar * constants.theta_SH * J_SOT / (2 * constants.e * constants.u0 * Ms_use * constants.tf)
    H_FL_SOT = R_SOT_FL_DL * H_DL_SOT

    dtheta_dt = gamma_red*(
        H_EFF[0]*(constants.alpha*cos(theta)*cos(phi) - sin(phi))
        + H_EFF[1]*(constants.alpha*cos(theta)*sin(phi) + cos(phi))
        - H_EFF[2]*(constants.alpha*sin(theta))
        + sin(theta)*(constants.alpha*H_FL_STT - H_DL_STT)
        - H_DL_SOT*(constants.alpha*sin(phi) + cos(theta)*cos(phi))
        + H_FL_SOT*(constants.alpha*cos(theta)*cos(phi)-sin(phi))
    )

    # Guard against the sin(theta) -> 0 singularity in the spherical
    # representation. Near the poles, phi is geometrically undefined and the
    # 1/sin(theta) terms amplify floating-point roundoff into spurious
    # precession. Freezing phi for one step is the standard remedy; the
    # magnetization continues to evolve through theta only, and phi recovers
    # naturally as soon as theta moves away from the pole.
    sin_theta = sin(theta)
    if abs(sin_theta) < 1e-8:
        dphi_dt = 0.0
    else:
        inv_sin = 1.0 / sin_theta
        dphi_dt = gamma_red*(
            inv_sin*(
                H_EFF[0]*(-constants.alpha*sin(phi) - cos(theta)*cos(phi))
                + H_EFF[1]*(constants.alpha*cos(phi)-cos(theta)*sin(phi))
                + H_EFF[2]*sin(theta)
            )
            - (constants.alpha*H_DL_STT + H_FL_STT)
            - H_DL_SOT*inv_sin*(constants.alpha*cos(theta)*cos(phi)-sin(phi))
            - H_FL_SOT*inv_sin*(constants.alpha*sin(phi)+cos(theta)*cos(phi))
        )

    theta_1 = theta + dtheta_dt * constants.t_step
    phi_1 = phi + dphi_dt * constants.t_step
    mz = cos(theta_1)
    return (mz, phi_1, theta_1)
