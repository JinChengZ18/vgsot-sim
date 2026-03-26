from math import cos, sin

from .anisotropy import field
from .configs import PhysicalConstantsConfig


def switching(V_MTJ, I_SOT, R_MTJ, theta, phi, ESTT, ESOT, VNV=1, NON=0, R_SOT_FL_DL=0.83, R_STT_FL_DL=0, constants: PhysicalConstantsConfig | None = None):
    """Dynamic-switching module and solving LLG analytically."""
    if constants is None:
        raise ValueError("constants must be provided.")

    n = 1
    ENE = 1
    H_EFF, _ = field(theta, phi, V_MTJ, n, NON, ENE, VNV, constants)

    I_MTJ = V_MTJ / R_MTJ
    J_STT = I_MTJ / constants.A1
    J_SOT = I_SOT / constants.A2

    gamma_red = constants.gamma / (1 + constants.alpha**2)
    H_DL_STT = ESTT * constants.h_bar * constants.P * J_STT / (2 * constants.e * constants.u0 * constants.Ms * constants.tf)
    H_FL_STT = R_STT_FL_DL * H_DL_STT
    H_DL_SOT = ESOT * constants.h_bar * constants.theta_SH * J_SOT / (2 * constants.e * constants.u0 * constants.Ms * constants.tf)
    H_FL_SOT = R_SOT_FL_DL * H_DL_SOT

    dtheta_dt = gamma_red*(
        H_EFF[0]*(constants.alpha*cos(theta)*cos(phi) - sin(phi))
        + H_EFF[1]*(constants.alpha*cos(theta)*sin(phi) + cos(phi))
        - H_EFF[2]*(constants.alpha*sin(theta))
        + sin(theta)*(constants.alpha*H_FL_STT - H_DL_STT)
        - H_DL_SOT*(constants.alpha*sin(phi) + cos(theta)*cos(phi))
        + H_FL_SOT*(constants.alpha*cos(theta)*cos(phi)-sin(phi))
    )

    dphi_dt = gamma_red*(
        1/sin(theta)*(
            H_EFF[0]*(-constants.alpha*sin(phi) - cos(theta)*cos(phi))
            + H_EFF[1]*(constants.alpha*cos(phi)-cos(theta)*sin(phi))
            + H_EFF[2]*sin(theta)
        )
        - (constants.alpha*H_DL_STT + H_FL_STT)
        - H_DL_SOT/sin(theta)*(constants.alpha*cos(theta)*cos(phi)-sin(phi))
        - H_FL_SOT/sin(theta)*(constants.alpha*sin(phi)-cos(theta)*cos(phi))
    )

    theta_1 = theta + dtheta_dt * constants.t_step
    phi_1 = phi + dphi_dt * constants.t_step
    mz = cos(theta_1)
    return (mz, phi_1, theta_1)
