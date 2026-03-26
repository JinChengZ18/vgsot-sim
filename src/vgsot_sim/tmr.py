from .configs import PhysicalConstantsConfig
from .initialize import init


def tmr(V_MTJ, mz, constants: PhysicalConstantsConfig):
    """Tunnel Magnetoresistance module."""
    PAP = 0
    R_p, _, _, _ = init(PAP, constants)
    return R_p * (1 + (V_MTJ / constants.Vh) ** 2 + constants.TMR) / (1 + (V_MTJ / constants.Vh) ** 2 + constants.TMR * (0.5 * (1 + mz)))
