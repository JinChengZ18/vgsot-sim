from math import sqrt, exp, cos, pi

from .configs import PhysicalConstantsConfig


def init(PAP, constants: PhysicalConstantsConfig):
    """
    To initialize the VGSOT MTJ

    input:
    PAP = 1 if Parallel, 0 if Antiparallel

    output:
    R_MTJ, theta, mz, phi
    """

    if PAP == 1:
        theta = pi - sqrt((constants.kb * constants.T) / (constants.u0 * constants.Ms * constants.Heff * constants.v))
    else:
        theta = sqrt((constants.kb * constants.T) / (constants.u0 * constants.Ms * constants.Heff * constants.v))
    phi = pi
    mz = cos(theta)

    F = (constants.tox / (constants.RA * (constants.phi_bar ** 0.5))) * exp((2 * constants.tox * (2 * constants.m * constants.e * constants.phi_bar) ** 0.5) / constants.h_bar)
    Rp = (constants.tox / (F * (constants.phi_bar ** 0.5) * constants.A1)) * exp((2 * constants.tox * (2 * constants.m * constants.e * constants.phi_bar) ** 0.5) / constants.h_bar)
    R_MTJ = Rp * (1 + constants.TMR / (constants.TMR + 2)) / (1 + constants.TMR * cos(theta) / (constants.TMR + 2))
    return (R_MTJ, theta, mz, phi)
