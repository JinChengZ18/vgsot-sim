from math import sqrt, sin, cos, pi
import numpy as np

from .configs import PhysicalConstantsConfig
from .stochastic import stochastic


def field(theta, phi, V_MTJ, n, NON, ENE, VNV, constants: PhysicalConstantsConfig):
    """
    Anisotropy module
    H_eff = H_PMA + H_D + H_TH + H_EX + H_VCMA, vectorwise
    """
    ex = np.array([1, 0, 0])
    ey = np.array([0, 1, 0])
    ez = np.array([0, 0, 1])

    mx = sin(theta)*cos(phi)*ex
    my = sin(theta)*sin(phi)*ey
    mz = cos(theta)*ez
    m = np.add.reduce([mx, my, mz])

    H_PMA = 2 * constants.Ki / (constants.u0 * constants.Ms * constants.tf) * mz
    H_VCMA = -2 * constants.beta * V_MTJ / (constants.u0 * constants.Ms * constants.tox * constants.tf) * mz

    Nx = Ny = pi * constants.tf / (4 * constants.D)
    Nz = 1 - 2 * Nx
    N = np.array([Nx, 0, 0, 0, Ny, 0, 0, 0, Nz]).reshape(3, 3)
    H_D = -constants.Ms * np.dot(m, N)

    Hx_ex = 0
    Hy_ex = -50 * 1000 / (4 * pi)
    Hz_ex = 0
    H_EX = Hx_ex*ex + Hy_ex*ey + Hz_ex*ez

    H_th_mag = sqrt(2 * constants.kb * constants.T * constants.alpha / (constants.u0 * constants.Ms * constants.gamma * constants.v * constants.t_step))
    sigma = stochastic(n)
    H_TH = np.dot(H_th_mag, sigma)

    H_eff = H_PMA + H_D + NON * H_TH + ENE * H_EX + VNV * H_VCMA
    H_eff_perpendicular = np.dot(H_PMA + VNV * H_VCMA + H_D, ez)
    return H_eff, H_eff_perpendicular
