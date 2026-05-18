from math import sqrt, sin, cos, pi
import numpy as np

from .configs import PhysicalConstantsConfig
from .stochastic import stochastic
from .demag import demag_factors


def field(theta, phi, V_MTJ, n, NON, ENE, VNV, constants: PhysicalConstantsConfig,
          demag_mode: str = "ellipsoid", Ki_T: float | None = None, Ms_T: float | None = None,
          rng=None):
    """
    Anisotropy module
    H_eff = H_PMA + H_D + H_TH + H_EX + H_VCMA, vectorwise

    Optional `Ki_T` and `Ms_T` overrides allow callers (e.g. the self-heating
    loop) to inject already-temperature-corrected material parameters. When
    omitted, the RT values from `constants` are used.

    `demag_mode` selects the demag tensor model — "ellipsoid" (default, exact
    oblate spheroid) or "thin_disk" (linearised, legacy). Both use D_elec
    for the in-plane diameter.

    `rng` (optional `np.random.Generator`) is forwarded to `stochastic()` so
    callers can fully control the thermal-noise stream for byte-reproducible
    simulations. When None (default) the legacy global np.random state is used.
    """
    ex = np.array([1, 0, 0])
    ey = np.array([0, 1, 0])
    ez = np.array([0, 0, 1])

    mx = sin(theta)*cos(phi)*ex
    my = sin(theta)*sin(phi)*ey
    mz = cos(theta)*ez
    m = np.add.reduce([mx, my, mz])

    Ki_use = constants.Ki if Ki_T is None else Ki_T
    Ms_use = constants.Ms if Ms_T is None else Ms_T

    H_PMA = 2 * Ki_use / (constants.u0 * Ms_use * constants.tf) * mz
    H_VCMA = -2 * constants.beta * V_MTJ / (constants.u0 * Ms_use * constants.tox * constants.tf) * mz

    Nx, Ny, Nz = demag_factors(constants, mode=demag_mode)
    N = np.array([Nx, 0, 0, 0, Ny, 0, 0, 0, Nz]).reshape(3, 3)
    H_D = -Ms_use * np.dot(m, N)

    H_EX = np.array([constants.h_ex_x, constants.h_ex_y, constants.h_ex_z])

    H_th_mag = sqrt(2 * constants.kb * constants.T * constants.alpha / (constants.u0 * Ms_use * constants.gamma * constants.v * constants.t_step))
    xi = stochastic(n, rng=rng)
    H_TH = H_th_mag * xi

    H_eff = H_PMA + H_D + NON * H_TH + ENE * H_EX + VNV * H_VCMA
    H_eff_perpendicular = np.dot(H_PMA + VNV * H_VCMA + H_D, ez)
    return H_eff, H_eff_perpendicular
