from math import sqrt, exp, cos, pi

import numpy as np

from .configs import PhysicalConstantsConfig


def compute_Rp(constants: PhysicalConstantsConfig):
    """Brinkman-Dynes-Rowell parallel-state resistance.

    Pure function of barrier and area parameters; no random component.
    Pulled out of init() so callers that only need R_P (e.g. tmr.py) do
    not pay the cost of an initial-angle draw and do not pick up its
    stochasticity in their resistance value.
    """
    F = (
        constants.tox
        / (constants.RA * (constants.phi_bar ** 0.5))
    ) * exp(
        (2 * constants.tox
         * (2 * constants.m * constants.e * constants.phi_bar) ** 0.5)
        / constants.h_bar
    )
    Rp = (
        constants.tox
        / (F * (constants.phi_bar ** 0.5) * constants.A1)
    ) * exp(
        (2 * constants.tox
         * (2 * constants.m * constants.e * constants.phi_bar) ** 0.5)
        / constants.h_bar
    )
    return Rp


def init(PAP, constants: PhysicalConstantsConfig, rng=None):
    """
    To initialize the VGSOT MTJ in a thermally equilibrated state near one
    of the two perpendicular minima.

    Near each pole, the polar angle θ has a Rayleigh distribution
        p(θ) ∝ θ · exp(-θ² / (2 σ²)),   σ² ≈ k_B T / (μ_0 M_s H_eff V) ≈ 1/(2Δ)
    while the azimuth φ is uniform on [0, 2π). Drawing each MC trial from
    this distribution (rather than from a single deterministic offset, as in
    the upstream port) gives the Monte-Carlo ensemble its proper thermal
    spread; using a fixed initial point biases short-pulse switching
    probabilities by suppressing the natural starting-angle variance.

    Parameters
    ----------
    PAP : int
        1 for the "parallel-pole" basin (θ near π), 0 for the other basin.
    constants : PhysicalConstantsConfig
    rng : np.random.Generator, optional
        Seedable RNG for reproducibility. Falls back to np.random global.

    Returns
    -------
    R_MTJ, theta, mz, phi
    """
    if rng is None:
        rng = np.random

    sigma_theta = sqrt(
        constants.kb * constants.T
        / (constants.u0 * constants.Ms * constants.Heff * constants.v)
    )
    # Rayleigh draw via two independent Gaussians; equivalent to
    # np.random.rayleigh(sigma) but exposes the rng argument uniformly.
    u = rng.normal(0.0, sigma_theta)
    v = rng.normal(0.0, sigma_theta)
    delta_theta = sqrt(u * u + v * v)
    phi = float(np.arctan2(v, u) % (2 * pi))

    if PAP == 1:
        theta = pi - delta_theta
    else:
        theta = delta_theta
    mz = cos(theta)

    Rp = compute_Rp(constants)
    R_MTJ = Rp * (1 + constants.TMR / (constants.TMR + 2)) / (1 + constants.TMR * cos(theta) / (constants.TMR + 2))
    return (R_MTJ, theta, mz, phi)
