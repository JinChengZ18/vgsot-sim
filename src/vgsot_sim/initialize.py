from math import sqrt, exp, cos, pi

import numpy as np

from .configs import PhysicalConstantsConfig


def resistance_area_bdr(constants: PhysicalConstantsConfig, m_eff_ratio: float = 1.0) -> float:
    """Simmons / Brinkman-Dynes-Rowell low-bias tunnel resistance-area product.

    Genuine barrier-physics PREDICTOR of RA from the rectangular-barrier
    Simmons (1963) low-voltage limit:

        RA(phi, t_ox) = (h^2 t_ox) / (q^2 sqrt(2 m* e phi))
                        * exp( 2 t_ox sqrt(2 m* e phi) / hbar )    [Ohm.m^2]

    so RA grows with BOTH the barrier height ``phi_bar`` (eV) and the barrier
    thickness ``tox`` — unlike the legacy expression this is genuinely a function
    of the barrier (it does not cancel to RA/A1). ``m_eff_ratio`` scales the
    free-electron mass to an effective mass (≈0.3–0.4 for MgO brings the nominal
    phi_bar=0.4 eV prediction down to the measured RA range).

    This is a PREDICTOR, not the model's working RA: the calibrated compact model
    uses the *measured* ``constants.RA`` (see :func:`compute_Rp`). At the nominal
    free-electron mass + phi_bar=0.4 eV this predicts RA ~10x the measured value,
    which is exactly why a measured/effective-mass-fitted RA is used directly
    rather than a first-principles barrier output.
    """
    h = 2.0 * pi * constants.h_bar
    m = constants.m * m_eff_ratio
    s = (2.0 * m * constants.e * constants.phi_bar) ** 0.5          # sqrt(2 m* e phi)  [kg.m/s]
    return (h ** 2 * constants.tox) / (constants.e ** 2 * s) \
        * exp(2.0 * constants.tox * s / constants.h_bar)


def compute_Rp(constants: PhysicalConstantsConfig) -> float:
    """Parallel-state MTJ resistance R_P = RA / A1.

    R_P is the calibrated resistance-area product ``constants.RA`` divided by the
    electrical area ``A1 = pi D_elec^2 / 4``. RA is the measured calibration
    handle (default 16.6e-12 Ohm.m^2 -> R_P ~ 5 kOhm at D_elec = 65 nm). The
    barrier-physics basis for RA (Simmons / Brinkman-Dynes-Rowell, where the
    barrier height and thickness genuinely enter) is provided by
    :func:`resistance_area_bdr`, a separate predictor.

    Pure function of RA and area; no random component, so callers needing only
    R_P (e.g. ``tmr.py``) do not pay an initial-angle draw or pick up its
    stochasticity.

    NOTE (corrected 2026): the previous body dressed R_P = RA/A1 in a
    BDR-looking expression whose barrier-tunneling factors (the WKB exponential,
    t_ox and sqrt(phi)) cancelled EXACTLY, so the result never depended on the
    barrier. It now returns the honest closed form, with the real BDR physics
    available via :func:`resistance_area_bdr`.
    """
    return constants.RA / constants.A1


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
