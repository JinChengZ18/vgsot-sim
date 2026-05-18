import numpy as np


def stochastic(n=1, rng=None):
    """Three independent unit Gaussian samples for the thermal field.

    Brown 1963 / FDT: each Cartesian component of H_TH is an independent
    N(0, 1) draw, scaled afterwards by sqrt(2 alpha kT / (mu0 Ms gamma V dt))
    in the caller. Returning a unit-length vector (previous behaviour) gives
    constant |H_TH| and underestimates the noise variance.

    Parameters
    ----------
    n : float
        Standard deviation passed to np.random.normal. Kept for backward
        compatibility; the physical scaling is applied by the caller, so the
        default n=1 is the only correct choice. Non-default values are kept
        only to avoid breaking older callers and will produce non-physical
        statistics.
    rng : np.random.Generator, optional
        If provided, used as the random source (enables seeding for
        reproducibility). Falls back to legacy np.random global state when
        None to preserve prior call sites.

    Returns
    -------
    np.ndarray, shape (3,)
        [xi_x, xi_y, xi_z], each ~ N(0, n^2), statistically independent.
    """
    if rng is None:
        return np.random.normal(0.0, n, 3)
    return rng.normal(0.0, n, 3)
