"""Demagnetisation tensor for the free-layer MTJ pillar.

Two forms are exposed (selected by `mode` in `demag_factors`):

- "thin_disk":  N_x = N_y = π t_f / (4 D),  N_z = 1 - 2 N_x
                Linearisation valid only when q = t_f / D << 1.

- "ellipsoid": exact oblate-spheroid analytic solution[*]
                N_x = N_y = q² / (2(1-q²)) · [arccos(q) / (q√(1-q²)) - 1]
                N_z = 1 - 2 N_x
                Reduces to the thin-disk form in the q → 0 limit.

The exact form is the default; switching to "thin_disk" is provided only
for backwards comparison with prior runs and with documents that quote the
linearised expression.

[*] Stoner–Wohlfarth oblate spheroid; see thesis §2.2.2.4 and the standard
reference: J. A. Osborn, Phys. Rev. 67, 351 (1945).
"""
from __future__ import annotations

from math import acos, pi, sqrt

from .configs import PhysicalConstantsConfig


def demag_factors(constants: PhysicalConstantsConfig, mode: str = "ellipsoid"):
    """Return (N_x, N_y, N_z) for the free-layer pillar.

    Uses constants.tf for the short axis and constants.D_elec (NOT D_phys)
    for the long-axis diameter, matching the magnetic-active region.

    Parameters
    ----------
    constants : PhysicalConstantsConfig
    mode : {"ellipsoid", "thin_disk"}
        Demag tensor model. Defaults to the exact ellipsoid form.

    Returns
    -------
    (Nx, Ny, Nz) : tuple of floats
        Diagonal entries of the demag tensor; N_x + N_y + N_z = 1 for any
        ellipsoid with the symmetry t_x = t_y < t_z (or here the in-plane
        diameter > out-of-plane thickness, so N_x = N_y < N_z is enforced).

    Notes
    -----
    For a perpendicular-MTJ free layer the z axis points out-of-plane.
    Because the geometry is oblate (in-plane diameter > thickness), the
    demag factor along the easy axis is larger than along the in-plane
    directions, penalising out-of-plane orientation by 0.5 μ_0 M_s²(N_z-N_x).
    """
    q = constants.tf / constants.D_elec
    if not 0 < q < 1:
        raise ValueError(
            f"Aspect ratio q = t_f/D_elec = {q!r} is outside (0, 1); "
            "ellipsoid demag is only valid for oblate pillars."
        )

    if mode == "thin_disk":
        Nx = pi * q / 4
    elif mode == "ellipsoid":
        # Oblate spheroid analytic result (Osborn 1945; thesis §2.2.2.4).
        sqrt_1_q2 = sqrt(1.0 - q * q)
        Nx = (q * q / (2.0 * (1.0 - q * q))) * (
            acos(q) / (q * sqrt_1_q2) - 1.0
        )
    else:
        raise ValueError(f"Unknown demag mode={mode!r}. Use 'ellipsoid' or 'thin_disk'.")
    Ny = Nx
    Nz = 1.0 - 2.0 * Nx
    return (Nx, Ny, Nz)
