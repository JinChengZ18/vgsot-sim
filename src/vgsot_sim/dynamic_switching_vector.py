"""
Vector-form sLLG step with **Cayley-transform** rotation update.

This module is the alternative integrator promised by §2.2.3.2 of the
chapter. It complements the legacy spherical-coordinate explicit-Euler
stepper in `dynamic_switching.py` with two improvements:

1. **Explicit spin-Hall polarisation `sigma_SH`** instead of the LLG being
   re-derived for one particular σ̂ direction. The user passes
   `sigma_SH` as a unit 3-vector; the same routine handles σ̂ = ±x̂, ±ŷ,
   or any oblique direction.

2. **Cayley midpoint rotation** that conserves `|m|=1` to machine
   precision. Writing the LL form as
   $$\\dot{\\mathbf{m}} = \\boldsymbol{\\omega}(\\mathbf{m}) \\times \\mathbf{m},$$
   the Cayley step
   $$\\mathbf{m}_{n+1} = (I - \\tfrac{\\Delta t}{2}[\\boldsymbol{\\omega}]_\\times)^{-1}(I + \\tfrac{\\Delta t}{2}[\\boldsymbol{\\omega}]_\\times)\\mathbf{m}_n$$
   is unconditionally norm-preserving and second-order accurate when
   `omega` is evaluated at the midpoint $\\mathbf{m}_{n+1/2}$.

The closed-form Cayley rotation for a 3-vector `omega` is
   $$\\mathbf{m}_{n+1} = \\mathbf{m}_n + \\frac{2s}{1 + s^2 |\\omega|^2}
                     \\Big(\\boldsymbol{\\omega}\\times\\mathbf{m}_n
                          + s\\,\\boldsymbol{\\omega}\\times(\\boldsymbol{\\omega}\\times\\mathbf{m}_n)\\Big),
   \\quad s \\equiv \\Delta t/2.$$
which avoids any matrix inversion at runtime.

For backward compatibility the spherical-Euler stepper in
`dynamic_switching.switching()` is unchanged. Choose between the two via
`run_piecewise_direct_excitation(..., integrator="cayley")` or
`= "euler_spherical"` (default).
"""
from __future__ import annotations

import numpy as np

from .anisotropy import field
from .configs import PhysicalConstantsConfig


def _omega_from_state(m, H_eff, sigma_SH, *,
                      H_DL_SOT, H_FL_SOT, H_DL_STT, H_FL_STT, sigma_STT,
                      alpha, gamma_red):
    """Return the angular-velocity vector ω such that dm/dt = ω × m.

    Derivation. The standard LL form is

        dm/dt = -γ_red m × H_eff − α γ_red m × (m × H_eff)
                + SOT_contribution + STT_contribution.

    Both H_eff terms and the SOT/STT torque terms are perpendicular to m
    (as required for |m|=1). For any D ⊥ m and |m|=1, taking ω = m × D
    yields ω × m = D (BAC–CAB identity). Applying this to each component:

      precession + damping:
          D_H = -γ_red m × H − αγ_red m × (m × H)
          ω_H = m × D_H = +γ_red H + αγ_red (m × H)        (modulo m-parallel)

      SOT/STT (both Gilbert torques after LL conversion):
          D_J = (α H_DL + H_FL) γ_red (m × σ̂)
               + (α H_FL − H_DL) γ_red (m × (m × σ̂))
          ω_J = m × D_J          (computed directly)

    The m-parallel components of ω are immaterial — they contribute nothing
    to ω × m — so we drop them where convenient.

    Returns a (3,) ndarray.
    """
    # Effective-field part (precession + damping). After ω = m × D_H and
    # dropping the m-parallel component, the contribution simplifies to
    # +γ_red H + αγ_red (m × H).
    H = np.asarray(H_eff, dtype=float)
    m = np.asarray(m, dtype=float)
    omega = gamma_red * H + alpha * gamma_red * np.cross(m, H)

    # SOT torques. The standard Gilbert form is
    #     tau_DL = -gamma * H_DL * m x (m x sigma)
    #     tau_FL = +gamma * H_FL * m x sigma
    # which after the LL conversion (1+alpha²) factor enters as
    #     dm/dt += (alpha H_DL + H_FL) gamma_red (m x sigma)
    #          + (alpha H_FL - H_DL) gamma_red (m x (m x sigma))
    # The first term is a torque a × m with a = -(alpha H_DL + H_FL) gamma_red sigma.
    # The second term is m × (m × A) = (m·A) m - A; this is NOT of form (a × m), but its
    # action on m is identical to a precession by ω = (-(alpha H_FL - H_DL) gamma_red) A_perp
    # where A_perp = A - (m·A) m. For unit m, (m × (m × A)) = -A_perp.
    # Since we want dm/dt = ω × m, and (ω × m) is automatically perpendicular to m, we can
    # absorb the parallel component freely. We use the equivalence
    #     m × (m × A) = -A_perp   =>   -A_perp × m / |m|² × m  but |m|=1
    # Easier: just compute ω so that ω × m matches the desired dm/dt.
    # Given dm/dt = D where D ⊥ m, we have ω = m × D + λm (parallel component arbitrary).
    # We pick the minimal-magnitude ω = m × D (perpendicular to m).
    #
    # Build D = full LL RHS, then take ω = m × D - (precession & damping contribute
    # naturally with the above derivation).
    if H_DL_SOT != 0.0 or H_FL_SOT != 0.0:
        sigma_S = np.asarray(sigma_SH, dtype=float)
        # (m × sigma) and (m × (m × sigma))
        mxs = np.cross(m, sigma_S)
        mxmxs = np.cross(m, mxs)
        # SOT contribution to dm/dt
        D_sot = (alpha * H_DL_SOT + H_FL_SOT) * gamma_red * mxs \
              + (alpha * H_FL_SOT - H_DL_SOT) * gamma_red * mxmxs
        # Convert to omega-contribution using ω_sot = m × D_sot (since D ⊥ m)
        omega = omega + np.cross(m, D_sot)

    if H_DL_STT != 0.0 or H_FL_STT != 0.0:
        sigma_T = np.asarray(sigma_STT, dtype=float)
        mxs = np.cross(m, sigma_T)
        mxmxs = np.cross(m, mxs)
        D_stt = (alpha * H_DL_STT + H_FL_STT) * gamma_red * mxs \
              + (alpha * H_FL_STT - H_DL_STT) * gamma_red * mxmxs
        omega = omega + np.cross(m, D_stt)

    return omega


def cayley_step(m, omega, dt):
    """Rotate `m` by the Cayley transform of (dt × [ω]_×).

    Returns the new magnetisation, preserving |m| to machine precision.
    Closed-form expression — no matrix inversion at runtime:

        m_new = m + (2s / (1 + s² |ω|²)) (ω × m + s ω × (ω × m))

    where s = dt/2.
    """
    s = 0.5 * dt
    cross1 = np.cross(omega, m)
    cross2 = np.cross(omega, cross1)
    omega_sq = float(np.dot(omega, omega))
    factor = (2.0 * s) / (1.0 + s * s * omega_sq)
    return m + factor * (cross1 + s * cross2)


def switching_vector(m, V_MTJ, I_SOT, R_MTJ, ESTT, ESOT, *,
                     VNV=1, NON=0,
                     R_SOT_FL_DL=0.83, R_STT_FL_DL=0.0,
                     sigma_SH=np.array([-1.0, 0.0, 0.0]),
                     sigma_STT=np.array([0.0, 0.0, +1.0]),
                     constants: PhysicalConstantsConfig | None = None,
                     Ki_T: float | None = None,
                     Ms_T: float | None = None,
                     demag_mode: str = "ellipsoid",
                     rng=None):
    """One Cayley-transform sLLG step in Cartesian form.

    Computes a single-step Cayley rotation of the unit magnetisation
    `m` (3-vector) under the combined precession + damping + SOT + STT
    dynamics, with the spin polarisation directions exposed via
    `sigma_SH` and `sigma_STT`.

    Compared to the spherical-Euler stepper in
    `dynamic_switching.switching()`:
      * `|m|=1` is preserved to ~1e-16 (no manual renormalisation needed)
      * the σ̂_SH direction is a runtime parameter rather than implicit
        in the closed-form expansion
      * no 1/sin(θ) pole near the easy axis

    `rng` (optional `np.random.Generator`) is forwarded to `field()` →
    `stochastic()` so the thermal-noise stream is byte-reproducible without
    touching `np.random`'s global state. When None (default) the legacy
    global state is used.

    Returns
    -------
    m_new : np.ndarray, shape (3,)
        Updated magnetisation. Unit length to machine precision.
    """
    if constants is None:
        raise ValueError("constants must be provided.")

    m = np.asarray(m, dtype=float)
    # Spherical coords are needed only to call field() (which is parameterised
    # by θ, φ). Use the current m to derive them.
    theta = float(np.arccos(np.clip(m[2], -1.0, 1.0)))
    phi = float(np.arctan2(m[1], m[0]))

    H_eff, _ = field(theta, phi, V_MTJ, n=1, NON=NON, ENE=1, VNV=VNV,
                     constants=constants, demag_mode=demag_mode,
                     Ki_T=Ki_T, Ms_T=Ms_T, rng=rng)
    H_eff = np.asarray(H_eff, dtype=float)

    # T-corrected M_s flows into the SOT/STT prefactor (matches the
    # spherical stepper).
    Ms_use = constants.Ms if Ms_T is None else Ms_T

    I_MTJ = V_MTJ / R_MTJ if R_MTJ != 0.0 else 0.0
    J_STT = I_MTJ / constants.A1
    J_SOT = I_SOT / constants.A2

    gamma_red = constants.gamma / (1.0 + constants.alpha ** 2)
    H_DL_STT = ESTT * constants.h_bar * constants.P * J_STT / (2 * constants.e * constants.u0 * Ms_use * constants.tf)
    H_FL_STT = R_STT_FL_DL * H_DL_STT
    H_DL_SOT = ESOT * constants.h_bar * constants.theta_SH * J_SOT / (2 * constants.e * constants.u0 * Ms_use * constants.tf)
    H_FL_SOT = R_SOT_FL_DL * H_DL_SOT

    omega = _omega_from_state(
        m, H_eff, sigma_SH,
        H_DL_SOT=H_DL_SOT, H_FL_SOT=H_FL_SOT,
        H_DL_STT=H_DL_STT, H_FL_STT=H_FL_STT,
        sigma_STT=sigma_STT,
        alpha=constants.alpha, gamma_red=gamma_red,
    )
    m_new = cayley_step(m, omega, constants.t_step)
    # Renormalise to remove the residual O(machine eps) drift.
    nrm = np.linalg.norm(m_new)
    if nrm > 0.0:
        m_new = m_new / nrm
    return m_new
