"""Four-parameter Sigmoid fit + η_c calibration vs Néel-Brown.

Given a measured (V, P_sw) trace (one device, one direction, one pulse
width), fit a 4-parameter logistic
    P_sw(V) = y0 + L / (1 + exp(-(V - V_th)/k))
and report the slope β_s = 1/k and the goodness R². The C2C calibration
factor η_c is computed against the analytic NB slope at the same
(Δ, V_c0, t_w, τ_0) operating point.

Used by 07 sigmoid_fig and downstream variability_sim.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import curve_fit

from . import nb_fit


def sigmoid4p(V, y0, L, Vth, k):
    """4-parameter logistic. `k` is the scale (mV per e-fold)."""
    return y0 + L / (1.0 + np.exp(-(V - Vth) / k))


def wilson(p, n: int = 100, z: float = 1.96):
    """Wilson 95% confidence interval for a binomial proportion.

    Vectorised over p. Handles edge cases p ∈ {0, 1} cleanly (interval
    degenerates to the point estimate).
    """
    p = np.clip(np.asarray(p, dtype=float), 0.0, 1.0)
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return centre - half, centre + half


@dataclass
class SigmoidFitResult:
    """One-curve Sigmoid fit summary.

    Attributes
    ----------
    y0, L, Vth, k : float
        Logistic parameters in fit units (V).
    beta : float
        Slope 1/k (V⁻¹).
    R2 : float
        Coefficient of determination of the fit residuals.
    """
    y0: float
    L: float
    Vth: float
    k: float
    beta: float
    R2: float


def fit_sigmoid(V, P, p0=None, bounds=None) -> SigmoidFitResult:
    """Fit `sigmoid4p` to (V, P) measurements via least squares."""
    V = np.asarray(V, dtype=float)
    P = np.asarray(P, dtype=float)
    if p0 is None:
        p0 = (0.0, 1.0, float(np.median(V)), 0.02)
    if bounds is None:
        bounds = ([-0.2, 0.2, V.min() - 0.2, 1e-4],
                  [0.2, 1.5, V.max() + 0.2, 0.5])
    popt, _ = curve_fit(sigmoid4p, V, P, p0=p0, bounds=bounds, maxfev=50000)
    y0, L, Vth, k = popt
    pred = sigmoid4p(V, *popt)
    ss_res = float(np.sum((P - pred) ** 2))
    ss_tot = float(np.sum((P - np.mean(P)) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return SigmoidFitResult(
        y0=float(y0), L=float(L), Vth=float(Vth), k=float(k),
        beta=float(1.0 / k), R2=R2,
    )


def eta_c(beta_meas: float, beta_nb: float) -> float:
    """C2C calibration factor η_c = β_meas / β_NB."""
    return float(beta_meas / beta_nb)


@dataclass
class NBReference:
    """Bundles the NB reference at a single operating point.

    Use both the analytic slope (closed form) and the curve-fit slope
    (Sigmoid fit to the NB curve over a V range), because they differ
    when the NB curve has the Gumbel skew.
    """
    Delta: float
    Vc0: float
    tau0_ns: float
    tw_ns: float
    Vth_NB: float
    beta_NB_analytic: float


def nb_reference(Delta: float, Vc0: float, tw_ns: float, tau0_ns: float = 1.0) -> NBReference:
    """Compute the NB threshold and analytic slope at a pulse width."""
    Vth_NB = float(nb_fit.vth_nb(tw_ns, Delta, Vc0, tau0=tau0_ns))
    beta_NB_a = nb_fit.beta_nb_analytic(Delta, Vc0)
    return NBReference(Delta=Delta, Vc0=Vc0, tau0_ns=tau0_ns, tw_ns=tw_ns,
                       Vth_NB=Vth_NB, beta_NB_analytic=beta_NB_a)
