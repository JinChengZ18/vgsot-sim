"""Néel-Brown thermal-activation fitting from measured hysteresis loops.

Given a set of (pulse-width, R-V loop) traces, extract critical voltages,
fit a log-linear V_c(t_w) trend, and invert it into the underlying
thermal-stability factor Δ and zero-temperature critical voltage V_c0.

Used by the 06 / 07 analysis scripts; consolidated here so the same fit
can be invoked from notebooks or other vgsot-sim cases.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def extract_Vc(V, R, V_guard: float = 0.05):
    """Return (V_th_minus, V_th_plus) from an R-V hysteresis trace.

    Definition: V_th_+ is the smallest positive V > V_guard at which R
    drops below the midpoint R_th = (R_min + R_max)/2 when sweeping
    backwards from high V; V_th_- is the analogous on the negative side.
    """
    V = np.asarray(V, dtype=float)
    R = np.asarray(R, dtype=float)
    R_th = 0.5 * (R.min() + R.max())
    neg = sorted([(v, r) for v, r in zip(V, R) if v < -V_guard], key=lambda x: x[0])
    V_th_minus = next((v for v, r in neg if r > R_th), np.nan)
    pos = sorted([(v, r) for v, r in zip(V, R) if v > V_guard], key=lambda x: -x[0])
    V_th_plus = next((v for v, r in pos if r < R_th), np.nan)
    return V_th_minus, V_th_plus


def loglin_fit(tw_arr, v_arr):
    """Fit V(t_w) = a - b·ln(t_w). Returns (a, b)."""
    tw_arr = np.asarray(tw_arr, dtype=float)
    v_arr = np.asarray(v_arr, dtype=float)
    slope, intercept = np.polyfit(np.log(tw_arr), v_arr, 1)
    return float(intercept), float(-slope)


def nb_params(a, b, tau0: float = 1.0):
    """Invert (a, b) into (Δ, V_c0). Assumes V_c0 = a - b·ln(τ_0·ln 2).

    With V(t_w) = V_c0 (1 - ln(t_w/(τ_0 ln 2))/Δ), the fit slope b equals
    V_c0/Δ. Returns (Δ, V_c0).
    """
    Vc0 = a - b * np.log(tau0 * np.log(2))
    return float(Vc0 / b), float(Vc0)


def tau_ret(tau0: float, Delta: float) -> float:
    """Néel-Brown retention time τ_ret = τ_0 · exp(Δ)."""
    return float(tau0 * np.exp(Delta))


def psw_nb(V, t_ns, Delta: float, Vc0: float, tau0: float = 1.0):
    """Néel-Brown switching probability, with V > V_c0 clipped to deterministic."""
    V = np.asarray(V, dtype=float)
    t_ns = np.asarray(t_ns, dtype=float)
    x = np.maximum(Delta * (1.0 - V / Vc0), 0.0)
    return 1.0 - np.exp(-(t_ns / tau0) * np.exp(-x))


def vth_nb(t_ns, Delta: float, Vc0: float, tau0: float = 1.0):
    """Median-switching voltage at given pulse width."""
    t_ns = np.asarray(t_ns, dtype=float)
    arg = np.maximum(t_ns / (tau0 * np.log(2.0)), 1.0 + 1e-14)
    return Vc0 * (1.0 - np.log(arg) / Delta)


def beta_nb_analytic(Delta: float, Vc0: float) -> float:
    """Analytic Sigmoid slope of the NB curve at its median: 2 Δ ln2 / V_c0 (V⁻¹)."""
    return float(2.0 * Delta * np.log(2.0) / Vc0)


@dataclass
class NBFitResult:
    """Container for the inverted Néel-Brown parameters of one direction.

    Attributes
    ----------
    a, b : float
        Log-linear coefficients V(t) = a - b·ln(t/ns).
    Delta : float
        Thermal-stability factor Δ = E_b / (k_B T).
    Vc0 : float
        Zero-temperature critical voltage (V).
    tau_ret_ns : float
        Retention time τ_0 · exp(Δ) at the assumed τ_0 (ns).
    """
    a: float
    b: float
    Delta: float
    Vc0: float
    tau_ret_ns: float


def fit_direction(tw_arr, v_arr, tau0: float = 1.0) -> NBFitResult:
    """Fit log-linear NB law to one (pulse width, critical voltage) series.

    Parameters
    ----------
    tw_arr : array-like
        Pulse widths in ns.
    v_arr : array-like
        Corresponding critical voltages in V (absolute magnitude).
    tau0 : float, default 1.0
        Attempt-time inverse, in ns.
    """
    a, b = loglin_fit(tw_arr, v_arr)
    Delta, Vc0 = nb_params(a, b, tau0=tau0)
    return NBFitResult(a=a, b=b, Delta=Delta, Vc0=Vc0,
                       tau_ret_ns=tau_ret(tau0, Delta))
