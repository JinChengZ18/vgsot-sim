"""Process-variability budget + D2D Monte-Carlo for sMTJ arrays.

Combines the Brinkman-decomposed CV(Δ) budget (from PDK mismatch params)
with a Monte-Carlo sweep over CV(Δ) values, fitting a Sigmoid to each
wafer-averaged P_sw(V) curve and returning the broadening factor F.

Used by 07/variability_sim-V3.0.py; the analytic skeleton lives here so
the same machinery can be invoked from a vgsot-sim case or notebook.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import curve_fit

from . import nb_fit
from .sigmoid_fit import sigmoid4p


@dataclass
class PDKBudgetInputs:
    """PDK-mismatch CV inputs to the Brinkman variance budget."""
    CV_RP: float = 0.07          # CV(R_P) from PDK Rp_mis
    CV_RSOT: float = 0.07        # CV(R_SOT) from PDK Rsot_mis  (informational only)
    CV_TMR: float = 0.04         # CV(TMR) used as proxy for CV(H_k)
    CV_TOX: float = 0.003        # CV(t_ox) — 300 mm process control
    CV_PHI: float = 0.005        # CV(phi_ox)
    CV_TF:  float = 0.003        # CV(t_f) — MBE precision
    CV_MS:  float = 0.02         # CV(M_s) — literature
    T_OX:   float = 1.0          # nominal t_ox (nm)
    PHI_BAR: float = 0.6         # nominal phi (eV)
    KAPPA:  float = 10.25        # BDR sensitivity: d ln(R·A)/d t_ox per √eV


@dataclass
class CVBudgetResult:
    """Decomposed CV(Δ) budget; component CVs all dimensionless."""
    CV_RA: float
    CV_D:  float
    CV_V:  float       # CV of magnetic volume
    CV_HK: float
    CV_MS: float
    CV_TF: float
    CV_Delta: float
    components: dict = field(default_factory=dict)


def cv_delta_budget(inp: PDKBudgetInputs | None = None) -> CVBudgetResult:
    """Brinkman-decomposed CV(Δ) from PDK mismatch parameters.

    The decomposition assumes:
      ln(R·A) ≈ κ t_ox √φ + const  (BDR low-bias)
    so that CV(R·A) follows from CV(t_ox) and CV(φ); the residual
    variance in CV(R_P) is attributed to area mismatch -> CV(D).
    Δ ∝ K_eff V / kT scales with H_k, M_s and volume V.
    """
    if inp is None:
        inp = PDKBudgetInputs()
    SENS_TOX = inp.KAPPA * np.sqrt(inp.PHI_BAR)
    SENS_PHI = inp.KAPPA * inp.T_OX / (2 * np.sqrt(inp.PHI_BAR))
    CV_RA = float(np.sqrt(
        (SENS_TOX * inp.CV_TOX * inp.T_OX) ** 2
        + (SENS_PHI * inp.CV_PHI * inp.PHI_BAR) ** 2
    ))
    CV_D  = float(np.sqrt(max(inp.CV_RP ** 2 - CV_RA ** 2, 0.0)) / 2.0)
    CV_V  = float(np.sqrt((2 * CV_D) ** 2 + inp.CV_TF ** 2))
    CV_HK = inp.CV_TMR
    CV_Delta = float(np.sqrt(CV_HK ** 2 + inp.CV_MS ** 2 + CV_V ** 2))
    return CVBudgetResult(
        CV_RA=CV_RA, CV_D=CV_D, CV_V=CV_V,
        CV_HK=CV_HK, CV_MS=inp.CV_MS, CV_TF=inp.CV_TF,
        CV_Delta=CV_Delta,
        components={
            "V_mag":    CV_V ** 2 / CV_Delta ** 2,
            "H_k":      CV_HK ** 2 / CV_Delta ** 2,
            "M_s":      inp.CV_MS ** 2 / CV_Delta ** 2,
            "t_f":      inp.CV_TF ** 2 / CV_Delta ** 2,
        },
    )


def wafer_average_psw(V_dense, cv: float, Delta: float, Vc0: float,
                      tw_ns: float, tau0_ns: float = 1.0,
                      N_samples: int = 20000, rng=None):
    """Monte-Carlo wafer-averaged P_sw(V) for a Gaussian Δ distribution.

    For cv == 0 returns the deterministic single-device curve; otherwise
    draws N_samples Δ ~ N(Delta, cv*Delta), clipped at Δ_min = 0.5 to
    avoid negative arguments.
    """
    V_dense = np.asarray(V_dense, dtype=float)
    if rng is None:
        rng = np.random.default_rng()
    if cv == 0.0:
        Ds = np.array([Delta])
    else:
        Ds = rng.normal(Delta, cv * Delta, N_samples)
        Ds = np.clip(Ds, 0.5, None)
    P_mean = np.zeros_like(V_dense)
    for j, V in enumerate(V_dense):
        x = Ds * (1.0 - V / Vc0)
        P_mean[j] = (1.0 - np.exp(-(tw_ns / tau0_ns) * np.exp(-x))).mean()
    return P_mean


def fit_sigmoid_to_average(V_dense, P_mean, p0=None, bounds=None):
    """Fit `1/(1+exp(-β(V-Vth)))` to wafer-averaged curve. Returns (Vth, β)."""
    V_dense = np.asarray(V_dense, dtype=float)
    P_mean = np.asarray(P_mean, dtype=float)

    def sigmoid(V, Vth, beta):
        return 1.0 / (1.0 + np.exp(-beta * (V - Vth)))

    if p0 is None:
        p0 = (float(np.median(V_dense)), 50.0)
    if bounds is None:
        bounds = ([0.0, 0.1], [1.5, 1000.0])
    popt, _ = curve_fit(sigmoid, V_dense, P_mean, p0=p0, bounds=bounds, maxfev=20000)
    return float(popt[0]), float(popt[1])


def transfer_function_F(cv_sweep, Delta: float, Vc0: float,
                        tw_ns: float, tau0_ns: float = 1.0,
                        V_dense=None, N_samples: int = 20000,
                        rng=None):
    """Sweep CV(Δ) values and return (beta_eff, vth_eff) arrays.

    F(CV) = beta_eff(CV) / beta_eff(CV=0) is the D2D broadening factor;
    callers can compute it from the returned beta_eff array.
    """
    cv_sweep = np.asarray(cv_sweep, dtype=float)
    if V_dense is None:
        V_dense = np.linspace(0.60, 1.10, 1500)
    if rng is None:
        rng = np.random.default_rng(seed=42)
    beta_eff = np.zeros_like(cv_sweep)
    vth_eff  = np.zeros_like(cv_sweep)
    for i, cv in enumerate(cv_sweep):
        P_mean = wafer_average_psw(V_dense, cv, Delta, Vc0,
                                   tw_ns=tw_ns, tau0_ns=tau0_ns,
                                   N_samples=N_samples, rng=rng)
        try:
            Vth, beta = fit_sigmoid_to_average(V_dense, P_mean)
        except RuntimeError:
            Vth, beta = np.nan, np.nan
        vth_eff[i]  = Vth
        beta_eff[i] = beta
    return beta_eff, vth_eff
