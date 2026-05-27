"""Process-variability budget + D2D Monte-Carlo for sMTJ arrays.

Combines the Brinkman-decomposed CV(Δ) budget (from PDK mismatch params)
with a Monte-Carlo sweep over CV(Δ) values, fitting a Sigmoid to each
wafer-averaged P_sw(V) curve and returning the broadening factor F.

Used by 07/variability_sim-V3.0.py; the analytic skeleton lives here so
the same machinery can be invoked from a vgsot-sim case or notebook.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import numpy as np
from scipy.optimize import curve_fit

from . import nb_fit
from .sigmoid_fit import sigmoid4p
from ..configs import PhysicalConstantsConfig


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


@dataclass
class MacrospinMismatchSample:
    """One sampled macrospin device plus the process scale factors used."""
    constants: PhysicalConstantsConfig
    scales: dict = field(default_factory=dict)


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
    total_var = CV_Delta ** 2

    def variance_fraction(component_var: float) -> float:
        return float(component_var / total_var) if total_var > 0.0 else 0.0

    return CVBudgetResult(
        CV_RA=CV_RA, CV_D=CV_D, CV_V=CV_V,
        CV_HK=CV_HK, CV_MS=inp.CV_MS, CV_TF=inp.CV_TF,
        CV_Delta=CV_Delta,
        components={
            # Split the volume term into lateral-area and thickness pieces.
            # CV_V already includes CV_TF, so reporting CV_V and CV_TF as
            # separate variance-budget bars would count thickness twice.
            "V_mag_area": variance_fraction((2 * CV_D) ** 2),
            "H_k":        variance_fraction(CV_HK ** 2),
            "M_s":        variance_fraction(inp.CV_MS ** 2),
            "t_f":        variance_fraction(inp.CV_TF ** 2),
        },
    )


def _positive_gaussian_scale(rng, cv: float, *, floor: float = 0.05,
                             z_score: float | None = None) -> float:
    """Draw a Gaussian multiplicative factor and keep physical values positive."""
    if cv <= 0.0:
        return 1.0
    z = rng.normal() if z_score is None else z_score
    return float(max(1.0 + cv * z, floor))


def sample_macrospin_process_constants(
    base: PhysicalConstantsConfig | None = None,
    inp: PDKBudgetInputs | None = None,
    *,
    rng=None,
    z_scores: dict | None = None,
) -> MacrospinMismatchSample:
    """Sample one macrospin device from the PDK process-mismatch budget.

    The NB-level variability figure reduces PDK mismatch to a CV(Delta)
    distribution.  Macrospin switching is sensitive to more knobs than
    Delta alone, so this helper maps the same budget onto solver parameters:

    - residual Rp mismatch after the Brinkman RA term -> D and D_elec
    - t_f and M_s process terms -> tf and Ms
    - the CV(TMR) anisotropy proxy used by the NB budget -> Ki
    - Rsot mismatch -> rho, hence R_W and voltage-to-current conversion
    - t_ox / phi terms -> tox, phi_bar and the first-order RA scale

    Actual TMR is scaled with the same CV(TMR) draw as the Ki proxy.  In a
    pure SOT, V_MTJ=0 sweep that transport scale is mostly diagnostic; Ki,
    D_elec, tf, Ms, and R_W carry the switching impact.  `z_scores` may
    override selected standard-normal draws (`D`, `tf`, `Ms`, `Ki_proxy`,
    `Rsot`, `tox`, `phi`) for stratified or antithetic populations.
    """
    base = base or PhysicalConstantsConfig()
    inp = inp or PDKBudgetInputs()
    if rng is None:
        rng = np.random.default_rng()

    budget = cv_delta_budget(inp)
    z_scores = z_scores or {}
    d_scale = _positive_gaussian_scale(rng, budget.CV_D, z_score=z_scores.get("D"))
    tf_scale = _positive_gaussian_scale(rng, inp.CV_TF, z_score=z_scores.get("tf"))
    ms_scale = _positive_gaussian_scale(rng, inp.CV_MS, z_score=z_scores.get("Ms"))
    hk_proxy_scale = _positive_gaussian_scale(
        rng, inp.CV_TMR, z_score=z_scores.get("Ki_proxy")
    )
    rsot_scale = _positive_gaussian_scale(rng, inp.CV_RSOT, z_score=z_scores.get("Rsot"))
    tox_scale = _positive_gaussian_scale(rng, inp.CV_TOX, z_score=z_scores.get("tox"))
    phi_scale = _positive_gaussian_scale(rng, inp.CV_PHI, z_score=z_scores.get("phi"))

    sens_tox = inp.KAPPA * np.sqrt(inp.PHI_BAR)
    sens_phi = inp.KAPPA * inp.T_OX / (2.0 * np.sqrt(inp.PHI_BAR))
    ln_ra = (
        sens_tox * inp.T_OX * (tox_scale - 1.0)
        + sens_phi * inp.PHI_BAR * (phi_scale - 1.0)
    )
    ra_scale = float(np.exp(ln_ra))

    constants = replace(
        base,
        D=base.D * d_scale,
        D_elec=base.D_elec * d_scale,
        tf=base.tf * tf_scale,
        Ms=base.Ms * ms_scale,
        Ki=base.Ki * hk_proxy_scale,
        H_k_eff_RT=base.H_k_eff_RT * hk_proxy_scale,
        rho=base.rho * rsot_scale,
        tox=base.tox * tox_scale,
        phi_bar=base.phi_bar * phi_scale,
        RA=base.RA * ra_scale,
        TMR=base.TMR * hk_proxy_scale,
    )
    return MacrospinMismatchSample(
        constants=constants,
        scales={
            "D": d_scale,
            "tf": tf_scale,
            "Ms": ms_scale,
            "Ki_proxy": hk_proxy_scale,
            "Rsot": rsot_scale,
            "tox": tox_scale,
            "phi": phi_scale,
            "RA": ra_scale,
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
