from .configs import PhysicalConstantsConfig
from .initialize import compute_Rp


def tmr_eff(V_MTJ, constants: PhysicalConstantsConfig):
    """Bias-dependent TMR(V_MTJ), selected by constants.tmr_model.

    Two forms are available:

      - "lorentzian": single-parameter Zhang-et-al Verilog-A form
            TMR(V) = TMR0 / (1 + (V/Vh)^2)
        Decays smoothly with V^2; useful when only V_h is known.

      - "pdk": three-parameter Hikstor SOT-MRAM PDK extraction
            TMR(V) = (TMR0 / k_TMR) * [1 / (a V^2 + b |V| + c) - 1]
        Matches measured non-symmetric decay around the write window.
        See doc Section 2.2.2.3 and parameter table 2.2.2.5.
    """
    model = constants.tmr_model
    if model == "lorentzian":
        return constants.TMR / (1.0 + (V_MTJ / constants.Vh) ** 2)
    if model == "pdk":
        denom = (constants.a_tmr * V_MTJ * V_MTJ
                 + constants.b_tmr * abs(V_MTJ)
                 + constants.c_tmr)
        return (constants.TMR / constants.k_tmr) * (1.0 / denom - 1.0)
    raise ValueError(
        f"Unknown tmr_model={model!r}. Use 'lorentzian' or 'pdk'."
    )


def tmr(V_MTJ, mz, constants: PhysicalConstantsConfig,
        include_series: bool = True):
    """MTJ resistance R(m_z, V_MTJ) via conductance interpolation.

        G(m_z) = (1+m_z)/2 * G_P + (1-m_z)/2 * G_AP(V)
        R(m_z, V) = R_P * (1 + TMR_eff(V)) / (1 + TMR_eff(V) * (1+m_z)/2)

    Code convention (consistent with the pre-existing Lorentzian form):
        m_z = +1 -> parallel (R = R_P, bias-independent)
        m_z = -1 -> antiparallel (R = R_P * (1 + TMR_eff(V)) = R_AP(V))

    TMR(V) form is chosen via constants.tmr_model. Equivalent to the prior
    hard-coded Lorentzian when tmr_model="lorentzian".

    Series-resistance toggle
    ------------------------
    When `include_series=True` (default) and `constants.R_series` is
    nonzero, the contact / lead parasitic R_series is added on top of the
    intrinsic tunnel resistance. Set `include_series=False` to recover
    the bare R_MTJ for model/test purposes. Setting `R_series=0` in the
    config gives byte-identical behaviour to the pre-toggle code, so
    existing simulations are unaffected.
    """
    R_p = compute_Rp(constants)
    Teff = tmr_eff(V_MTJ, constants)
    R_intrinsic = R_p * (1.0 + Teff) / (1.0 + Teff * 0.5 * (1.0 + mz))
    R_extra = float(getattr(constants, "R_series", 0.0)) if include_series else 0.0
    return R_intrinsic + R_extra
