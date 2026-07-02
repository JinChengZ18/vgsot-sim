"""sLLG ↔ RTN bridge — does the macrospin engine reproduce the 2-state telegraph?

Stage 2 of the RC roadmap (``scripts/10_rtn_reservoir/README.md``). The RTN model
(:mod:`vgsot_sim.rtn.telegraph`) is a *reduced* two-state abstraction of a
low-barrier sMTJ. This module drives the repo's OWN stochastic LLG engine
(:func:`vgsot_sim.dynamic_switching_vector.switching_vector`) in a free-running,
low-barrier regime and checks whether the emergent statistics match the RTN
closed forms:

  (a) ``<m_z>(bias) → tanh``     — the input→state nonlinearity,
  (b) dwell-time distribution → exponential (Poisson escape),
  (c) power spectral density → Lorentzian (the two-state RTN signature),

and lets us calibrate the device attempt time ``tau0`` and the ``bias → V_rtn``
map that the RTN model leaves phenomenological.

Getting the barrier right
-------------------------
The perpendicular thermal-stability barrier is

    Delta = K_u^eff · v / (kB·T),   K_u^eff = Ki/tf − 0.5·mu0·Ms²·(N_z − N_x)

(see :func:`vgsot_sim.material_temperature.k_u_eff_of_T`). The shape-anisotropy
term is **not** negligible here: for the default pillar ``N_z − N_x ≈ 0.96``,
worth ~208 ``kB·T``. Setting ``Ki`` tiny to fake a "PMA-only" low ``Delta`` drives
``K_u^eff`` negative → the easy axis goes **in-plane** and there is no
perpendicular telegraph at all. :func:`ki_for_delta` inverts the *correct*
relation, keeping the well perpendicular.

The bias knob
-------------
The symmetry-breaking bias used here is a longitudinal field ``h_ex_z`` (a Zeeman
tilt of the two wells) — the transparent analog of the RTN's antisymmetric tilt.
In the two-state limit the occupancy gives

    <m_z> = tanh( mu0·Ms·h_ex_z·v / (kB·T) ),

so the RTN's ``V`` maps *linearly* to ``h_ex_z`` with slope :func:`tilt_per_field`.
This is the honest minimal mapping: ``V_MTJ`` enters the LLG only through the
*symmetric* VCMA term (it cannot tilt the wells) and the SOT current drives via a
damping torque rather than a simple energy tilt (see the audit and
``scripts/10_rtn_reservoir/README.md``). Mapping the RTN ``V`` onto a real terminal drive is
therefore a separate calibration, not a one-to-one identity.

NumPy-only; the evolution is inherently sequential.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional, Sequence

import numpy as np

from ..configs import PhysicalConstantsConfig
from ..demag import demag_factors
from ..dynamic_switching_vector import switching_vector
from ..material_temperature import k_u_eff_of_T


# ---------------------------------------------------------------------------
# Barrier ↔ Ki and bias ↔ tilt
# ---------------------------------------------------------------------------
def delta_of_ki(Ki: float, constants: PhysicalConstantsConfig) -> float:
    """Perpendicular thermal-stability factor ``Delta = K_u^eff·v/(kB·T)`` for ``Ki``."""
    c = replace(constants, Ki=float(Ki))
    return k_u_eff_of_T(c.T, c) * c.v / (c.kb * c.T)


def ki_for_delta(target_delta: float, constants: PhysicalConstantsConfig) -> float:
    """Interface anisotropy ``Ki`` [J/m²] giving a perpendicular barrier ``target_delta``.

    Inverts ``Delta = (Ki/tf − 0.5·mu0·Ms²·(N_z−N_x))·v/(kB·T)`` including the demag
    term, so the resulting well stays perpendicular (unlike a PMA-only estimate).
    """
    Nx, _, Nz = demag_factors(constants, mode="ellipsoid")
    K_eff = target_delta * constants.kb * constants.T / constants.v
    return (K_eff + 0.5 * constants.u0 * constants.Ms ** 2 * (Nz - Nx)) * constants.tf


def tilt_per_field(constants: PhysicalConstantsConfig) -> float:
    """Dimensionless Zeeman tilt per unit longitudinal field: ``mu0·Ms·v/(kB·T)`` [m/A]."""
    return constants.u0 * constants.Ms * constants.v / (constants.kb * constants.T)


def low_barrier_constants(target_delta: float,
                          constants: Optional[PhysicalConstantsConfig] = None,
                          *, h_ex_z: float = 0.0) -> PhysicalConstantsConfig:
    """A copy of ``constants`` retuned to perpendicular barrier ``target_delta``.

    In-plane bias fields are zeroed (free evolution along the easy axis) and the
    longitudinal tilt field ``h_ex_z`` is set; everything else is inherited.
    """
    c = constants or PhysicalConstantsConfig()
    return replace(c, Ki=ki_for_delta(target_delta, c),
                   h_ex_x=0.0, h_ex_y=0.0, h_ex_z=float(h_ex_z))


# ---------------------------------------------------------------------------
# Free-running LLG driver
# ---------------------------------------------------------------------------
def free_run(constants: PhysicalConstantsConfig, n_steps: int, *,
             h_ex_z: float = 0.0, seed: Optional[int] = None,
             m0: Optional[np.ndarray] = None) -> np.ndarray:
    """Free-running (``I_SOT=0``, ``V_MTJ=0``, thermal noise on) m_z trace.

    Advances :func:`switching_vector` for ``n_steps`` under longitudinal bias
    ``h_ex_z`` and returns the ``m_z(t)`` array (length ``n_steps``), sampled every
    ``constants.t_step`` seconds. VCMA is off; the only drive is thermal noise plus
    the optional Zeeman tilt.
    """
    if n_steps <= 0:
        raise ValueError(f"n_steps must be positive, got {n_steps}")
    cc = replace(constants, h_ex_z=float(h_ex_z))
    rng = np.random.default_rng(seed)
    if m0 is None:
        m = np.array([0.02, 0.0, 0.9998])
    else:
        m = np.asarray(m0, dtype=float)
    m = m / np.linalg.norm(m)
    mz = np.empty(int(n_steps), dtype=float)
    for i in range(int(n_steps)):
        m = switching_vector(m, 0.0, 0.0, 5000.0, 0, 0, VNV=0, NON=1,
                             constants=cc, rng=rng)
        mz[i] = m[2]
    return mz


# ---------------------------------------------------------------------------
# Statistics: dwell times, <m_z>, PSD
# ---------------------------------------------------------------------------
def _hysteretic_state(mz: np.ndarray, thr: float) -> np.ndarray:
    """Forward-filled ±1 state from m_z with a hysteretic dead-band ``|m_z|<thr``."""
    state = np.where(mz > thr, 1, np.where(mz < -thr, -1, 0)).astype(np.int8)
    last = 0
    filled = np.empty_like(state)
    for i, s in enumerate(state):
        if s != 0:
            last = s
        filled[i] = last
    return filled


@dataclass
class DwellStats:
    n_flips: int
    dwell_ns: np.ndarray   # per-dwell durations
    tau_mean_ns: float     # mean dwell (≈ retention time tau_ret)
    cv: float              # coefficient of variation (→1 for exponential/Poisson)


def dwell_times(mz: np.ndarray, t_step: float, *, thr: float = 0.5) -> DwellStats:
    """Extract dwell-time statistics via a hysteretic ``|m_z|>thr`` state.

    Using a hysteresis threshold (not ``sign(m_z)``) avoids counting equator
    jitter as escapes — the failure mode flagged in the audit. ``t_step`` in
    seconds; dwell durations are returned in ns.
    """
    dt_ns = t_step * 1e9
    filled = _hysteretic_state(np.asarray(mz, float), thr)
    change = np.where(np.diff(filled) != 0)[0]
    if len(change) < 2:
        return DwellStats(int(len(change)), np.array([]), float("nan"), float("nan"))
    dwell = np.diff(change) * dt_ns
    return DwellStats(int(len(change)), dwell, float(dwell.mean()),
                      float(dwell.std() / dwell.mean()))


def infer_tau0(tau_mean_ns: float, delta: float) -> float:
    """Device attempt time from the measured retention: ``tau0 = tau_ret / exp(Delta)`` [ns]."""
    return float(tau_mean_ns / np.exp(delta))


@dataclass
class PSDFit:
    f: np.ndarray
    P: np.ndarray
    f_corner_hz: float
    tau_corner_ns: float


def psd_lorentzian(mz: np.ndarray, t_step: float, *, n_seg: int = 8,
                   f_max_hz: float = 5e7) -> PSDFit:
    """Welch-averaged PSD of ``m_z`` and a best-effort Lorentzian corner.

    A two-state RTN has a Lorentzian spectrum ``P(f) = P0/(1+(f/f_c)²)`` with
    ``f_c = 1/(2π·tau_c)``, ``tau_c`` the autocorrelation time. The corner is
    estimated as the half-power crossing of the smoothed, band-limited PSD
    (``f < f_max_hz``, to stay below the tens-of-GHz intra-well precession band).

    Caveat: the macrospin spectrum hits a white-noise floor a few × ``f_c`` above
    the corner (per-step thermal jitter), so the direct PSD corner is only
    accurate to ~2× and depends on ``f_max_hz``. For a quantitative corner use the
    dwell-time route — for a symmetric telegraph ``tau_c = tau_dwell/2`` from
    :func:`dwell_times` — which is robust. The PSD here is mainly for confirming the
    Lorentzian *shape*.
    """
    x = np.asarray(mz, float) - np.mean(mz)
    n = len(x)
    seg = max(1, n // n_seg)
    acc = None
    k = 0
    for s in range(0, n - seg + 1, seg):
        X = np.fft.rfft(x[s:s + seg])
        P = (np.abs(X) ** 2) / seg
        acc = P if acc is None else acc + P
        k += 1
    P = acc / max(k, 1)
    f = np.fft.rfftfreq(seg, d=t_step)
    f_c = float("nan")
    band = (f > 0) & (f <= f_max_hz)
    if np.count_nonzero(band) >= 8:
        fb, Pb = f[band], P[band]
        w = max(3, len(Pb) // 50)
        Ps = np.convolve(Pb, np.ones(w) / w, mode="same")
        plateau = float(np.median(Ps[:max(2, len(Ps) // 20)]))
        below = np.where(Ps < plateau / 2.0)[0]
        if len(below):
            f_c = float(fb[below[0]])
    tau_c = (1.0 / (2.0 * np.pi * f_c) * 1e9) if f_c and np.isfinite(f_c) else float("nan")
    return PSDFit(f, P, f_c, tau_c)


def mean_mz_vs_bias(constants: PhysicalConstantsConfig, hz_list: Sequence[float],
                    n_steps: int, *, seeds: Sequence[int] = (0, 1, 2),
                    m0: Optional[np.ndarray] = None):
    """Ensemble-averaged ``<m_z>`` vs longitudinal bias, with the tanh prediction.

    For each ``h_ex_z`` runs one trajectory per seed (seeds split between the two
    wells so the ensemble is symmetric at zero bias) and averages the time-mean
    ``m_z``. Returns ``(hz, tilt, mz_mean, tanh_pred)`` arrays, where
    ``tilt = tilt_per_field·h_ex_z`` and ``tanh_pred = tanh(tilt)``.

    NB: faithful equilibrium sampling needs each run ≫ the retention time; at high
    ``Delta`` (long dwell) use a large ``n_steps`` or this will be init-biased.
    """
    slope = tilt_per_field(constants)
    hz_arr, tilt_arr, mz_arr, tanh_arr = [], [], [], []
    for hz in hz_list:
        vals = []
        for j, sd in enumerate(seeds):
            start = np.array([0.02, 0.0, 0.9998]) if j % 2 == 0 else np.array([0.02, 0.0, -0.9998])
            vals.append(free_run(constants, n_steps, h_ex_z=hz, seed=sd,
                                  m0=(m0 if m0 is not None else start)).mean())
        tilt = slope * hz
        hz_arr.append(float(hz)); tilt_arr.append(float(tilt))
        mz_arr.append(float(np.mean(vals))); tanh_arr.append(float(np.tanh(tilt)))
    return (np.array(hz_arr), np.array(tilt_arr), np.array(mz_arr), np.array(tanh_arr))
