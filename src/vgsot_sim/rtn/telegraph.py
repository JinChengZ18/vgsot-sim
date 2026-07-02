"""Random-telegraph-noise (RTN) free-running superparamagnetic sMTJ model.

Where the LLG solver (`dynamic_switching`, `ser_cases`) answers *"given a write
pulse, did the cell switch?"* — a pulse-driven, one-shot question — this module
answers *"how does a low-barrier junction hop between its two states in
continuous time under a sustained bias?"*. That free-running, memory-bearing
dynamics is a candidate device-physics **node** primitive for reservoir computing,
which the pulse-switching cases do not provide. NB: this module models a SINGLE
node only — it is not a reservoir (no input weights, no node coupling, no readout
training; see ``scripts/10_rtn_reservoir/``).

Scope of ``V``
--------------
``V`` here is the *phenomenological effective bias / barrier-tilt* of the two-state
model, NOT literally the junction voltage ``V_MTJ``. In the package's own LLG path
``V_MTJ`` enters only through the *symmetric* VCMA term (it changes the barrier
*height*, not the inter-well *asymmetry*), and the real asymmetric drive is the SOT
*current* ``I_SOT`` (a Slonczewski/Sun quadratic law), so neither maps one-to-one
onto this linear antisymmetric tilt. The linear-in-``V`` rate law below is therefore
a *small-bias linearisation* valid for ``|V| << Vc0``; mapping a real device drive
(``I_SOT`` / VCMA) onto ``V`` requires a calibration that this package does not yet
provide (planned: ``rtn/bridge.py``, stage 2 of the research plan).

Physics
-------
A superparamagnetic junction sits in a double well and hops between states
``s in {-1, +1}`` as a continuous-time two-state Markov process. A bias tilts the
landscape, raising one Néel-Brown escape barrier and lowering the other:

    r_up(V) = (1/tau0) * exp[ -Delta * (1 - V/Vc0) ]    # -1 -> +1
    r_dn(V) = (1/tau0) * exp[ -Delta * (1 + V/Vc0) ]    # +1 -> -1

so ``r_up`` is the usual Néel-Brown rate and ``r_dn`` is the same law with the
bias reversed. The barrier argument is floored at 0 — exactly as
`vgsot_sim.analysis.nb_fit.psw_nb` does — so neither rate can exceed the attempt
frequency ``1/tau0`` (an escape rate faster than the attempt rate is unphysical).
The floor activates at ``|V| = Vc0`` where the barrier vanishes; for ``|V| > Vc0``
the well no longer exists, the two-state picture breaks down and the device is
deterministically pinned. **Treat ``|V| < Vc0`` as the physical domain.** Two closed
forms follow — the two knobs reservoir computing cares about:

* **Stationary mean** (input -> state nonlinearity)

      <s>_inf(V) = (r_up - r_dn)/(r_up + r_dn) = tanh(Delta * V / Vc0)

* **Correlation / relaxation time** (fading memory)

      tau(V) = 1 / (r_up + r_dn),

  maximal at V=0 (tau_max = tau0*exp(Delta)/2) and shrinking under strong drive,
  so bias trades memory (small |V|, long tau) against nonlinearity (large |V|,
  saturated tanh) — the canonical reservoir memory/nonlinearity tradeoff,
  grounded in device physics.

The per-step update uses the two-state propagator, which is EXACT when ``V`` is
held constant across the step (piecewise-constant input), so ``dt`` need not be
small relative to ``tau``:

    P(s_{t+dt} = +1 | s_t) = p_inf + (1[s_t=+1] - p_inf) * exp(-dt/tau),
    p_inf = r_up / (r_up + r_dn).

For a bias that varies *within* a step it becomes an adiabatic approximation,
accurate when ``dt`` is small vs. the input's variation timescale.

Units
-----
**Nanoseconds throughout**, matching :mod:`vgsot_sim.analysis.nb_fit` (whose
``tau0`` default is 1.0 ns and whose `fit_direction` returns ``Delta``/``Vc0``).
A :class:`NBFitResult` from `nb_fit.fit_direction` therefore drops straight into
:meth:`TelegraphParams.from_nb_fit`. Rates are in 1/ns and ``tau`` is in ns.

NumPy-only by design (the evolution is inherently sequential); no dependency on
the LLG kernels.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional, Union

import numpy as np

ArrayLike = Union[float, np.ndarray]


# ---------------------------------------------------------------------------
# Rate law (ns units) -- mirrors analysis.nb_fit's Neel-Brown convention
# ---------------------------------------------------------------------------
def neel_brown_rate(V: ArrayLike, *, tau0: float = 1.0,
                    Delta: ArrayLike = 5.15, Vc0: ArrayLike = 0.884) -> np.ndarray:
    """Néel-Brown escape rate W(V) = (1/tau0) exp[-max(Delta (1 - V/Vc0), 0)]  [1/ns].

    The barrier argument is floored at 0 (the same clip as
    :func:`vgsot_sim.analysis.nb_fit.psw_nb`) so the rate can never exceed the
    attempt frequency ``1/tau0`` — a thermal-activation escape faster than the
    attempt rate is unphysical. The floor activates at ``V = Vc0`` where the
    barrier vanishes; for ``V > Vc0`` the well no longer exists and the two-state
    picture breaks down (deterministic pinning). ``|V| < Vc0`` is the physical
    domain of this model.
    """
    V = np.asarray(V, dtype=np.float64)
    x = np.maximum(np.asarray(Delta) * (1.0 - V / np.asarray(Vc0)), 0.0)
    return (1.0 / tau0) * np.exp(-x)


def up_down_rates(V: ArrayLike, *, tau0: float = 1.0,
                  Delta: ArrayLike = 5.15, Vc0: ArrayLike = 0.884):
    """Escape rates ``(r_up, r_dn)`` [1/ns]; ``r_dn`` is ``r_up`` with V -> -V."""
    V = np.asarray(V, dtype=np.float64)
    r_up = neel_brown_rate(V, tau0=tau0, Delta=Delta, Vc0=Vc0)
    r_dn = neel_brown_rate(-V, tau0=tau0, Delta=Delta, Vc0=Vc0)
    return np.asarray(r_up, dtype=np.float64), np.asarray(r_dn, dtype=np.float64)


def stationary_mean(V: ArrayLike, *, Delta: ArrayLike = 5.15,
                    Vc0: ArrayLike = 0.884) -> np.ndarray:
    """Time-averaged state ``<s>_inf = tanh(Delta V / Vc0)`` in [-1, 1].

    Closed form valid in the physical domain ``|V| < Vc0``. Beyond ``Vc0`` one
    well has vanished and ``<s>`` is ``sign(V)`` deterministically (dynamics
    frozen, no hopping); the tanh form already returns ~±1 there. Note that beyond
    ``Vc0`` this closed form no longer equals ``(r_up - r_dn)/(r_up + r_dn)``
    computed from the *clipped* rates — the equality holds only on ``|V| < Vc0``.
    """
    V = np.asarray(V, dtype=np.float64)
    return np.tanh(np.asarray(Delta) * V / np.asarray(Vc0))


def relaxation_time(V: ArrayLike, *, tau0: float = 1.0,
                    Delta: ArrayLike = 5.15, Vc0: ArrayLike = 0.884) -> np.ndarray:
    """Correlation time ``tau(V) = 1/(r_up + r_dn)`` [ns] (the fading memory)."""
    r_up, r_dn = up_down_rates(V, tau0=tau0, Delta=Delta, Vc0=Vc0)
    return 1.0 / (r_up + r_dn)


def tau_max(*, tau0: float = 1.0, Delta: float = 5.15) -> float:
    """Zero-bias correlation time ``tau(0) = tau0 * exp(Delta) / 2`` [ns]."""
    return float(tau0 * np.exp(Delta) / 2.0)


# ---------------------------------------------------------------------------
@dataclass
class TelegraphParams:
    """Néel-Brown parameters (ns units) for the two-state RTN model.

    Defaults are the vgsot-sim §2.3.3 Device-A Néel-Brown fit (``Delta=5.15``,
    ``Vc0=0.884``); ``Vc0`` is the pulse-width-extracted SOT-equivalent critical
    voltage, and the model's ``V`` is the effective barrier-tilt it scales (not
    literally ``V_MTJ`` — see module docstring). For reservoir operation lower
    ``Delta`` (faster hopping, a few–tens of ns) — e.g. ``Delta≈3.8`` gives
    ``tau_max≈22 ns``; at the default ``Delta=5.15`` the tanh saturates by
    ``|V|≈0.3 V`` so the usable nonlinear/memory window is narrow.
    """
    tau0: float = 1.0     # attempt time [ns]
    Delta: float = 5.15   # thermal stability factor (dimensionless)
    Vc0: float = 0.884    # zero-thermal critical voltage [V]

    @classmethod
    def from_nb_fit(cls, fit, tau0: float = 1.0) -> "TelegraphParams":
        """Build from a :class:`vgsot_sim.analysis.nb_fit.NBFitResult`.

        ``fit.Delta`` / ``fit.Vc0`` come straight out of `fit_direction`; ``tau0``
        is the assumed attempt time in ns (same convention as nb_fit).
        """
        return cls(tau0=float(tau0), Delta=float(fit.Delta), Vc0=float(fit.Vc0))


class TelegraphArray:
    """A stateful population of ``n`` two-state superparamagnetic sMTJs.

    Each device is an independent two-state Markov node; one :meth:`step` advances
    the whole population by ``dt`` (ns) under a per-device bias and returns the new
    ``{-1, +1}`` states. This is the reservoir's pool of dynamical nodes.

    ``Delta`` / ``Vc0`` may be passed as per-device arrays of shape ``(n,)`` for a
    heterogeneous reservoir; otherwise the scalar ``params`` values are broadcast.
    """

    def __init__(self, n: int, params: Optional[TelegraphParams] = None, *,
                 Delta: Optional[np.ndarray] = None, Vc0: Optional[np.ndarray] = None,
                 seed: Optional[int] = None):
        if n <= 0:
            raise ValueError(f"n must be positive, got {n}")
        self.n = int(n)
        self.params = params or TelegraphParams()
        self.tau0 = float(self.params.tau0)
        self.Delta = (np.full(n, self.params.Delta, dtype=np.float64)
                      if Delta is None else np.asarray(Delta, dtype=np.float64))
        self.Vc0 = (np.full(n, self.params.Vc0, dtype=np.float64)
                    if Vc0 is None else np.asarray(Vc0, dtype=np.float64))
        if self.Delta.shape != (n,) or self.Vc0.shape != (n,):
            raise ValueError("Delta and Vc0 must have shape (n,) when supplied.")
        self.rng = np.random.default_rng(seed)
        self.reset()

    def reset(self, state: Optional[np.ndarray] = None) -> None:
        """Reset the population to ``state`` (default: random ``{-1, +1}``)."""
        if state is None:
            self.s = self.rng.choice(np.array([-1.0, 1.0]), size=self.n).astype(np.float64)
        else:
            self.s = np.asarray(state, dtype=np.float64).copy()

    def step(self, V: ArrayLike, dt: float) -> np.ndarray:
        """Advance all devices by ``dt`` [ns] under per-device bias ``V`` [V].

        The propagator is exact when ``V`` is held constant across the step
        (piecewise-constant input), so ``dt`` need not be small relative to
        ``tau``; for a bias that varies within the step it is an adiabatic
        approximation. ``V`` is the model's phenomenological effective bias/tilt
        (see module docstring), physical only for ``|V| < Vc0``. Returns a copy of
        the new ``{-1, +1}`` state vector ``(n,)``.
        """
        if dt <= 0:
            raise ValueError(f"dt must be positive, got {dt}")
        V = np.broadcast_to(np.asarray(V, dtype=np.float64), (self.n,))
        if np.any(np.abs(V) > self.Vc0):
            warnings.warn(
                "TelegraphArray.step: |V| > Vc0 is outside the two-state model's "
                "physical domain (barrier vanished, rate clipped to 1/tau0); the "
                "telegraph statistics there are not physical.",
                RuntimeWarning, stacklevel=2,
            )
        r_up, r_dn = up_down_rates(V, tau0=self.tau0, Delta=self.Delta, Vc0=self.Vc0)
        k = r_up + r_dn
        p_inf = r_up / k                       # stationary P(+1)
        decay = np.exp(-k * dt)
        is_plus = self.s > 0.0
        p_plus_next = np.where(is_plus, p_inf + (1.0 - p_inf) * decay,
                               p_inf * (1.0 - decay))
        u = self.rng.random(self.n)
        self.s = np.where(u < p_plus_next, 1.0, -1.0)
        return self.s.copy()

    @property
    def state(self) -> np.ndarray:
        """Current ``{-1, +1}`` state vector (no copy)."""
        return self.s


def simulate_trace(V_t: np.ndarray, dt: float, *, n: int = 1,
                   params: Optional[TelegraphParams] = None,
                   seed: Optional[int] = None) -> np.ndarray:
    """Run ``n`` RTN nodes through a bias time-series ``V_t`` [V] at step ``dt`` [ns].

    Returns the state trace, shape ``(len(V_t), n)`` in ``{-1, +1}``. The
    population-mean column ``out.mean(axis=1)`` is the reservoir read-out signal.
    """
    arr = TelegraphArray(n, params=params, seed=seed)
    V_t = np.asarray(V_t, dtype=np.float64)
    out = np.empty((V_t.shape[0], n), dtype=np.float64)
    for i, v in enumerate(V_t):
        out[i] = arr.step(v, dt)
    return out
