"""RTN physical reservoir — input projection + heterogeneous fading-memory nodes + ridge readout.

Stage 3 of the RC roadmap (``scripts/10_rtn_reservoir/README.md``). This is the layer the audit found missing
(A1–A6): the telegraph module alone is a single node, and broadcasting one bias to
identical independent nodes collapses to a single low-pass tanh (A2). A reservoir
needs every node to see a DISTINCT projection of the input and to span a range of
fading-memory timescales.

Construction
------------
- **Input projection.** Node ``j`` sees bias ``V_j(t) = a_in·w_j·u(t) + b_j`` with
  random input weights ``w_j`` and per-node offsets ``b_j`` (the W_in masking that
  A2 said was missing). Inputs are scaled so ``|V_j| < Vc0`` (the physical domain).
- **Heterogeneous nodes.** Per-node ``(Delta_j, Vc0_j)`` give distinct nonlinear
  gains ``Delta_j/Vc0_j`` and — crucially — a spread of fading-memory time constants
  ``tau_j(V) = 1/(r↑+r↓)`` (A3: heterogeneity alone is not enough without W_in, but
  together they span a high-dimensional state space).
- **Readout.** A trained linear ridge map from the node-state design matrix.

State modes
-----------
- ``"meanfield"`` (default): deterministic expected-state leaky integrator
  ``x_j(t) = s∞_j(t) + (x_j(t-1) − s∞_j(t))·exp(−dt/tau_j(t))`` with
  ``s∞_j = tanh(Delta_j V_j/Vc0_j)``. This is the mean field of the RTN node
  (= averaging many physical devices / long time), clean and reproducible — best for
  *measuring* capacity.
- ``"stochastic"``: actual binary :class:`TelegraphArray` states, optionally
  ``n_replicas`` physical devices per logical node, averaged — exposes the device
  sampling-noise penalty on capacity.

All timescales are in ns (the telegraph convention); ``dt`` is the input step in ns.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .telegraph import TelegraphArray, TelegraphParams, relaxation_time, stationary_mean


@dataclass
class ReservoirConfig:
    n_nodes: int = 100
    tau0: float = 1.0          # attempt time [ns]
    delta_range: tuple = (0.5, 4.0)   # per-node Delta ~ U(lo, hi) -> spread of memory timescales
    vc0: float = 0.884
    a_in: float = 0.5          # input scaling (keeps |V| < Vc0 for u in [-1,1])
    bias_spread: float = 0.25  # per-node offset b_j ~ U(-spread, spread) [V]
    w_binary: bool = True      # input weights in {-1,+1} (else N(0,1))
    dt: float = 1.0            # input time step [ns] (best linear-MC operating point)
    # ── optional inter-node coupling (mean-field mode only) ─────────────
    # Each node's bias gains a recurrent term (W_res @ x)_j, ESN-style
    # (0 = no coupling, the independent filter bank). Hardware reading: node
    # j's drive includes weighted read-outs of other nodes' states — a
    # routing/drive overhead a physical design must pay.
    # Topologies (measured, see scripts/10_rtn_reservoir/README.md):
    #   "random": sparse Gaussian rescaled to spectral radius `coupling_radius`.
    #             HURTS linear MC at every (radius, scale) tested.
    #   "ring":   simple cycle i -> i+1 with weight `coupling_radius`; combine
    #             with input_mode="single" and homogeneous low Delta. BEATS the
    #             filter-bank MC ceiling ~4x (MC~32 at n=100) and MC grows
    #             with n. Best when the per-hop small-signal gain
    #             g = coupling_radius*coupling_scale_v*Delta/Vc0 sits just
    #             below 1 (lossless-ish propagation without saturation).
    coupling_radius: float = 0.0
    coupling_topology: str = "random"   # "random" | "ring"
    coupling_density: float = 0.1       # (random topology only)
    coupling_scale_v: float = 0.5   # volts of bias per unit recurrent activation
    input_mode: str = "all"             # "all" | "single" (input enters node 0 only)


class Reservoir:
    """A heterogeneous bank of RTN fading-memory nodes with random input projection."""

    def __init__(self, cfg: Optional[ReservoirConfig] = None, *, seed: Optional[int] = None):
        self.cfg = cfg or ReservoirConfig()
        n = self.cfg.n_nodes
        rng = np.random.default_rng(seed)
        self.rng = rng
        lo, hi = self.cfg.delta_range
        self.Delta = rng.uniform(lo, hi, size=n)
        self.Vc0 = np.full(n, self.cfg.vc0, dtype=float)
        self.tau0 = float(self.cfg.tau0)
        self.dt = float(self.cfg.dt)
        if self.cfg.w_binary:
            self.W_in = rng.choice(np.array([-1.0, 1.0]), size=n)
        else:
            self.W_in = rng.standard_normal(n)
        if self.cfg.input_mode == "single":       # input enters node 0 only
            self.W_in = np.zeros(n)
            self.W_in[0] = 1.0
        elif self.cfg.input_mode != "all":
            raise ValueError(f"unknown input_mode={self.cfg.input_mode!r}")
        self.bias = rng.uniform(-self.cfg.bias_spread, self.cfg.bias_spread, size=n)
        self.a_in = float(self.cfg.a_in)
        self.n = n
        # optional recurrent coupling matrix
        self.W_res = None
        if self.cfg.coupling_radius > 0.0:
            if self.cfg.coupling_topology == "ring":
                W = np.zeros((n, n))
                for i in range(n):
                    W[(i + 1) % n, i] = self.cfg.coupling_radius
                self.W_res = W
            elif self.cfg.coupling_topology == "random":
                W = rng.standard_normal((n, n))
                W *= (rng.random((n, n)) < self.cfg.coupling_density)
                np.fill_diagonal(W, 0.0)
                rho = float(np.max(np.abs(np.linalg.eigvals(W))))
                if rho > 0:
                    self.W_res = W * (self.cfg.coupling_radius / rho)
            else:
                raise ValueError(
                    f"unknown coupling_topology={self.cfg.coupling_topology!r}")

    def _bias_of(self, u_t: float) -> np.ndarray:
        """Per-node bias V_j for scalar input u_t, clipped to the physical domain."""
        V = self.a_in * self.W_in * u_t + self.bias
        lim = 0.98 * self.Vc0
        return np.clip(V, -lim, lim)

    def run(self, u, *, mode: str = "meanfield", washout: int = 0,
            n_replicas: int = 1, seed: Optional[int] = None) -> np.ndarray:
        """Drive the reservoir with input series ``u`` (T,); return states (T−washout, n).

        ``mode="meanfield"`` integrates the deterministic expected state;
        ``mode="stochastic"`` samples binary device states (``n_replicas`` averaged).
        """
        u = np.asarray(u, dtype=float)
        T = u.shape[0]
        X = np.empty((T, self.n), dtype=float)
        if mode == "meanfield":
            x = np.zeros(self.n)
            lim = 0.98 * self.Vc0
            for t in range(T):
                V = self._bias_of(u[t])
                if self.W_res is not None:
                    V = np.clip(V + self.cfg.coupling_scale_v * (self.W_res @ x),
                                -lim, lim)
                s_inf = stationary_mean(V, Delta=self.Delta, Vc0=self.Vc0)
                tau = relaxation_time(V, tau0=self.tau0, Delta=self.Delta, Vc0=self.Vc0)
                decay = np.exp(-self.dt / tau)
                x = s_inf + (x - s_inf) * decay
                X[t] = x
        elif mode == "stochastic":
            if self.W_res is not None:
                raise NotImplementedError(
                    "inter-node coupling is implemented for the mean-field mode only")
            R = max(1, int(n_replicas))
            Delta_t = np.tile(self.Delta, R)
            Vc0_t = np.tile(self.Vc0, R)
            arr = TelegraphArray(self.n * R, TelegraphParams(tau0=self.tau0, Vc0=self.cfg.vc0),
                                 Delta=Delta_t, Vc0=Vc0_t, seed=seed)
            for t in range(T):
                V = np.tile(self._bias_of(u[t]), R)
                s = arr.step(V, self.dt)
                X[t] = s.reshape(R, self.n).mean(axis=0)
        else:
            raise ValueError(f"unknown mode={mode!r}")
        return X[washout:]


# ---------------------------------------------------------------------------
# Ridge readout
# ---------------------------------------------------------------------------
def ring_reservoir(n_nodes: int = 100, *, delta: float = 1.0, hop_gain: float = 0.87,
                   w: float = 0.9, seed: Optional[int] = None) -> Reservoir:
    """Tuned simple-cycle (delay-line) reservoir — the measured way past the
    filter-bank MC ceiling (see scripts/10_rtn_reservoir/README.md).

    Homogeneous low-``delta`` nodes on a unidirectional ring, input injected into
    node 0 only. ``hop_gain`` is the ring-weight x bias-scale product; the per-hop
    small-signal gain ``hop_gain*delta/Vc0`` should sit just below 1 (stable,
    near-lossless propagation). Defaults give MC ~ 32 at n=100 (seed=3).
    """
    cfg = ReservoirConfig(n_nodes=n_nodes, delta_range=(delta, delta),
                          bias_spread=0.0, coupling_radius=w,
                          coupling_topology="ring", coupling_scale_v=hop_gain / w,
                          input_mode="single")
    return Reservoir(cfg, seed=seed)


def ridge_fit(X: np.ndarray, Y: np.ndarray, alpha: float = 1e-6) -> np.ndarray:
    """Closed-form ridge weights ``W`` for ``Y ≈ [X, 1] @ W`` (bias column appended)."""
    Xb = np.hstack([X, np.ones((X.shape[0], 1))])
    A = Xb.T @ Xb + alpha * np.eye(Xb.shape[1])
    return np.linalg.solve(A, Xb.T @ Y)


def ridge_predict(X: np.ndarray, W: np.ndarray) -> np.ndarray:
    Xb = np.hstack([X, np.ones((X.shape[0], 1))])
    return Xb @ W


def _r2_corr(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Squared Pearson correlation (the memory-capacity / IPC metric)."""
    yt = y_true - y_true.mean()
    yp = y_pred - y_pred.mean()
    denom = np.sqrt((yt @ yt) * (yp @ yp))
    if denom == 0:
        return 0.0
    return float((yt @ yp) ** 2 / (denom ** 2))


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------
@dataclass
class MemoryCapacityResult:
    mc_total: float
    mc_k: np.ndarray          # per-delay MC_k
    delays: np.ndarray


def memory_capacity(reservoir: Reservoir, *, n_samples: int = 2500, max_delay: int = 40,
                    washout: int = 200, train_frac: float = 0.6, alpha: float = 1e-6,
                    mode: str = "meanfield", n_replicas: int = 1,
                    seed: int = 0) -> MemoryCapacityResult:
    """Linear short-term memory capacity (Jaeger): MC = Σ_k corr²(û[t−k], u[t−k]).

    i.i.d. uniform input; for each delay ``k`` a separate ridge readout reconstructs
    ``u[t−k]`` and ``MC_k`` is the test-set squared correlation. ``MC`` is bounded by
    the number of linearly independent reservoir signals (≤ n_nodes).
    """
    rng = np.random.default_rng(seed)
    u = rng.uniform(-1.0, 1.0, size=n_samples + washout)
    X = reservoir.run(u, mode=mode, washout=washout, n_replicas=n_replicas, seed=seed)
    u = u[washout:]
    T = X.shape[0]
    n_train = int(train_frac * T)
    mc_k = np.zeros(max_delay)
    for k in range(1, max_delay + 1):
        Xk, yk = X[k:], u[:-k]                       # state at t vs input at t−k
        Xtr, ytr = Xk[:n_train], yk[:n_train]
        Xte, yte = Xk[n_train:], yk[n_train:]
        W = ridge_fit(Xtr, ytr[:, None], alpha=alpha)
        yp = ridge_predict(Xte, W)[:, 0]
        mc_k[k - 1] = _r2_corr(yte, yp)
    return MemoryCapacityResult(float(mc_k.sum()), mc_k, np.arange(1, max_delay + 1))


def narma10(u: np.ndarray) -> np.ndarray:
    """NARMA-10 target series for input ``u`` ∈ [0, 0.5]. Standard 10th-order benchmark."""
    u = np.asarray(u, float)
    y = np.zeros_like(u)
    for t in range(len(u)):
        if t < 10:
            y[t] = 0.0
            continue
        y[t] = (0.3 * y[t - 1] + 0.05 * y[t - 1] * np.sum(y[t - 10:t])
                + 1.5 * u[t - 10] * u[t - 1] + 0.1)
    return y


@dataclass
class TaskResult:
    nrmse: float
    r2: float


def narma10_task(reservoir: Reservoir, *, n_samples: int = 3000, washout: int = 200,
                 train_frac: float = 0.6, alpha: float = 1e-6, mode: str = "meanfield",
                 n_replicas: int = 1, seed: int = 0) -> TaskResult:
    """Train+test the reservoir on NARMA-10; return NRMSE and R²."""
    rng = np.random.default_rng(seed)
    u = rng.uniform(0.0, 0.5, size=n_samples + washout)
    y = narma10(u)
    X = reservoir.run(u, mode=mode, washout=washout, n_replicas=n_replicas, seed=seed)
    y = y[washout:]
    T = X.shape[0]
    n_train = int(train_frac * T)
    W = ridge_fit(X[:n_train], y[:n_train, None], alpha=alpha)
    yp = ridge_predict(X[n_train:], W)[:, 0]
    yte = y[n_train:]
    nrmse = float(np.sqrt(np.mean((yp - yte) ** 2)) / (np.std(yte) + 1e-12))
    return TaskResult(nrmse, _r2_corr(yte, yp))


def parity_task(reservoir, *, n_bits: int = 3, n_samples: int = 3000, washout: int = 200,
                train_frac: float = 0.6, alpha: float = 1e-6, mode: str = "meanfield",
                n_replicas: int = 1, seed: int = 0) -> float:
    """Delayed N-bit parity (a hard nonlinear+memory task); returns test accuracy.

    Input is a random ±1 bit stream; the target at ``t`` is the product (parity) of
    the last ``n_bits`` inputs. A linear readout cannot solve parity without the
    reservoir's nonlinear mixing, so accuracy ≫ 0.5 demonstrates nonlinear memory.
    """
    rng = np.random.default_rng(seed)
    bits = rng.choice(np.array([-1.0, 1.0]), size=n_samples + washout)
    X = reservoir.run(bits, mode=mode, washout=washout, n_replicas=n_replicas, seed=seed)
    bits = bits[washout:]
    y = np.ones(len(bits))
    for k in range(n_bits):
        y[k:] *= bits[:len(bits) - k]               # product of last n_bits inputs
    X, y = X[n_bits:], y[n_bits:]
    T = X.shape[0]
    n_train = int(train_frac * T)
    W = ridge_fit(X[:n_train], y[:n_train, None], alpha=alpha)
    yp = np.sign(ridge_predict(X[n_train:], W)[:, 0])
    return float(np.mean(yp == y[n_train:]))


def _legendre_norm(u: np.ndarray, degree: int) -> np.ndarray:
    """Orthonormal Legendre polynomial on Uniform[-1,1]: degree 1 → √3 u, 2 → √5(3u²−1)/2."""
    if degree == 1:
        return np.sqrt(3.0) * u
    if degree == 2:
        return np.sqrt(5.0) * (3.0 * u ** 2 - 1.0) / 2.0
    raise ValueError("degree must be 1 or 2")


@dataclass
class IPCResult:
    total: float
    deg1: float            # = linear memory capacity
    deg2: float            # quadratic (squares + pairwise products)


def information_processing_capacity(reservoir: Reservoir, *, max_delay: int = 15,
                                    pair_delay: int = 8, n_samples: int = 4000,
                                    washout: int = 200, train_frac: float = 0.6,
                                    alpha: float = 1e-6, threshold: float = 0.02,
                                    mode: str = "meanfield", seed: int = 0) -> IPCResult:
    """Approximate Dambre-2012 IPC over degree-1 and degree-2 Legendre targets.

    Sums test-set R² over an orthonormal basis of past-input functions: degree-1
    (= linear MC), degree-2 squares ``ℓ2(u[t−k])``, and pairwise products
    ``ℓ1(u[t−i])·ℓ1(u[t−j])``. A ``threshold`` cutoff (drop R² below it) suppresses
    the finite-sample noise floor — a practical surrogate for Dambre's shuffle test,
    so this is an approximate IPC, not the exact value.
    """
    rng = np.random.default_rng(seed)
    u = rng.uniform(-1.0, 1.0, size=n_samples + washout)
    X = reservoir.run(u, mode=mode, washout=washout, seed=seed)
    u = u[washout:]
    T = X.shape[0]
    n_train = int(train_frac * T)

    def cap(target: np.ndarray, shift: int) -> float:
        Xs, ys = X[shift:], target[:len(target) - shift] if shift else target
        W = ridge_fit(Xs[:n_train], ys[:n_train, None], alpha=alpha)
        yp = ridge_predict(Xs[n_train:], W)[:, 0]
        r2 = _r2_corr(ys[n_train:], yp)
        return r2 if r2 >= threshold else 0.0

    l1 = _legendre_norm(u, 1)
    l2 = _legendre_norm(u, 2)
    deg1 = sum(cap(np.roll(l1, k), 0) for k in range(1, max_delay + 1))
    deg2 = sum(cap(np.roll(l2, k), 0) for k in range(1, max_delay + 1))
    for i in range(1, pair_delay + 1):                # pairwise products ℓ1[t-i]·ℓ1[t-j]
        for j in range(i + 1, pair_delay + 1):
            deg2 += cap(np.roll(l1, i) * np.roll(l1, j), 0)
    return IPCResult(float(deg1 + deg2), float(deg1), float(deg2))


def kernel_quality(reservoir, *, n_streams: int = 40, t_len: int = 60,
                   washout: int = 150, sv_threshold: float = 1e-4,
                   generalization: bool = False, seed: int = 0) -> int:
    """Kernel rank (separation) or generalization rank of the reservoir.

    Drives ``n_streams`` input sequences and returns the effective rank
    (singular values above ``sv_threshold`` x the largest) of the matrix of final
    states. ``generalization=False``: fully independent random streams — higher
    rank = richer separation (bounded by min(n_streams, n_nodes)).
    ``generalization=True``: streams share the last ``t_len//4`` inputs and differ
    only in the distant past — LOWER rank = better generalization (the reservoir
    should forget the irrelevant far past). Legenstein & Maass's kernel /
    generalization pair.
    """
    rng = np.random.default_rng(seed)
    common_tail = rng.uniform(-1.0, 1.0, size=t_len // 4)
    finals = np.empty((n_streams, reservoir.n))
    for k in range(n_streams):
        u = rng.uniform(-1.0, 1.0, size=washout + t_len)
        if generalization:
            u[-len(common_tail):] = common_tail
        X = reservoir.run(u, washout=washout)
        finals[k] = X[-1]
    sv = np.linalg.svd(finals - finals.mean(axis=0), compute_uv=False)
    if sv[0] <= 0:
        return 0
    return int(np.sum(sv > sv_threshold * sv[0]))


def esp_convergence(reservoir, *, t_len: int = 400, seed: int = 0) -> float:
    """Echo-state-property check: distance between two state trajectories started
    from different initial conditions under the SAME input, at the final step.
    Values ~0 mean the initial condition is forgotten (ESP holds).

    Mean-field mode only (deterministic given the input)."""
    rng = np.random.default_rng(seed)
    u = rng.uniform(-1.0, 1.0, size=t_len)
    lim = 0.98 * reservoir.Vc0
    xs = []
    for x0 in (np.full(reservoir.n, -0.9), np.full(reservoir.n, +0.9)):
        x = x0.copy()
        for t in range(t_len):
            V = reservoir._bias_of(u[t])
            if reservoir.W_res is not None:
                V = np.clip(V + reservoir.cfg.coupling_scale_v * (reservoir.W_res @ x),
                            -lim, lim)
            s_inf = stationary_mean(V, Delta=reservoir.Delta, Vc0=reservoir.Vc0)
            tau = relaxation_time(V, tau0=reservoir.tau0, Delta=reservoir.Delta,
                                  Vc0=reservoir.Vc0)
            x = s_inf + (x - s_inf) * np.exp(-reservoir.dt / tau)
        xs.append(x)
    return float(np.linalg.norm(xs[0] - xs[1]) / np.sqrt(reservoir.n))


class DelayReservoir:
    """Single-node time-delay reservoir (Appeltant 2011) built on one RTN node.

    One physical nonlinear node is time-multiplexed into ``n_virtual`` virtual nodes
    by a random input mask, with delayed feedback that couples the virtual chain —
    the recurrence that the independent-node :class:`Reservoir` (a filter bank)
    lacks.

    Measured trade-off (see ``scripts/10_rtn_reservoir/README.md``): with feedback +
    per-node bias this recurrence yields *higher nonlinear* capacity than the filter
    bank (IPC deg-2 ≈ 3.5 vs ≈ 1.6) and solves parity-2/3/4, but its *linear* memory
    is actually **lower** (MC ≈ 5.5 vs ≈ 7.6) and does not scale with ``n_virtual``.
    So it is a different operating point, not a strict improvement; genuinely beating
    the linear-MC ceiling needs a larger/better-tuned coupled network (future work).

    Per input step the virtual nodes update sequentially:
        V_i = a_in·mask_i·u + fb·x_{i-1},   x_{i-1} = previous virtual node (wrap →
              last node of the previous cycle),
        x_i = s∞(V_i) + (x_i^prev − s∞(V_i))·exp(−θ/tau(V_i)),
    with the same tanh nonlinearity / fading memory as the RTN node.
    """

    def __init__(self, n_virtual: int = 100, *, tau0: float = 1.0, Delta: float = 2.0,
                 vc0: float = 0.884, a_in: float = 0.5, fb_gain: float = 0.3,
                 theta: float = 1.0, bias_spread: float = 0.2, seed: Optional[int] = None):
        rng = np.random.default_rng(seed)
        self.n = int(n_virtual)
        self.mask = rng.choice(np.array([-1.0, 1.0]), size=self.n)
        # per-virtual-node bias offset breaks the odd-tanh symmetry (needed for even
        # parity / quadratic capacity), as in the spatial Reservoir.
        self.bias = rng.uniform(-bias_spread, bias_spread, size=self.n)
        self.tau0 = float(tau0)
        self.Delta = float(Delta)
        self.vc0 = float(vc0)
        self.a_in = float(a_in)
        self.fb_gain = float(fb_gain)
        self.theta = float(theta)

    def run(self, u, *, washout: int = 0, **_ignored) -> np.ndarray:
        """Drive the delay reservoir; extra kwargs (mode/n_replicas/seed) are accepted
        and ignored so it is drop-in compatible with the benchmark helpers."""
        u = np.asarray(u, float)
        T = u.shape[0]
        X = np.empty((T, self.n), float)
        x = np.zeros(self.n)
        lim = 0.98 * self.vc0
        for nstep in range(T):
            x_prev = x                      # states one delay-cycle ago
            new = np.empty(self.n)
            fb = x_prev[-1]                 # wrap-around delayed feedback
            for i in range(self.n):
                V = self.a_in * self.mask[i] * u[nstep] + self.fb_gain * fb + self.bias[i]
                V = lim if V > lim else (-lim if V < -lim else V)
                s_inf = np.tanh(self.Delta * V / self.vc0)
                tau = relaxation_time(V, tau0=self.tau0, Delta=self.Delta, Vc0=self.vc0)
                xi = s_inf + (x_prev[i] - s_inf) * np.exp(-self.theta / tau)
                new[i] = xi
                fb = xi                     # intra-cycle chain coupling
            x = new
            X[nstep] = x
        return X[washout:]
