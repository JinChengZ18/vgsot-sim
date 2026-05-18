"""Binomial sampling reliability + MC sampling-size sensitivity.

Two related questions answered here:

1. **Coverage of a sample mean**: given a true success probability p and K
   independent Bernoulli trials, what is P(|k/K − p| < ε)?  Exact closed
   form via the scipy binomial CDF + asymptotic CLT comparison are both
   provided.

2. **Sampling-size sensitivity for D2D MC**: how many MC samples N are
   needed before the estimator F_hat(N, CV_Δ) of the D2D broadening
   factor is within a target tolerance of its asymptotic value?

Used by both 08 scripts (`hw_sampling_reliability` and
`sampling_sensitivity_sim`); centralising it here lets the chapter
numbers be regenerated from a single source.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

from .variability import wafer_average_psw, fit_sigmoid_to_average


# ── 1. Binomial coverage / K_req ──────────────────────────────────────────

def exact_coverage(K: int, p: float, eps: float, strict: bool = True) -> float:
    """Exact P(|k/K − p| < eps) (strict) or ≤ eps (closed) for k ~ Bin(K, p).

    Returns a single probability. Clamps the success-count window to
    [0, K] for edge cases p ≈ 0 or p ≈ 1.
    """
    K = int(K)
    if strict:
        k_lo = int(np.floor(K * (p - eps))) + 1
        k_hi = int(np.ceil(K * (p + eps))) - 1
    else:
        k_lo = int(np.ceil(K * (p - eps)))
        k_hi = int(np.floor(K * (p + eps)))
    k_lo = max(k_lo, 0)
    k_hi = min(k_hi, K)
    if k_lo > k_hi:
        return 0.0
    return float(stats.binom.cdf(k_hi, K, p) - stats.binom.cdf(k_lo - 1, K, p))


def clt_K_required(p: float, eps: float, confidence: float = 0.95) -> float:
    """CLT-asymptotic K such that P(|k/K - p| < eps) ≥ confidence."""
    z = stats.norm.ppf(0.5 + confidence / 2.0)
    return float(z * z * p * (1 - p) / (eps * eps))


def K_required_exact(p: float, eps: float, confidence: float = 0.95,
                     K_max: int = 200_000) -> int:
    """Smallest K (1..K_max) such that exact_coverage(K, p, eps) ≥ confidence.

    Uses monotone bisection. Raises if K_max is insufficient. Note that
    binomial coverage is only weakly monotone in K (lattice effects cause
    small dips), so `K_required_persistent` is usually preferred for
    publishing protocol thresholds.
    """
    lo, hi = 1, K_max
    if exact_coverage(hi, p, eps) < confidence:
        raise ValueError(
            f"K_max={K_max} too small for p={p}, eps={eps}, conf={confidence}"
        )
    while lo < hi:
        mid = (lo + hi) // 2
        if exact_coverage(mid, p, eps) >= confidence:
            hi = mid
        else:
            lo = mid + 1
    return int(lo)


def K_required_persistent(p: float, eps: float, confidence: float = 0.95) -> int:
    """Smallest K* such that coverage(K') ≥ confidence for **all** K' ≥ K*.

    Uses the asymptotic CLT estimate × 4 as an upper-bound K_hi, doubles
    until coverage(K_hi) ≥ confidence, then enumerates 1..K_hi and applies
    a reverse-cumulative minimum to find the lower edge of the persistent
    tail. Strictly more conservative than `K_required_exact`.
    """
    z = stats.norm.ppf((1 + confidence) / 2)
    K_hi = max(int(np.ceil(z * z * p * (1 - p) / (eps * eps)) * 4), 200)
    while exact_coverage(K_hi, p, eps) < confidence:
        K_hi *= 2
    covs = np.array([exact_coverage(K, p, eps) for K in range(1, K_hi + 1)])
    min_cov = np.minimum.accumulate(covs[::-1])[::-1]
    ok = np.where(min_cov >= confidence)[0]
    if len(ok) == 0:
        raise RuntimeError("Persistence target unreachable")
    return int(ok[0] + 1)


def mc_coverage(K: int, p: float, eps: float, M: int, rng) -> float:
    """Empirical coverage from M Bin(K, p) replicates.

    Uses integer comparison `k_lo ≤ k ≤ k_hi` rather than `|k/K - p| < eps`
    to avoid the float-precision edge-bias at p = 0.5, eps = 0.10.
    """
    k = rng.binomial(n=K, p=p, size=M)
    k_lo = int(np.floor(K * (p - eps))) + 1
    k_hi = int(np.ceil(K * (p + eps))) - 1
    return float(np.mean((k >= k_lo) & (k <= k_hi)))


def K_required_persistent_mc(p: float, eps: float, M: int, rng,
                             confidence: float = 0.95,
                             K_hi_cap: int = 10_000_000) -> int | None:
    """MC analogue of `K_required_persistent` using empirical coverage."""
    z = stats.norm.ppf((1 + confidence) / 2)
    K_hi = max(int(np.ceil(z * z * p * (1 - p) / (eps * eps)) * 4), 200)
    while mc_coverage(K_hi, p, eps, M, rng) < confidence + 0.01:
        K_hi *= 2
        if K_hi > K_hi_cap:
            return None
    covs = np.array([mc_coverage(K, p, eps, M, rng) for K in range(1, K_hi + 1)])
    min_cov = np.minimum.accumulate(covs[::-1])[::-1]
    ok = np.where(min_cov >= confidence)[0]
    return int(ok[0] + 1) if len(ok) else None


# ── 2. Sampling-size sensitivity for D2D MC ───────────────────────────────

def mc_beta_eff(N: int, cv: float, Delta: float, Vc0: float,
                tw_ns: float, tau0_ns: float = 1.0, rng=None,
                V_dense=None) -> float:
    """Single-CV Monte-Carlo β_eff estimator.

    Draws N Gaussian Δ samples, computes the wafer-averaged P_sw(V),
    fits a Sigmoid, returns its slope β. Used inside `f_curve_sweep` and
    `recommended_N`.
    """
    if V_dense is None:
        V_dense = np.linspace(0.60, 1.10, 600)
    P_mean = wafer_average_psw(V_dense, cv, Delta, Vc0,
                               tw_ns=tw_ns, tau0_ns=tau0_ns,
                               N_samples=N, rng=rng)
    _, beta = fit_sigmoid_to_average(V_dense, P_mean)
    return beta


def f_ref(cv: float, Delta: float, Vc0: float, tw_ns: float,
          tau0_ns: float = 1.0,
          N_ref: int = 50_000, n_replicas: int = 5,
          base_seed: int = 1234) -> tuple[float, float]:
    """Reference F = β(CV)/β(0) at large N, with `n_replicas` independent runs.

    Returns (mean, std) of F across replicas. Each replica uses a distinct
    seed derived from `base_seed + cv*10000 + r*7919` to break the seed
    sharing across CVs that occurs when one seed is reused.
    """
    base = int(round(base_seed + cv * 10000))
    rng0 = np.random.default_rng(seed=base)
    beta_at_zero = mc_beta_eff(N_ref, 0.0, Delta, Vc0, tw_ns, tau0_ns, rng=rng0)
    replicas = []
    for r in range(n_replicas):
        rng = np.random.default_rng(seed=base + r * 7919)
        beta = mc_beta_eff(N_ref, cv, Delta, Vc0, tw_ns, tau0_ns, rng=rng)
        replicas.append(beta / beta_at_zero)
    return float(np.mean(replicas)), float(np.std(replicas, ddof=1))


def recommended_N(cv: float, Delta: float, Vc0: float, tw_ns: float,
                  tau0_ns: float = 1.0,
                  N_grid=(100, 200, 500, 1000, 2000, 5000, 10_000),
                  R: int = 40, tol: float = 0.02, confidence: float = 0.95,
                  base_seed: int = 2026):
    """For each N in `N_grid`, run R independent F estimates and return
    (N, mean F̂, std F̂, fraction within `tol` relative error). The
    recommended N is the smallest N at which that fraction ≥ confidence.
    """
    master = np.random.default_rng(seed=base_seed + int(round(cv * 10000)))
    F_ref_val, _ = f_ref(cv, Delta, Vc0, tw_ns, tau0_ns,
                          N_ref=50_000, n_replicas=3,
                          base_seed=base_seed + 999)
    rows = []
    rec_N = None
    for N in N_grid:
        F_rep = np.empty(R)
        for r in range(R):
            seed = int(master.integers(0, 2 ** 31 - 1))
            rng = np.random.default_rng(seed=seed)
            beta = mc_beta_eff(int(N), cv, Delta, Vc0, tw_ns, tau0_ns, rng=rng)
            rng0 = np.random.default_rng(seed=seed + 1)
            beta0 = mc_beta_eff(int(N), 0.0, Delta, Vc0, tw_ns, tau0_ns, rng=rng0)
            F_rep[r] = beta / beta0
        rel_err = np.abs(F_rep - F_ref_val) / F_ref_val
        within = float((rel_err < tol).mean())
        rows.append((int(N), float(F_rep.mean()), float(F_rep.std(ddof=1)), within))
        if rec_N is None and within >= confidence:
            rec_N = int(N)
    return F_ref_val, rows, rec_N
