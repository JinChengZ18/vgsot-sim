r"""
Experiment C (thesis §2.2.3.2) — Norm preservation & renormalisation
decorativeness.

Tests SC1 / C9: the Cayley step preserves |m|=1 unconditionally; therefore the
defensive per-step renormalisation in `switching_vector()`
(dynamic_switching_vector.py:208-211) is *decorative*, not load-bearing.

This script imports the *public* API only (no edits to src/). It has three
parts:

  Part A  — Algebraic identity on a (dt, |omega|) grid. Calls `cayley_step`
            directly (NO renorm) and compares the closed-form rotation against
            the matrix form  (I - s Omega)^-1 (I + s Omega)  via
            np.linalg.solve, wrapped in try/except to record singular/overflow
            corners. Two grids:
              * PHYSICAL  : dt in [1e-15, 1e-12], |omega| in [1e8, 1e12]
              * ADVERTISED: dt in [1e-15, 1e6 ], |omega| in [1   , 1e14]
            Records max ||m|-1| over each grid and the closed-vs-matrix
            mismatch (the |omega|->1e14 corner exposes the 'unconditional /
            any-dt' wording that must be softened).

  Part B  — LOAD-BEARING raw-trajectory test. Drives a RAW Cartesian state with
            _omega_from_state + cayley_step over 0.75 ns with NO renorm and NO
            theta/phi spherical round-trip. The shipped driver
            (time_series_cases.py:326-327) silently re-imposes |m|=1 by
            reconstructing m from (theta, phi) every step, which would mask any
            drift; here we keep the raw Cartesian vector so any norm growth is
            visible. Reports max_t ||m(t)|-1|.

  Part C  — Shipped-driver renorm ON/OFF. This ONLY probes pipeline
            INSENSITIVITY (the driver re-imposes |m|=1 via spherical
            reconstruction regardless of the renorm line), so it is explicitly
            NOT load-bearing. Kept as a documented sanity control.

Decisive criterion (from the experiment plan, §2.2.3.2 Exp C):
  * physical grid   max ||m|-1| <= 8.9e-16
  * raw trajectory  max ||m|-1| <  1e-12        => renorm decorative
  * advertised grid |omega|->1e14 closed-vs-matrix mismatch (matrix singular /
    overflow)                                    => 'unconditional / any-dt'
                                                    wording must be softened.

Run:  PYTHONPATH=src python scripts/02_integrator/test_norm_preservation.py
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.anisotropy import field
from vgsot_sim.dynamic_switching_vector import (
    _omega_from_state,
    cayley_step,
    switching_vector,
)

RESULT_DIR = Path(__file__).resolve().parents[2] / "result" / "sec_2_2_3_2" / "C"
RESULT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _skew(w):
    """Cross-product matrix [w]_x  (so [w]_x @ v == w x v)."""
    wx, wy, wz = float(w[0]), float(w[1]), float(w[2])
    return np.array([[0.0, -wz, wy],
                     [wz, 0.0, -wx],
                     [-wy, wx, 0.0]])


def cayley_matrix_step(m, omega, dt):
    """Matrix-form Cayley step  m_new = (I - s Omega)^-1 (I + s Omega) m.

    Returns (m_new, status) where status is "ok", "singular", or "overflow".
    Wrapped so the advertised-grid extreme corners (|omega|->1e14, dt->1e6)
    that drive (I - s Omega) singular / overflowing are recorded, not crashed.
    """
    s = 0.5 * dt
    Omega = _skew(omega)
    with warnings.catch_warnings():
        warnings.simplefilter("error", category=RuntimeWarning)
        try:
            A = np.eye(3) - s * Omega
            B = np.eye(3) + s * Omega
            rhs = B @ np.asarray(m, float)
            if not np.all(np.isfinite(A)) or not np.all(np.isfinite(rhs)):
                return None, "overflow"
            m_new = np.linalg.solve(A, rhs)
            if not np.all(np.isfinite(m_new)):
                return None, "overflow"
            return m_new, "ok"
        except np.linalg.LinAlgError:
            return None, "singular"
        except (RuntimeWarning, FloatingPointError, OverflowError):
            return None, "overflow"


def omega_at(m, cc, *, V_MTJ, I_SOT, R_MTJ, ESOT, ESTT, R_SOT_FL_DL,
             VNV, NON=0, h_th_ext=None, ENE=1, T=None):
    """Verified omega path (from harness): build H_eff then _omega_from_state."""
    m = np.asarray(m, float)
    theta = float(np.arccos(np.clip(m[2], -1.0, 1.0)))
    phi = float(np.arctan2(m[1], m[0]))
    H_eff, _ = field(theta, phi, V_MTJ, n=1, NON=NON, ENE=ENE, VNV=VNV,
                     constants=cc, demag_mode="ellipsoid",
                     h_th_ext=h_th_ext, T=T)
    H_eff = np.asarray(H_eff, float)
    Ms_use = cc.Ms
    I_MTJ = V_MTJ / R_MTJ if R_MTJ else 0.0
    J_STT = I_MTJ / cc.A1
    J_SOT = I_SOT / cc.A2
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
    H_DL_STT = ESTT * cc.h_bar * cc.P * J_STT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_DL_SOT = ESOT * cc.h_bar * cc.theta_SH * J_SOT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_FL_SOT = R_SOT_FL_DL * H_DL_SOT
    return _omega_from_state(
        m, H_eff, np.array([-1.0, 0.0, 0.0]),
        H_DL_SOT=H_DL_SOT, H_FL_SOT=H_FL_SOT,
        H_DL_STT=H_DL_STT, H_FL_STT=0.0,
        sigma_STT=np.array([0.0, 0.0, 1.0]),
        alpha=cc.alpha, gamma_red=gamma_red,
    )


# ---------------------------------------------------------------------------
# Part A — (dt, |omega|) grid: closed-form cayley_step (no renorm) vs matrix
# ---------------------------------------------------------------------------
def part_A_grid(name, dt_list, omega_mag_list, n_dirs=8, seed=0):
    """Sweep (dt, |omega|). For each cell, draw a few unit m and unit omega-hat,
    scale omega-hat to |omega|, run closed-form cayley_step (no renorm) and the
    matrix form. Record:
      * max ||m_new|-1| over the closed-form grid (norm preservation)
      * max closed-vs-matrix relative mismatch (where matrix is ok)
      * count of singular / overflow corners
    """
    rng = np.random.default_rng(seed)
    max_norm_err = 0.0
    max_norm_err_cell = None
    max_mismatch = 0.0
    max_mismatch_cell = None
    n_ok = n_singular = n_overflow = 0
    worst_corner_mag = None
    cell_records = []

    for dt in dt_list:
        for wmag in omega_mag_list:
            cell_norm_err = 0.0
            cell_mismatch = 0.0
            cell_status = "ok"
            for _ in range(n_dirs):
                # random unit m
                m = rng.normal(size=3)
                m = m / np.linalg.norm(m)
                # random unit omega direction, scaled to wmag
                wdir = rng.normal(size=3)
                wdir = wdir / np.linalg.norm(wdir)
                omega = wmag * wdir

                m_closed = cayley_step(m, omega, dt)
                ne = abs(np.linalg.norm(m_closed) - 1.0)
                if np.isfinite(ne):
                    cell_norm_err = max(cell_norm_err, ne)

                m_mat, status = cayley_matrix_step(m, omega, dt)
                if status == "ok":
                    n_ok += 1
                    denom = np.linalg.norm(m_closed)
                    if denom > 0 and np.all(np.isfinite(m_mat)):
                        mismatch = np.linalg.norm(m_mat - m_closed) / denom
                        if np.isfinite(mismatch):
                            cell_mismatch = max(cell_mismatch, mismatch)
                elif status == "singular":
                    n_singular += 1
                    cell_status = "singular"
                else:
                    n_overflow += 1
                    cell_status = "overflow"

            if cell_norm_err > max_norm_err:
                max_norm_err = cell_norm_err
                max_norm_err_cell = (dt, wmag)
            if cell_mismatch > max_mismatch:
                max_mismatch = cell_mismatch
                max_mismatch_cell = (dt, wmag)
            if cell_status != "ok":
                if worst_corner_mag is None or wmag > worst_corner_mag:
                    worst_corner_mag = wmag
            cell_records.append({
                "dt": float(dt), "omega_mag": float(wmag),
                "norm_err": float(cell_norm_err),
                "closed_vs_matrix_mismatch": float(cell_mismatch),
                "matrix_status": cell_status,
            })

    return {
        "name": name,
        "n_cells": len(dt_list) * len(omega_mag_list),
        "n_dirs_per_cell": n_dirs,
        "max_norm_err": float(max_norm_err),
        "max_norm_err_cell": [float(x) for x in max_norm_err_cell] if max_norm_err_cell else None,
        "max_closed_vs_matrix_mismatch": float(max_mismatch),
        "max_mismatch_cell": [float(x) for x in max_mismatch_cell] if max_mismatch_cell else None,
        "n_matrix_ok": n_ok,
        "n_matrix_singular": n_singular,
        "n_matrix_overflow": n_overflow,
        "worst_failing_omega_mag": (float(worst_corner_mag)
                                    if worst_corner_mag is not None else None),
        "cells": cell_records,
    }


# ---------------------------------------------------------------------------
# Part B — LOAD-BEARING raw Cartesian trajectory (no renorm, no round-trip)
# ---------------------------------------------------------------------------
def part_B_raw_trajectory(cc, i_sot_A, *, t_horizon=0.75e-9, dt=1e-12,
                          m0=None, ene=1, non=0, label=""):
    """Integrate a RAW Cartesian m with _omega_from_state + cayley_step.

    NO renorm. NO theta/phi spherical reconstruction. Records max_t ||m|-1|
    over the whole 0.75 ns drive. This is the sole load-bearing norm test: if
    |m| stays unit without any projection, the per-step renorm in
    switching_vector() is decorative.
    """
    if m0 is None:
        # start slightly off the +z pole, as in the AP->P drive
        th0 = 0.05
        m = np.array([np.sin(th0), 0.0, np.cos(th0)], float)
    else:
        m = np.asarray(m0, float).copy()

    n_steps = int(round(t_horizon / dt))
    max_norm_err = 0.0
    norm_trace = np.empty(n_steps + 1)
    norm_trace[0] = np.linalg.norm(m)
    mz_trace = np.empty(n_steps + 1)
    mz_trace[0] = m[2]

    for k in range(n_steps):
        omega = omega_at(
            m, cc, V_MTJ=0.0, I_SOT=i_sot_A, R_MTJ=0.0,
            ESOT=1, ESTT=0, R_SOT_FL_DL=0.83, VNV=0, NON=non, ENE=ene,
        )
        m = cayley_step(m, omega, dt)          # NO renorm, NO round-trip
        nrm = np.linalg.norm(m)
        norm_trace[k + 1] = nrm
        mz_trace[k + 1] = m[2]
        ne = abs(nrm - 1.0)
        if ne > max_norm_err:
            max_norm_err = ne

    return {
        "label": label,
        "i_sot_A": float(i_sot_A),
        "dt": float(dt),
        "n_steps": int(n_steps),
        "t_horizon_ns": float(t_horizon * 1e9),
        "max_norm_err": float(max_norm_err),
        "final_norm_err": float(abs(norm_trace[-1] - 1.0)),
        "mz_start": float(mz_trace[0]),
        "mz_end": float(mz_trace[-1]),
        "switched": bool(mz_trace[0] * mz_trace[-1] < 0),
        "norm_trace": norm_trace,
        "mz_trace": mz_trace,
    }


# ---------------------------------------------------------------------------
# Part C — shipped driver renorm ON vs OFF (pipeline insensitivity, NOT
#          load-bearing).  We cannot toggle the src renorm without editing
#          src/, so instead we run the shipped switching_vector() step (renorm
#          ON) and a renorm-OFF replica built from the SAME public sub-calls,
#          BOTH through the spherical (theta, phi) reconstruction that the
#          driver uses at time_series_cases.py:326-327. Because the driver
#          re-imposes |m|=1 every step via that reconstruction, the two
#          pipelines must agree to ~machine precision regardless of the renorm
#          line — which is exactly what makes this control NON-load-bearing.
# ---------------------------------------------------------------------------
def _shipped_step_renorm_off(m, cc, *, i_sot_A):
    """One switching_vector-equivalent step WITHOUT the final renorm, using the
    same public omega path the shipped step uses."""
    omega = omega_at(
        m, cc, V_MTJ=0.0, I_SOT=i_sot_A, R_MTJ=0.0,
        ESOT=1, ESTT=0, R_SOT_FL_DL=0.83, VNV=1, NON=0, ENE=1,
    )
    return cayley_step(m, omega, cc.t_step)


def part_C_pipeline_insensitivity(cc, i_sot_A, *, t_horizon=0.75e-9):
    """Run two driver-style pipelines (renorm ON via switching_vector, renorm
    OFF via replica) BOTH with the spherical (theta,phi) round-trip every step,
    and report the max divergence in m_z. Small divergence == renorm decorative
    *for the shipped pipeline*; but see the docstring: this only proves pipeline
    insensitivity, not norm preservation."""
    dt = cc.t_step
    n_steps = int(round(t_horizon / dt))
    th0 = 0.05
    m_on = np.array([np.sin(th0), 0.0, np.cos(th0)], float)
    m_off = m_on.copy()
    max_div = 0.0
    for _ in range(n_steps):
        # renorm ON: shipped step (normalises internally), then driver round-trip
        m_on_new = switching_vector(
            m_on, 0.0, i_sot_A, 0.0, 0, 1,
            VNV=1, NON=0, R_SOT_FL_DL=0.83,
            sigma_SH=np.array([-1.0, 0.0, 0.0]),
            constants=cc,
        )
        th = float(np.arccos(np.clip(m_on_new[2], -1.0, 1.0)))
        ph = float(np.arctan2(m_on_new[1], m_on_new[0]))
        m_on = np.array([np.sin(th) * np.cos(ph),
                         np.sin(th) * np.sin(ph), np.cos(th)])

        # renorm OFF replica: same omega path, NO renorm, then driver round-trip
        m_off_raw = _shipped_step_renorm_off(m_off, cc, i_sot_A=i_sot_A)
        th = float(np.arccos(np.clip(m_off_raw[2], -1.0, 1.0)))
        ph = float(np.arctan2(m_off_raw[1], m_off_raw[0]))
        m_off = np.array([np.sin(th) * np.cos(ph),
                          np.sin(th) * np.sin(ph), np.cos(th)])

        div = float(np.linalg.norm(m_on - m_off))
        if div > max_div:
            max_div = div
    return {
        "i_sot_A": float(i_sot_A),
        "t_horizon_ns": float(t_horizon * 1e9),
        "n_steps": int(n_steps),
        "max_pipeline_divergence_renorm_on_off": float(max_div),
        "note": ("pipeline INSENSITIVITY only; both pipelines re-impose |m|=1 "
                 "via spherical reconstruction (time_series_cases.py:326-327), "
                 "so this is NOT a norm-preservation test"),
    }


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    cc = PhysicalConstantsConfig()
    cc.t_step = 1e-12

    print("=" * 76)
    print("Experiment C — norm preservation & renormalisation decorativeness")
    print("=" * 76)

    # ---- Part A: physical grid ----------------------------------------
    dt_phys = np.geomspace(1e-15, 1e-12, 7)
    omega_phys = np.geomspace(1e8, 1e12, 9)
    A_phys = part_A_grid("physical", dt_phys, omega_phys, n_dirs=8, seed=1)

    # ---- Part A: advertised grid --------------------------------------
    dt_adv = np.geomspace(1e-15, 1e6, 12)
    omega_adv = np.geomspace(1.0, 1e14, 15)
    A_adv = part_A_grid("advertised", dt_adv, omega_adv, n_dirs=8, seed=2)

    print("\n[Part A] (dt, |omega|) grid — cayley_step (no renorm) vs matrix")
    for tag, A in (("PHYSICAL", A_phys), ("ADVERTISED", A_adv)):
        print(f"  {tag:11s}: cells={A['n_cells']:3d}  "
              f"max||m|-1|={A['max_norm_err']:.3e}  "
              f"closed-vs-matrix={A['max_closed_vs_matrix_mismatch']:.3e}  "
              f"matrix(ok/sing/ovf)={A['n_matrix_ok']}/"
              f"{A['n_matrix_singular']}/{A['n_matrix_overflow']}")
    print(f"  ADVERTISED worst failing |omega| corner: "
          f"{A_adv['worst_failing_omega_mag']}")

    # ---- Part B: LOAD-BEARING raw trajectory --------------------------
    # Drive currents spanning sub-/super-threshold (AP->P, negative I_SOT).
    i_list_uA = [-600, -1300, -2000]
    B_runs = []
    for i_uA in i_list_uA:
        B = part_B_raw_trajectory(cc, i_uA * 1e-6, t_horizon=0.75e-9,
                                  dt=1e-12, label=f"I_SOT={i_uA}uA")
        B_runs.append(B)
    B_max_over_runs = max(b["max_norm_err"] for b in B_runs)

    print("\n[Part B] LOAD-BEARING raw Cartesian trajectory (no renorm, no "
          "round-trip), 0.75 ns")
    for b in B_runs:
        print(f"  {b['label']:14s}: steps={b['n_steps']}  "
              f"max||m|-1|={b['max_norm_err']:.3e}  "
              f"mz {b['mz_start']:+.3f}->{b['mz_end']:+.3f}  "
              f"switched={b['switched']}")
    print(f"  -> max over runs: {B_max_over_runs:.3e}")

    # ---- Part C: pipeline insensitivity (NOT load-bearing) ------------
    C = part_C_pipeline_insensitivity(cc, -1300e-6, t_horizon=0.75e-9)
    print("\n[Part C] shipped-driver renorm ON/OFF (pipeline insensitivity, "
          "NOT load-bearing)")
    print(f"  max pipeline divergence (renorm on vs off): "
          f"{C['max_pipeline_divergence_renorm_on_off']:.3e}")

    # ---- criteria -----------------------------------------------------
    crit_phys = A_phys["max_norm_err"] <= 8.9e-16
    crit_raw = B_max_over_runs < 1e-12
    crit_advertised_softening = (A_adv["n_matrix_singular"]
                                 + A_adv["n_matrix_overflow"]) > 0

    print("\n[Criteria]")
    print(f"  physical grid max||m|-1| <= 8.9e-16 : "
          f"{A_phys['max_norm_err']:.3e}  -> {crit_phys}")
    print(f"  raw trajectory max||m|-1| < 1e-12   : "
          f"{B_max_over_runs:.3e}  -> {crit_raw}  "
          f"(=> renorm decorative)")
    print(f"  advertised |omega|->1e14 matrix fails: "
          f"sing+ovf={A_adv['n_matrix_singular'] + A_adv['n_matrix_overflow']}"
          f"  -> {crit_advertised_softening}  (=> soften 'unconditional/any-dt')")

    # ---- figures ------------------------------------------------------
    _make_figures(A_phys, A_adv, B_runs)

    # ---- json ---------------------------------------------------------
    out = {
        "experiment": "sec_2_2_3_2_C_norm_preservation",
        "env": {"note": "PYTHONPATH=src python; numpy matrix form via solve"},
        "config": {"t_step": cc.t_step, "theta_SH": cc.theta_SH,
                   "alpha": cc.alpha, "Ms": cc.Ms, "gamma": cc.gamma},
        "partA_physical": _strip_cells(A_phys),
        "partA_advertised": _strip_cells(A_adv),
        "partB_raw_trajectory": [
            {k: v for k, v in b.items()
             if k not in ("norm_trace", "mz_trace")} for b in B_runs
        ],
        "partB_max_norm_err_over_runs": float(B_max_over_runs),
        "partC_pipeline_insensitivity": C,
        "criteria": {
            "physical_grid_le_8p9e-16": bool(crit_phys),
            "physical_grid_value": float(A_phys["max_norm_err"]),
            "raw_traj_lt_1e-12_renorm_decorative": bool(crit_raw),
            "raw_traj_value": float(B_max_over_runs),
            "advertised_matrix_fails_soften_wording": bool(crit_advertised_softening),
            "advertised_n_singular_plus_overflow": int(
                A_adv["n_matrix_singular"] + A_adv["n_matrix_overflow"]),
        },
        "verdict": {
            "SC1_norm_preserved": bool(crit_phys and crit_raw),
            "renorm_decorative": bool(crit_raw),
            "wording_must_soften": bool(crit_advertised_softening),
        },
    }
    json_path = RESULT_DIR / "norm_preservation_results.json"
    json_path.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {json_path}")

    return out


def _strip_cells(A):
    """Keep summary + full cell table (cells are small)."""
    return A


def _make_figures(A_phys, A_adv, B_runs):
    # Figure 1: (dt,|omega|) grids — norm error + matrix status heatmaps
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    for col, (tag, A) in enumerate((("Physical grid", A_phys),
                                    ("Advertised grid", A_adv))):
        cells = A["cells"]
        dts = sorted(set(c["dt"] for c in cells))
        ws = sorted(set(c["omega_mag"] for c in cells))
        di = {d: i for i, d in enumerate(dts)}
        wi = {w: j for j, w in enumerate(ws)}

        norm_grid = np.full((len(ws), len(dts)), np.nan)
        status_grid = np.zeros((len(ws), len(dts)))  # 0 ok, 1 singular, 2 overflow
        for c in cells:
            i = di[c["dt"]]; j = wi[c["omega_mag"]]
            ne = c["norm_err"]
            norm_grid[j, i] = np.log10(ne) if ne > 0 else -17
            status_grid[j, i] = {"ok": 0, "singular": 1, "overflow": 2}[c["matrix_status"]]

        ax = axes[0, col]
        im = ax.pcolormesh(np.arange(len(dts) + 1), np.arange(len(ws) + 1),
                           norm_grid, cmap="viridis", shading="flat")
        ax.set_title(f"{tag}\nlog10 ||m|-1|  (closed-form cayley_step, no renorm)")
        ax.set_xlabel("dt index (small->large)")
        ax.set_ylabel("|omega| index (small->large)")
        fig.colorbar(im, ax=ax, fraction=0.046)

        ax = axes[1, col]
        im = ax.pcolormesh(np.arange(len(dts) + 1), np.arange(len(ws) + 1),
                           status_grid, cmap="RdYlGn_r", vmin=0, vmax=2,
                           shading="flat")
        ax.set_title("matrix form status (0 ok / 1 singular / 2 overflow)")
        ax.set_xlabel("dt index (small->large)")
        ax.set_ylabel("|omega| index (small->large)")
        fig.colorbar(im, ax=ax, fraction=0.046, ticks=[0, 1, 2])

    fig.suptitle("Exp C Part A — Cayley norm preservation & matrix-form "
                 "breakdown on (dt, |omega|) grids", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    f1 = RESULT_DIR / "partA_grid_norm_preservation.png"
    fig.savefig(f1, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {f1}")

    # Figure 2: raw trajectory ||m(t)|-1| flat line + m_z(t)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for b in B_runs:
        t = np.arange(len(b["norm_trace"])) * b["dt"] * 1e9
        ax1.plot(t, np.abs(b["norm_trace"] - 1.0), lw=1.4, label=b["label"])
        ax2.plot(t, b["mz_trace"], lw=1.4, label=b["label"])
    ax1.set_yscale("log")
    ax1.axhline(1e-12, color="crimson", ls="--", lw=1.0, label="1e-12 criterion")
    ax1.set_ylabel(r"$\,||\mathbf{m}(t)|-1|$  (raw, no renorm)")
    ax1.set_title("Exp C Part B (load-bearing) — raw Cartesian trajectory, "
                  "0.75 ns, no renorm / no round-trip")
    ax1.legend(fontsize=9)
    ax1.grid(True, which="both", alpha=0.3)
    ax2.set_ylabel(r"$m_z(t)$")
    ax2.set_xlabel("Time (ns)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)
    fig.tight_layout()
    f2 = RESULT_DIR / "partB_raw_trajectory_norm.png"
    fig.savefig(f2, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {f2}")


if __name__ == "__main__":
    main()
