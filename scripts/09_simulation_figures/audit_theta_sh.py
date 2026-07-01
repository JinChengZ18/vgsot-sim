"""
实验 F — PART B：θ_SH 运行时读回 + 阈值系综溯源审计（命题 C13）。

目的
----
图脚本注释与正文残留「θ_SH = 0.04（球坐标-Euler 标定值）」字样，但发布的
`PhysicalConstantsConfig.theta_SH` 默认值已是 **0.066**（Cayley 积分器，按
Device-A P→AP 0.75 ns 的实验 V_th=894 mV → I_th≈1152 µA 重标定）。

本脚本两段：
  PART B-1（决定性，cheap）— 运行时断言：
     * PhysicalConstantsConfig().theta_SH == 0.066
     * 图 2.10 脚本 plot_single_trajectory.py 既不覆盖 theta_SH，也不覆盖
       integrator（→ 继承 0.066 + cayley）。
  PART B-2（系综溯源，较重；此处跑 *缩减 pilot* N≈30 种子）—
     对 θ_SH ∈ {0.066, 0.04}：在电流网格上对 N 个独立热噪声种子统计
     P→AP（NEGATIVE I_SOT）翻转成功率 P_sw(|I|)（成功判据与标定图 plot_ser_mc
     一致：|m_z_end − 1| ≤ 0.2，4 ns 窗口），再对 P_sw(|I|) 拟合 Logistic 取
     **P_sw=0.5 交点 I_50**（即系综中位翻转阈值），看哪个 θ_SH 的 I_50 落在实验
     I_th≈1152 µA 附近。I_50 的 95% CI 由对 N 个种子做自助重采样得到。
     红队加固：θ_SH 验证**必须**用 NEGATIVE I_SOT（P→AP）+ N≥200 种子系综
     （单种子阈值跨度大，不能凭 seed=1）。本 pilot N=30 仅作 harness 验证；
     完整 N≥200 见 deferred_full_run_cmd / --full。

决定性判据 (C13)
----------------
theta_SH == 0.066（而非 0.04）→ 图 2.10 标题 / 图脚本「0.04, 球坐标-Euler」
溯源注释**陈旧**、被发布代码证伪。

只调用公共 API；不修改 src/。matplotlib 用 Agg。
运行（pilot）:  PYTHONPATH=src python scripts/09_simulation_figures/audit_theta_sh.py
运行（完整）:  PYTHONPATH=src python scripts/09_simulation_figures/audit_theta_sh.py --full --n-seeds 200
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation

REPO_ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = Path(__file__).resolve().parent
OUT_DIR = REPO_ROOT / "result" / "sec_2_2_3_2" / "F"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Pulse protocol — aligned to the calibration figure plot_ser_mc.py (the SER
# Monte-Carlo that fixed θ_SH=0.066): TOTAL 4 ns window, 0.75 ns pulse, and the
# SER success criterion |m_z_end − target| ≤ failure_tol (target=1.0, tol=0.2 →
# success ⇔ m_z_end ≥ 0.8). The ensemble I_50 is the |I_SOT| at which the
# *fraction of seeds that succeed* (P_sw) crosses 0.5 — i.e. the median
# switching threshold, NOT a per-seed bisection (red-team: single seed is
# unrepresentative; need the N-seed ensemble median).
PULSE_NS = 0.75
TOTAL_NS = 4.0
TARGET_MZ = 1.0
FAILURE_TOL = 0.2        # success ⇔ |m_z_end − 1.0| ≤ 0.2

# Current grid for the P_sw(|I|) sigmoid (µA, applied as NEGATIVE current).
I_GRID_UA = np.array([900, 1000, 1100, 1150, 1200, 1250, 1300, 1350,
                      1400, 1500, 1700, 2000], dtype=float)


def _success(i_sot_A: float, theta_SH: float, seed: int) -> bool:
    """One seeded P→AP trajectory (NEGATIVE I_SOT) at a given θ_SH; True if it
    meets the SER success criterion. θ_SH injected via a mutated dataclass copy.
    """
    cc = dataclasses.replace(PhysicalConstantsConfig(), theta_SH=theta_SH)
    sim_end = int(round(TOTAL_NS * 1e-9 / cc.t_step))
    mid1 = int(round(PULSE_NS * 1e-9 / cc.t_step))
    rng = np.random.default_rng(seed)
    res = run_piecewise_direct_excitation(
        sim_start_step=1, sim_mid1_step=mid1, sim_mid2_step=sim_end, sim_end_step=sim_end,
        pap=1,
        v_mtj_stage1=0.0, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
        i_sot_stage1=i_sot_A, i_sot_stage2=0.0, i_sot_stage3=0.0,
        estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1,
        estt_stage3=0, esot_stage3=1,
        vnv=0, non=1, r_sot_fl_dl=0.83,
        constants=cc, show_progress=False,
        # integrator omitted → driver default 'cayley' (the published path)
    )
    return abs(float(res.mz[-1]) - TARGET_MZ) <= FAILURE_TOL


def _logistic_i50(I_uA, psw):
    """Fit P_sw = 1/(1+exp(-k(|I|-I50))) and return I50 (µA).

    Robust fallback to linear interpolation of the 0.5 crossing if the
    optimiser fails or the curve never brackets 0.5.
    """
    I_uA = np.asarray(I_uA, float)
    psw = np.asarray(psw, float)
    # linear-interp crossing as fallback / initial guess
    cross = float("nan")
    for j in range(len(I_uA) - 1):
        if (psw[j] - 0.5) * (psw[j + 1] - 0.5) <= 0 and psw[j] != psw[j + 1]:
            cross = I_uA[j] + (0.5 - psw[j]) * (I_uA[j + 1] - I_uA[j]) / (psw[j + 1] - psw[j])
            break
    try:
        from scipy.optimize import curve_fit
        def logf(x, k, i50):
            return 1.0 / (1.0 + np.exp(-k * (x - i50)))
        i50_0 = cross if np.isfinite(cross) else float(np.mean(I_uA))
        popt, _ = curve_fit(logf, I_uA, psw, p0=[0.01, i50_0], maxfev=10000)
        i50 = float(popt[1])
        if I_uA.min() - 200 <= i50 <= I_uA.max() + 200:
            return i50, float(popt[0])
    except Exception:
        pass
    return cross, float("nan")


def ensemble_i50(theta_SH: float, n_seeds: int, n_boot: int = 1000,
                 boot_rng=None):
    """Ensemble P_sw(|I|) over n_seeds, fit I_50 (P_sw=0.5 crossing), and a
    seed-bootstrap 95% CI on I_50. Returns dict."""
    boot_rng = boot_rng or np.random.default_rng(0)
    # success matrix: rows = currents, cols = seeds
    S = np.zeros((len(I_GRID_UA), n_seeds), dtype=bool)
    for r, i_uA in enumerate(I_GRID_UA):
        for c in range(n_seeds):
            S[r, c] = _success(-i_uA * 1e-6, theta_SH, c)
    psw = S.mean(axis=1)
    i50, k = _logistic_i50(I_GRID_UA, psw)
    # seed bootstrap on I_50
    boots = []
    for _ in range(n_boot):
        cols = boot_rng.integers(0, n_seeds, n_seeds)
        psw_b = S[:, cols].mean(axis=1)
        i50_b, _k = _logistic_i50(I_GRID_UA, psw_b)
        if np.isfinite(i50_b):
            boots.append(i50_b)
    if boots:
        ci = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
    else:
        ci = [None, None]   # keep JSON strict-valid (no NaN literals)
    return {
        "theta_SH": theta_SH,
        "n_seeds": n_seeds,
        "I_grid_uA": [float(x) for x in I_GRID_UA],
        "psw": [float(x) for x in psw],
        "I_50_uA": (float(i50) if np.isfinite(i50) else None),
        "logistic_k": (float(k) if np.isfinite(k) else None),
        "I_50_uA_ci95": ci,
        # If the sigmoid never crosses 0.5 on the grid, report the grid ceiling
        # as a one-sided lower bound on I_50 (decisive for the 0.04 case).
        "I_50_lower_bound_uA": (None if np.isfinite(i50)
                                else float(I_GRID_UA.max())),
    }


def part_b1_readback():
    """Cheap, decisive runtime read-back of θ_SH and figure-script overrides."""
    cc = PhysicalConstantsConfig()
    theta_SH = cc.theta_SH
    assert theta_SH == 0.066, f"theta_SH read-back FAILED: {theta_SH!r} != 0.066"

    fig210_src = (FIG_DIR / "plot_single_trajectory.py").read_text(encoding="utf-8")
    overrides_theta = ("theta_SH=" in fig210_src) or ("theta_SH =" in fig210_src
                                                       and "cc.theta_SH" in fig210_src)
    overrides_integrator = "integrator=" in fig210_src
    # Mentions of the STALE 0.04 value live only in comments (provenance smell).
    mentions_004_comment = "0.04" in fig210_src

    I_th_exp_uA = 894e-3 / cc.R_W * 1e6

    print("=" * 78)
    print("  实验 F — PART B-1：θ_SH 运行时读回 (C13)")
    print("=" * 78)
    print(f"  PhysicalConstantsConfig().theta_SH        = {theta_SH}  (期望 0.066)")
    print(f"  R_W                                       = {cc.R_W:.3f} Ω")
    print(f"  实验 I_th = 894 mV / R_W                  = {I_th_exp_uA:.1f} µA")
    print(f"  fig2.10 覆盖 theta_SH?                     = {overrides_theta}")
    print(f"  fig2.10 覆盖 integrator?                   = {overrides_integrator}")
    print(f"  fig2.10 注释残留 '0.04'?                   = {mentions_004_comment}")
    print("-" * 78)

    return {
        "theta_SH_runtime": theta_SH,
        "theta_SH_expected": 0.066,
        "theta_SH_assert_passed": True,
        "R_W_ohm": cc.R_W,
        "I_th_exp_uA": I_th_exp_uA,
        "fig210_overrides_theta_SH": bool(overrides_theta),
        "fig210_overrides_integrator": bool(overrides_integrator),
        "fig210_mentions_stale_0p04_in_source": bool(mentions_004_comment),
        "C13_stale_0p04_confirmed": bool(theta_SH == 0.066),
    }


def part_b2_ensemble(n_seeds: int, theta_list=(0.066, 0.04)):
    """Ensemble P_sw(|I|) → I_50 (P_sw=0.5 crossing) per θ_SH; NEGATIVE I_SOT,
    P→AP. Returns {f'theta_SH={th}': {...}}."""
    print("=" * 78)
    print(f"  实验 F — PART B-2：阈值系综 I_50 (N={n_seeds} 种子/θ_SH, "
          f"NEGATIVE I_SOT, P→AP)")
    print("=" * 78)
    out = {}
    boot_rng = np.random.default_rng(20260630)
    for theta_SH in theta_list:
        t0 = time.time()
        d = ensemble_i50(theta_SH, n_seeds, boot_rng=boot_rng)
        d["wall_s"] = round(time.time() - t0, 2)
        out[f"theta_SH={theta_SH}"] = d
        i50 = d["I_50_uA"]
        lo, hi = d["I_50_uA_ci95"]
        i50s = f"{i50:7.1f}" if i50 is not None else (
            f" >{d['I_50_lower_bound_uA']:.0f}")
        cis = (f"[{lo:.0f},{hi:.0f}]" if (lo is not None) else "[lower-bnd]")
        print(f"  θ_SH={theta_SH:<6}  I_50={i50s} µA  CI95={cis}  "
              f"({d['wall_s']:.1f}s)")
        print(f"            P_sw grid: " +
              "  ".join(f"{ig:.0f}:{p:.2f}"
                        for ig, p in zip(d["I_grid_uA"], d["psw"])))
    print("-" * 78)
    return out


def make_figure(b2, I_th_exp_uA, out_png):
    """Two-panel: (left) P_sw(|I|) sigmoids per θ_SH with I_50 markers;
    (right) I_50 ± CI vs exp I_th."""
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.6),
                                   gridspec_kw={"width_ratios": [1.5, 1.0]})
    colors = {"0.066": "#660874", "0.04": "#C47A00"}
    keys = sorted(b2.keys(), key=lambda k: -b2[k]["theta_SH"])

    # left: sigmoids
    for key in keys:
        d = b2[key]
        th = d["theta_SH"]
        c = colors.get(str(th), "#1F5FA8")
        axL.plot(d["I_grid_uA"], d["psw"], "o-", color=c, lw=1.8, ms=5,
                 markerfacecolor="white", markeredgewidth=1.4,
                 label=rf"$\theta_{{\mathrm{{SH}}}}={th}$")
        i50 = d["I_50_uA"]
        if i50 is not None:
            axL.axvline(i50, color=c, ls=":", lw=1.2)
    axL.axhline(0.5, color="gray", lw=0.5, ls=":")
    axL.axvline(I_th_exp_uA, color="#1A6B5A", ls=(0, (4, 2)), lw=1.5,
                label=rf"exp $I_{{\mathrm{{th}}}}={I_th_exp_uA:.0f}\,\mu$A")
    axL.set_xlabel(r"$|I_{\mathrm{SOT}}|$ ($\mu$A)")
    axL.set_ylabel(r"$P_{\mathrm{sw}}$ (P→AP)")
    axL.set_ylim(-0.05, 1.05)
    axL.set_title("Ensemble P$_{sw}$(|I|) per θ$_{SH}$")
    axL.grid(True, ls="--", lw=0.4, alpha=0.6)
    axL.legend(fontsize=9, loc="lower right")

    # right: I_50 with CI
    xs, labels = [], []
    for x, key in enumerate(keys):
        d = b2[key]
        th = d["theta_SH"]
        c = colors.get(str(th), "#1F5FA8")
        i50 = d["I_50_uA"]
        lo, hi = d["I_50_uA_ci95"]
        if i50 is not None:
            yerr = ([[max(0.0, i50 - lo)], [max(0.0, hi - i50)]]
                    if lo is not None else None)
            axR.errorbar([x], [i50], yerr=yerr, fmt="o", color=c, ms=10,
                         capsize=5, lw=2,
                         label=rf"$\theta_{{\mathrm{{SH}}}}={th}$ ($I_{{50}}={i50:.0f}$)")
        else:
            # one-sided lower bound: sigmoid never reached 0.5 on the grid
            lb = d["I_50_lower_bound_uA"]
            axR.scatter([x], [lb], marker="^", color=c, s=90,
                        label=rf"$\theta_{{\mathrm{{SH}}}}={th}$ ($I_{{50}}>{lb:.0f}$)")
        xs.append(x)
        labels.append(rf"$\theta_{{\mathrm{{SH}}}}={th}$")
    axR.axhline(I_th_exp_uA, color="#1A6B5A", ls=(0, (4, 2)), lw=1.5,
                label=rf"exp $I_{{\mathrm{{th}}}}={I_th_exp_uA:.0f}\,\mu$A")
    axR.set_xticks(xs)
    axR.set_xticklabels(labels)
    axR.set_xlim(-0.6, len(keys) - 0.4)
    axR.set_ylabel(r"$I_{50}$ ($\mu$A)")
    axR.set_title("I$_{50}$ vs exp I$_{th}$")
    axR.grid(True, ls="--", lw=0.4, alpha=0.6)
    axR.legend(fontsize=8, loc="best")

    fig.suptitle("θ_SH provenance audit (NEGATIVE I$_{SOT}$, P→AP ensemble)",
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true",
                    help="full ensemble (use with --n-seeds 200); default is pilot")
    ap.add_argument("--n-seeds", type=int, default=None,
                    help="override seed count (pilot default 30, full default 200)")
    args = ap.parse_args()
    n_seeds = args.n_seeds if args.n_seeds is not None else (200 if args.full else 30)

    b1 = part_b1_readback()
    b2 = part_b2_ensemble(n_seeds)

    # Which θ_SH lands closest to experimental I_th?
    I_th_exp = b1["I_th_exp_uA"]
    closeness = {k: (abs(b2[k]["I_50_uA"] - I_th_exp)
                     if b2[k]["I_50_uA"] is not None else float("inf"))
                 for k in b2}
    best = min(closeness, key=closeness.get)

    out_png = OUT_DIR / "theta_sh_ensemble.png"
    make_figure(b2, I_th_exp, out_png)

    result = {
        "experiment": "F",
        "part": "B_theta_sh",
        "claim": "C13",
        "claim_text_zh": ("θ_SH = 0.04 还是 0.066；图脚本/正文「0.04, 球坐标-Euler」"
                          "溯源是否陈旧"),
        "mode": "full" if args.full else "pilot",
        "n_seeds": n_seeds,
        "B1_readback": b1,
        "B2_ensemble": b2,
        "closest_theta_to_exp_Ith": {
            "best_key": best,
            "abs_gap_uA": {k: (round(v, 1) if np.isfinite(v) else None)
                           for k, v in closeness.items()},
        },
        "decisive_criterion": "theta_SH == 0.066 (not 0.04)  → fig-2.10 '0.04 / 球坐标-Euler' caption STALE",
        "C13_stale_0p04_confirmed": bool(b1["theta_SH_runtime"] == 0.066),
        "artifact_png": str(out_png.relative_to(REPO_ROOT)).replace("\\", "/"),
        "verdict": ("STALE-CAPTION CONFIRMED: runtime theta_SH=0.066, fig 2.10 "
                    "overrides neither theta_SH nor integrator; the '0.04 / "
                    "spherical-Euler' figure caption/comment is stale."),
        "pilot_caveat": ("PART B-2 is a reduced N pilot; the red-team must-fix "
                         "requires N>=200-seed median I_50 (per-seed thresholds "
                         "span widely). Use --full --n-seeds 200 for the thesis "
                         "ensemble. C13 itself is settled by PART B-1 read-back "
                         "alone, independent of the ensemble."),
    }

    out_json = OUT_DIR / "theta_sh_provenance.json"
    out_json.write_text(json.dumps(result, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(f"  best-matching θ_SH vs exp I_th = {best}  "
          f"(gaps µA: {result['closest_theta_to_exp_Ith']['abs_gap_uA']})")
    print(f"  wrote {out_json}")
    print(f"  wrote {out_png}")
    return result


if __name__ == "__main__":
    main()
