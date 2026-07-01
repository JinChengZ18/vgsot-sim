"""
实验 D-1 (C11/C12 动机) —— 显式 Euler 切向更新的模长爆炸 vs Cayley 保模长。

§2.2.3.2 的强声称 SC6/SC7：显式 Euler `m_{n+1}=m_n+Δt(ω×m_n)` 不保模长，
而 Cayley 变换 `cayley_step` 在任意 Δt 下精确保 |m|=1。本脚本用**单步 RHS 探针**
（无时间积分）量化二者，并按 |ω| 量级分两块独立设阈以纠正正文的两处夸大措辞。

三条曲线：
  - euler_tangential(m,w,dt) = m + dt*cross(w,m)   —— 显式切向 Euler（不保模）
  - cayley_step(m,w,dt)                            —— 仓库已发布的闭式 Cayley（精确保模）
  - Rodrigues 旋转参考 R(ω̂, |ω_⊥|dt)·m            —— 真实测地旋转的标尺（精确保模）

两块独立设阈的 |ω| 面板：
  (a) 人工 |ω|=1e12（m0 ⟂ ω，预期 Δt=1e-12→|m|-1=+0.414, Δt=1e-11→+9.05）；
  (b) 真实 |ω|≈3e10（由 field()+SOT 在近 AP 态实测得到），发散点在 Δt≳3e-11
      （即 Δt·|ω|≳1），**而非** 1e-11。

红队修订：
  - 增长律按 |ω_⊥|=|ω×m| 表述（解析律 √(1+(Δt|ω_⊥|)²)），不按裸 |ω|；
    当 m0 ⟂ ω 时 |ω_⊥|=|ω|，二者重合（人工面板）；真实态 m 与 ω 非正交，
    |ω_⊥|<|ω|，故按 |ω_⊥| 表述才与解析律逐位吻合。
  - 阈值按 |ω| 分面板，避免把人工 1e12 的「1e-12 发散」误植到真实物理步长。

仅依赖公共 API（dynamic_switching_vector.cayley_step / _omega_from_state、
anisotropy.field），不修改 src/。matplotlib Agg 后端。

运行：  PYTHONPATH=src python scripts/09_simulation_figures/c12_norm_stability_sweep.py
产物：  result/sec_2_2_3_2/D/c12_norm_stability_sweep.png
        result/sec_2_2_3_2/D/c12_norm_stability_sweep.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.anisotropy import field
from vgsot_sim.dynamic_switching_vector import _omega_from_state, cayley_step

# ── Style (matches the 09_simulation_figures palette) ────────────────────
THU_DEEP = "#660874"
CRIMSON, NAVY, TEAL, AMBER = "#A82038", "#1F5FA8", "#1A6B5A", "#C47A00"
CHARCOAL = "#2B2B2B"
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12.5,
    "legend.fontsize": 9.0,
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False,
    "axes.linewidth": 0.9,
})

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "result" / "sec_2_2_3_2" / "D"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ── Integrators / references (single-step RHS probes, no integration) ────
def euler_tangential(m, w, dt):
    """Explicit tangential-Euler update m + dt (ω × m). NOT norm-preserving."""
    return np.asarray(m, float) + dt * np.cross(np.asarray(w, float), np.asarray(m, float))


def rodrigues_rotate(m, w, dt):
    """Exact geodesic rotation of m about ω̂ by angle |ω_⊥|·dt (norm-preserving).

    The physical precession rate of m is set by ω_⊥ = ω × m (the component of ω
    that actually rotates m); we use |ω_⊥|·dt as the rotation angle so the
    reference matches the analytic growth law's argument exactly.
    """
    m = np.asarray(m, float)
    w = np.asarray(w, float)
    wn = np.linalg.norm(w)
    if wn == 0.0:
        return m.copy()
    k = w / wn
    w_perp = np.cross(w, m)
    angle = np.linalg.norm(w_perp) * dt
    # Rodrigues: m cosθ + (k×m) sinθ + k (k·m)(1-cosθ)
    return (m * np.cos(angle)
            + np.cross(k, m) * np.sin(angle)
            + k * np.dot(k, m) * (1.0 - np.cos(angle)))


def omega_at(m, cc, *, V_MTJ, I_SOT, R_MTJ, ESOT, ESTT, R_SOT_FL_DL,
             VNV, NON=0, h_th_ext=None, ENE=1, T=None):
    """Verified omega path — reuses the published _omega_from_state."""
    m = np.asarray(m, float)
    theta = float(np.arccos(np.clip(m[2], -1.0, 1.0)))
    phi = float(np.arctan2(m[1], m[0]))
    H_eff, _ = field(theta, phi, V_MTJ, n=1, NON=NON, ENE=ENE, VNV=VNV, constants=cc,
                     demag_mode="ellipsoid", h_th_ext=h_th_ext, T=T)
    H_eff = np.asarray(H_eff, float)
    Ms_use = cc.Ms
    I_MTJ = V_MTJ / R_MTJ if R_MTJ else 0.0
    J_STT = I_MTJ / cc.A1
    J_SOT = I_SOT / cc.A2
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
    H_DL_STT = ESTT * cc.h_bar * cc.P * J_STT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_DL_SOT = ESOT * cc.h_bar * cc.theta_SH * J_SOT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_FL_SOT = R_SOT_FL_DL * H_DL_SOT
    return _omega_from_state(m, H_eff, np.array([-1.0, 0.0, 0.0]),
                             H_DL_SOT=H_DL_SOT, H_FL_SOT=H_FL_SOT,
                             H_DL_STT=H_DL_STT, H_FL_STT=0.0,
                             sigma_STT=np.array([0.0, 0.0, 1.0]),
                             alpha=cc.alpha, gamma_red=gamma_red)


def analytic_norm(dt, w_perp_mag):
    """Closed-form |m| after one tangential-Euler step: √(1+(Δt|ω_⊥|)²).

    Derivation: m1 = m + dt (ω×m); since (ω×m) ⟂ m and |ω×m|=|ω_⊥|,
    |m1|² = |m|² + (dt|ω_⊥|)² = 1 + (dt|ω_⊥|)².
    """
    return np.sqrt(1.0 + (dt * w_perp_mag) ** 2)


def sweep_panel(m0, w, dt_list):
    """Run all three integrators across dt; return per-dt |m|-1 dicts."""
    m0 = np.asarray(m0, float)
    w = np.asarray(w, float)
    w_perp_mag = float(np.linalg.norm(np.cross(w, m0)))
    rows = []
    for dt in dt_list:
        n_eu = float(np.linalg.norm(euler_tangential(m0, w, dt)))
        n_ca = float(np.linalg.norm(cayley_step(m0, w, dt)))
        n_ro = float(np.linalg.norm(rodrigues_rotate(m0, w, dt)))
        n_an = float(analytic_norm(dt, w_perp_mag))
        rows.append({
            "dt": float(dt),
            "dt_times_omega": float(dt * np.linalg.norm(w)),
            "dt_times_omega_perp": float(dt * w_perp_mag),
            "euler_norm_minus1": n_eu - 1.0,
            "cayley_norm_minus1": n_ca - 1.0,
            "rodrigues_norm_minus1": n_ro - 1.0,
            "analytic_norm_minus1": n_an - 1.0,
        })
    return {"omega_mag": float(np.linalg.norm(w)),
            "omega_perp_mag": w_perp_mag,
            "rows": rows}


def main():
    cc = PhysicalConstantsConfig()

    # ── Panel (a): ARTIFICIAL |ω|=1e12, m0 ⟂ ω  (|ω_⊥| = |ω|) ───────────
    w_art = np.array([0.0, 0.0, 1.0e12])
    m_art = np.array([1.0, 0.0, 0.0])           # perpendicular to ω
    dt_art = np.array([1e-14, 1e-13, 1e-12, 3e-12, 1e-11, 3e-11])
    res_art = sweep_panel(m_art, w_art, dt_art)

    # ── Panel (b): REALISTIC |ω|≈3e10 from field()+SOT (near-AP state) ──
    m_real = np.array([0.0, 0.1, np.sqrt(1.0 - 0.1 ** 2)])
    m_real = m_real / np.linalg.norm(m_real)
    w_real = omega_at(m_real, cc, V_MTJ=0.0, I_SOT=-1.5e-3, R_MTJ=5000.0,
                      ESOT=1, ESTT=0, R_SOT_FL_DL=0.83, VNV=0, NON=0)
    dt_real = np.array([1e-13, 1e-12, 1e-11, 3e-11, 5e-11, 1e-10, 3e-10])
    res_real = sweep_panel(m_real, w_real, dt_real)

    # ── Pull out the decisive criterion numbers for the artificial panel ─
    def get(rows, dt):
        # rtol-only match (default atol=1e-8 is far too loose for dt~1e-12)
        for r in rows:
            if np.isclose(r["dt"], dt, rtol=1e-6, atol=0.0):
                return r
        return None
    art_1em12 = get(res_art["rows"], 1e-12)
    art_1em11 = get(res_art["rows"], 1e-11)

    # realistic divergence onset: first dt with dt*|ω| >= 1
    real_onset = None
    for r in res_real["rows"]:
        if r["dt_times_omega"] >= 1.0:
            real_onset = r["dt"]
            break

    # ── Figure: two independently-thresholded |ω| panels ────────────────
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))

    for ax, res, title, dt_arr in (
        (axes[0], res_art, r"(a) 人工 $|\omega|=10^{12}$ rad/s ($\mathbf{m}_0\perp\omega$)", dt_art),
        (axes[1], res_real, rf"(b) 真实 $|\omega|\approx{res_real['omega_mag']:.2e}$ rad/s (field+SOT)", dt_real),
    ):
        rows = res["rows"]
        dts = np.array([r["dt"] for r in rows])
        eu = np.array([abs(r["euler_norm_minus1"]) for r in rows])
        ca = np.array([abs(r["cayley_norm_minus1"]) for r in rows])
        ro = np.array([abs(r["rodrigues_norm_minus1"]) for r in rows])
        an = np.array([abs(r["analytic_norm_minus1"]) for r in rows])
        # clamp Cayley/Rodrigues machine-zero for log plot
        floor = 1e-17
        ca = np.clip(ca, floor, None)
        ro = np.clip(ro, floor, None)
        ax.loglog(dts, np.clip(eu, floor, None), "o-", color=CRIMSON, lw=1.8, ms=6,
                  label="显式 Euler 切向  $|m|-1$")
        ax.loglog(dts, np.clip(an, floor, None), "--", color=CHARCOAL, lw=1.3,
                  label=r"解析律 $\sqrt{1+(\Delta t|\omega_\perp|)^2}-1$")
        ax.loglog(dts, ca, "s-", color=NAVY, lw=1.4, ms=5,
                  label="Cayley  $|m|-1$ (保模长)")
        ax.loglog(dts, ro, "^:", color=TEAL, lw=1.2, ms=5,
                  label="Rodrigues 参考  $|m|-1$")
        # mark dt where dt*|omega| = 1
        dt_cross = 1.0 / res["omega_mag"]
        ax.axvline(dt_cross, color=AMBER, ls="-.", lw=1.1,
                   label=rf"$\Delta t\,|\omega|=1$  ($\Delta t={dt_cross:.2e}$)")
        ax.set_xlabel(r"时间步长 $\Delta t$ (s)")
        ax.set_ylabel(r"$\left|\,|m|-1\,\right|$")
        ax.set_title(title)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(loc="upper left", framealpha=0.9)

    # annotate the artificial decisive values
    axes[0].annotate(rf"$\Delta t=10^{{-12}}$: $|m|-1={art_1em12['euler_norm_minus1']:+.3f}$"
                     "\n"
                     rf"$\Delta t=10^{{-11}}$: $|m|-1={art_1em11['euler_norm_minus1']:+.3f}$",
                     xy=(0.97, 0.04), xycoords="axes fraction",
                     ha="right", va="bottom", fontsize=8.5,
                     bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=THU_DEEP, alpha=0.9))

    fig.suptitle("显式 Euler 切向更新模长爆炸 vs Cayley 保模长（单步 RHS 探针，无时间积分）",
                 fontsize=13, color=THU_DEEP, y=1.00)
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    png_path = OUT_DIR / "c12_norm_stability_sweep.png"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ── Decisive-criterion checks ───────────────────────────────────────
    checks = {
        "artificial_dt1e-12_norm_minus1": art_1em12["euler_norm_minus1"],
        "artificial_dt1e-12_expect_+0.414": bool(abs(art_1em12["euler_norm_minus1"] - 0.414214) < 1e-4),
        "artificial_dt1e-11_norm_minus1": art_1em11["euler_norm_minus1"],
        "artificial_dt1e-11_expect_+9.05": bool(abs(art_1em11["euler_norm_minus1"] - 9.049876) < 1e-3),
        "artificial_euler_matches_analytic": bool(all(
            np.isclose(r["euler_norm_minus1"], r["analytic_norm_minus1"], rtol=1e-9, atol=1e-12)
            for r in res_art["rows"])),
        "artificial_cayley_max_abs_norm_dev": float(max(
            abs(r["cayley_norm_minus1"]) for r in res_art["rows"])),
        "artificial_cayley_below_1e-14": bool(max(
            abs(r["cayley_norm_minus1"]) for r in res_art["rows"]) < 1e-14),
        "realistic_omega_mag": res_real["omega_mag"],
        "realistic_dt_cross_omega1": float(1.0 / res_real["omega_mag"]),
        "realistic_divergence_onset_dt": real_onset,
        "realistic_onset_is_3e-11_not_1e-11": bool(
            real_onset is not None and real_onset >= 3e-11 - 1e-13),
        "realistic_cayley_max_abs_norm_dev": float(max(
            abs(r["cayley_norm_minus1"]) for r in res_real["rows"])),
    }

    payload = {
        "experiment": "D1_norm_stability_sweep",
        "section": "2.2.3.2",
        "claims_tested": ["SC6", "SC7", "C11", "C12"],
        "description": ("Single-step RHS probe of explicit tangential-Euler norm "
                        "blow-up vs Cayley/Rodrigues norm preservation, two "
                        "independently-thresholded |omega| panels."),
        "growth_law": "|m| = sqrt(1 + (dt*|omega_perp|)^2), omega_perp = omega x m",
        "redteam_fixes_applied": [
            "growth expressed via |omega_perp|=|omega x m| not bare |omega|",
            "thresholds split per-panel by |omega| (artificial 1e12 vs realistic ~3e10)",
            "realistic divergence onset reported at dt*|omega|>=1 (~3e-11), not 1e-11",
        ],
        "artificial_panel": res_art,
        "realistic_panel": res_real,
        "decisive_checks": checks,
        "config": {
            "t_step_default": cc.t_step, "theta_SH": cc.theta_SH,
            "I_SOT_realistic": -1.5e-3, "R_MTJ": 5000.0,
        },
    }
    json_path = OUT_DIR / "c12_norm_stability_sweep.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # ── Console report (real numbers) ───────────────────────────────────
    print("=== D1: Euler norm blow-up vs Cayley/Rodrigues ===")
    print(f"[artificial |omega|=1e12, m0 perp omega]")
    print(f"  dt=1e-12: euler |m|-1 = {art_1em12['euler_norm_minus1']:+.6f}  "
          f"(expect +0.414214, match={checks['artificial_dt1e-12_expect_+0.414']})")
    print(f"  dt=1e-11: euler |m|-1 = {art_1em11['euler_norm_minus1']:+.6f}  "
          f"(expect +9.049876, match={checks['artificial_dt1e-11_expect_+9.05']})")
    print(f"  euler==analytic (all dt): {checks['artificial_euler_matches_analytic']}")
    print(f"  cayley max|m|-1 = {checks['artificial_cayley_max_abs_norm_dev']:.3e}  "
          f"(<1e-14: {checks['artificial_cayley_below_1e-14']})")
    print(f"[realistic |omega|={res_real['omega_mag']:.4e}]")
    print(f"  dt*|omega|=1 at dt = {checks['realistic_dt_cross_omega1']:.4e}")
    print(f"  divergence onset (first dt with dt*|omega|>=1) = {real_onset:.3e}  "
          f"(>=3e-11 not 1e-11: {checks['realistic_onset_is_3e-11_not_1e-11']})")
    print(f"  cayley max|m|-1 = {checks['realistic_cayley_max_abs_norm_dev']:.3e}")
    print(f"PNG : {png_path}")
    print(f"JSON: {json_path}")


if __name__ == "__main__":
    main()
