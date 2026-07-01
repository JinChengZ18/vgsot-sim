"""
实验 D-2 (C11/SC7 动机) —— 球坐标 1/sinθ 极点奇点的隔离与正文措辞纠偏。

§2.2.3.2 / §2.2.5 称 Cayley 向量积分器「无 1/sinθ 极点」。本脚本用**单步 RHS 探针**
（无时间积分）把这句话拆成两半验证：

  ① 球坐标 Euler 步进 `dynamic_switching.switching()` 的 dφ/dt 在 θ→0 处确实
     按 1/sinθ 发散（log-log 斜率 -1.00）；|dφ/dt|·sinθ 仅在 θ≲1e-3 才趋于常数
     （θ=0.1 时 cosθ≈1 的 O(θ²) 修正尚未可忽略）。
  ② 仓库已发布的极点保护（strict `|sinθ|<1e-8` → dφ/dt=0）实际只在 **θ≤1e-9**
     才触发，而非 θ≤1e-8。我们用 GUARD-ON（直接调 published switching()）对照
     GUARD-OFF（去掉保护的逐式重写）展示触发边界。

红队修订（关键）：**不**对比 Cayley 路径的 dφ_cay/dt——dφ/dt 是坐标伪量，
在极点附近本就无定义。Cayley 步进内部 (dynamic_switching_vector.py:178-179) 仍用
arccos/arctan2 把 m 转成 (θ,φ) 去调 field()，但 **field() 内部没有 1/sinθ 除法**
（H_eff 全部按 m 的笛卡尔分量组装，见 anisotropy.py:35-62）。因此我们改为展示
**笛卡尔 |dm/dt|_cay 的正则性**（θ→0 时收敛到有限值），并明确指出：

  「无 1/sinθ 极点」只**部分**成立——正确表述应为
  「内核场组装 field() 规避了 1/sinθ 除法（笛卡尔 dm/dt 处处有限），
    而 switching_vector 仍执行 arccos/arctan2 球坐标往返（仅用于喂 field 的角参数，
    不进入 1/sinθ 分母）」。

按 SPEC 加入一般面内场 h_ex_x=4e4 A/m（打破 φ 对称、让 dφ/dt 在极点附近非零，
否则若 H_eff 恰好与子午面对齐 dφ/dt 可能偶然为 0 而掩盖奇点）。

仅依赖公共 API，不修改 src/。matplotlib Agg 后端。

运行：  PYTHONPATH=src python scripts/11_pole_singularity/probe_pole_singularity.py
产物：  result/sec_2_2_3_2/D/probe_pole_singularity.png
        result/sec_2_2_3_2/D/probe_pole_singularity.json
"""
from __future__ import annotations

import dataclasses
import json
from math import sin, cos
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.anisotropy import field
from vgsot_sim.dynamic_switching import switching
from vgsot_sim.dynamic_switching_vector import _omega_from_state

# ── Style ────────────────────────────────────────────────────────────────
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

# Probe excitation (SOT pulse on, near the easy axis)
PHI0 = 0.3
V_MTJ = 0.0
I_SOT = -1.5e-3
R_MTJ = 5000.0
ESOT, ESTT = 1, 0
R_SOT_FL_DL = 0.83
GUARD_EPS = 1e-8           # strict |sin(theta)| < 1e-8 in switching()


def _torque_prefactors(cc):
    Ms_use = cc.Ms
    I_MTJ = V_MTJ / R_MTJ if R_MTJ else 0.0
    J_STT = I_MTJ / cc.A1
    J_SOT = I_SOT / cc.A2
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)
    H_DL_STT = ESTT * cc.h_bar * cc.P * J_STT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_FL_STT = 0.0
    H_DL_SOT = ESOT * cc.h_bar * cc.theta_SH * J_SOT / (2 * cc.e * cc.u0 * Ms_use * cc.tf)
    H_FL_SOT = R_SOT_FL_DL * H_DL_SOT
    return gamma_red, H_DL_STT, H_FL_STT, H_DL_SOT, H_FL_SOT


def dphi_dt_guard_off(theta, phi, cc):
    """Re-implementation of switching()'s dphi/dt WITHOUT the pole guard.

    Byte-identical to dynamic_switching.switching() lines 117-127 except the
    `if abs(sin_theta) < 1e-8` branch is removed, so 1/sin(theta) is always
    used. This exposes the bare 1/sin(theta) pole.
    """
    H_EFF, _ = field(theta, phi, V_MTJ, 1, 0, 1, 0, cc, demag_mode="ellipsoid")
    gamma_red, H_DL_STT, H_FL_STT, H_DL_SOT, H_FL_SOT = _torque_prefactors(cc)
    inv_sin = 1.0 / sin(theta)                  # GUARD OFF: no clamp
    dphi_dt = gamma_red * (
        inv_sin * (
            H_EFF[0] * (-cc.alpha * sin(phi) - cos(theta) * cos(phi))
            + H_EFF[1] * (cc.alpha * cos(phi) - cos(theta) * sin(phi))
            + H_EFF[2] * sin(theta))
        - (cc.alpha * H_DL_STT + H_FL_STT)
        - H_DL_SOT * inv_sin * (cc.alpha * cos(theta) * cos(phi) - sin(phi))
        - H_FL_SOT * inv_sin * (cc.alpha * sin(phi) + cos(theta) * cos(phi))
    )
    return dphi_dt


def dphi_dt_guard_on(theta, phi, cc):
    """dphi/dt as actually produced by the published switching() (guard ON).

    Read back from (phi_new - phi0)/t_step of the real stepper, so the guard
    (abs(sin_theta) < 1e-8 -> dphi_dt = 0) is exercised exactly as shipped.
    """
    _mz, phi_new, _theta_new = switching(
        V_MTJ, I_SOT, R_MTJ, theta, phi, ESTT, ESOT,
        VNV=0, R_SOT_FL_DL=R_SOT_FL_DL, constants=cc)
    return (phi_new - phi) / cc.t_step


def dmdt_cartesian_cayley(theta, phi, cc):
    """|dm/dt| in Cartesian form via the published _omega_from_state path.

    This is the RED-TEAM-correct quantity: the Cayley/vector stepper's RHS is
    dm/dt = omega x m, assembled from field() (no 1/sin division). It stays
    finite as theta -> 0, unlike the coordinate-artifact dphi/dt.
    """
    m = np.array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])
    H_eff, _ = field(theta, phi, V_MTJ, 1, 0, 1, 0, cc, demag_mode="ellipsoid")
    H_eff = np.asarray(H_eff, float)
    gamma_red, H_DL_STT, H_FL_STT, H_DL_SOT, H_FL_SOT = _torque_prefactors(cc)
    omega = _omega_from_state(
        m, H_eff, np.array([-1.0, 0.0, 0.0]),
        H_DL_SOT=H_DL_SOT, H_FL_SOT=H_FL_SOT,
        H_DL_STT=H_DL_STT, H_FL_STT=H_FL_STT,
        sigma_STT=np.array([0.0, 0.0, 1.0]),
        alpha=cc.alpha, gamma_red=gamma_red)
    dmdt = np.cross(omega, m)
    return float(np.linalg.norm(dmdt))


def main():
    # SPEC: add in-plane h_ex_x = 4e4 A/m to break phi symmetry near the pole.
    cc = dataclasses.replace(PhysicalConstantsConfig(), h_ex_x=4.0e4)

    thetas = np.array([10.0 ** (-k) for k in range(1, 11)])   # 1e-1 .. 1e-10

    rows = []
    for th in thetas:
        d_off = dphi_dt_guard_off(th, PHI0, cc)
        d_on = dphi_dt_guard_on(th, PHI0, cc)
        guard_fires = bool(abs(sin(th)) < GUARD_EPS)
        dmdt_cart = dmdt_cartesian_cayley(th, PHI0, cc)
        rows.append({
            "theta": float(th),
            "sin_theta": float(sin(th)),
            "dphi_dt_guard_off": float(d_off),
            "dphi_dt_guard_off_times_sin": float(abs(d_off) * sin(th)),
            "dphi_dt_guard_on": float(d_on),
            "guard_fires": guard_fires,
            "dmdt_cartesian_cayley": dmdt_cart,
        })

    # ── slope of log|dphi/dt| vs log theta (1/sin => slope -1) ───────────
    th_arr = np.array([r["theta"] for r in rows])
    dphi_off = np.array([abs(r["dphi_dt_guard_off"]) for r in rows])
    mask_small = th_arr <= 1e-3
    slope = float(np.polyfit(np.log10(th_arr[mask_small]),
                             np.log10(dphi_off[mask_small]), 1)[0])

    # constancy of |dphi/dt|*sin(theta). The product reaches its 1/sin plateau
    # only for theta <= 1e-3 (at theta=0.1 the cos(theta)~1 O(theta^2) correction
    # is still ~96% off). Report CV over both the SPEC window (theta<=1e-3, which
    # still includes the partly-settled boundary point) and the deep-pole window
    # (theta<=1e-4, where it is flat to 6 digits) — honest about the boundary.
    prod_all = np.array([r["dphi_dt_guard_off_times_sin"] for r in rows])
    prod_le_1e3 = prod_all[mask_small]
    cv_le_1e3 = float(np.std(prod_le_1e3) / np.mean(prod_le_1e3))
    mask_deep = th_arr <= 1e-4
    prod_deep = prod_all[mask_deep]
    cv_deep = float(np.std(prod_deep) / np.mean(prod_deep))
    plateau = float(np.mean(prod_deep))           # deep-pole plateau value
    prod_at_1em1 = rows[0]["dphi_dt_guard_off_times_sin"]
    reldev_1em1 = abs(prod_at_1em1 - plateau) / plateau

    # guard firing boundary
    first_fire_theta = None
    for r in rows:
        if r["guard_fires"]:
            first_fire_theta = r["theta"]
            break

    # Cartesian regularity: ratio of max to min over the whole sweep (≈1 => flat)
    dmdt_cart = np.array([r["dmdt_cartesian_cayley"] for r in rows])
    cart_span = float(dmdt_cart.max() / dmdt_cart.min())

    # ── Figure: 3 panels ────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4))

    # (a) |dphi/dt| guard-off vs guard-on, with 1/sin reference
    ax = axes[0]
    ax.loglog(th_arr, dphi_off, "o-", color=CRIMSON, lw=1.8, ms=6,
              label=r"GUARD-OFF $|d\phi/dt|$")
    dphi_on = np.array([abs(r["dphi_dt_guard_on"]) for r in rows])
    # guard-on collapses to 0 where it fires; clamp for log display
    dphi_on_disp = np.clip(dphi_on, 1e-2, None)
    ax.loglog(th_arr, dphi_on_disp, "s--", color=NAVY, lw=1.4, ms=5,
              label=r"GUARD-ON (published) $|d\phi/dt|$")
    ref = dphi_off[mask_small][0] * (th_arr / th_arr[mask_small][0]) ** (-1.0)
    ax.loglog(th_arr[mask_small], dphi_off[mask_small][0]
              * (th_arr[mask_small] / th_arr[mask_small][0]) ** (-1.0),
              ":", color=CHARCOAL, lw=1.2, label=r"$\propto 1/\sin\theta$ 参考 (slope $-1$)")
    for r in rows:
        if r["guard_fires"]:
            ax.axvline(r["theta"], color=AMBER, ls="-.", lw=1.0, alpha=0.7)
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel(r"$|d\phi/dt|$ (rad/s)")
    ax.set_title(rf"(a) 球坐标 $d\phi/dt$ 极点 (slope ${slope:.2f}$)")
    ax.invert_xaxis()
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="upper right", framealpha=0.9)

    # (b) |dphi/dt|*sin(theta): constant only for theta<=1e-3
    ax = axes[1]
    prod = np.array([r["dphi_dt_guard_off_times_sin"] for r in rows])
    ax.semilogx(th_arr, prod, "o-", color=TEAL, lw=1.8, ms=6)
    ax.axhline(plateau, color=CHARCOAL, ls="--", lw=1.1,
               label=rf"plateau $={plateau:.3e}$")
    ax.axvline(1e-3, color=AMBER, ls="-.", lw=1.1,
               label=r"$\theta=10^{-3}$ 恒定起点")
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel(r"$|d\phi/dt|\cdot\sin\theta$ (rad/s)")
    ax.set_title(r"(b) $|d\phi/dt|\cdot\sin\theta$ 仅 $\theta \leq 10^{-3}$ 恒定")
    ax.invert_xaxis()
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="upper right", framealpha=0.9)

    # (c) RED-TEAM FIX: Cartesian |dm/dt| (Cayley path) stays finite
    ax = axes[2]
    ax.loglog(th_arr, dmdt_cart, "D-", color=THU_DEEP, lw=1.8, ms=6,
              label=r"笛卡尔 $|dm/dt|_{\rm cay}$ (有限)")
    ax.loglog(th_arr, dphi_off, "o:", color=CRIMSON, lw=1.0, ms=4, alpha=0.6,
              label=r"对照 $|d\phi/dt|$ (坐标伪量, 发散)")
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel(r"$|\,\cdot\,|$ (rad/s)")
    ax.set_title("(c) 笛卡尔 RHS 正则 — field() 无 $1/\\sin\\theta$ 除法")
    ax.invert_xaxis()
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="center right", framealpha=0.9)

    fig.suptitle("球坐标 $1/\\sin\\theta$ 极点的隔离与「无极点」措辞纠偏（单步 RHS 探针）",
                 fontsize=13, color=THU_DEEP, y=1.00)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    png_path = OUT_DIR / "probe_pole_singularity.png"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ── Decisive-criterion checks ───────────────────────────────────────
    checks = {
        "guard_off_slope": slope,
        "guard_off_slope_within_-1.00pm0.05": bool(abs(slope - (-1.0)) <= 0.05),
        "product_cv_for_theta_le_1e-3": cv_le_1e3,
        "product_cv_for_theta_le_1e-4": cv_deep,
        "product_plateau_value": plateau,
        "product_constant_from_theta_le_1e-3_scale": bool(cv_le_1e3 < 5e-3),
        "product_constant_deep_pole_theta_le_1e-4": bool(cv_deep < 1e-3),
        "product_reldev_at_theta_0.1": reldev_1em1,
        "product_NOT_constant_at_theta_0.1": bool(reldev_1em1 > 0.05),
        "guard_first_fires_at_theta": first_fire_theta,
        "guard_fires_at_theta_le_1e-9_not_1e-8": bool(
            first_fire_theta is not None and first_fire_theta <= 1e-9 + 1e-12),
        "cartesian_dmdt_span_max_over_min": cart_span,
        "cartesian_dmdt_finite_regular": bool(cart_span < 2.0),
    }

    payload = {
        "experiment": "D2_pole_singularity_probe",
        "section": "2.2.3.2",
        "claims_tested": ["SC7", "C11"],
        "description": ("Single-step RHS probe isolating the spherical 1/sin(theta) "
                        "pole (guard-on switching() vs guard-off reimplementation) "
                        "and the Cartesian dm/dt regularity that motivates the "
                        "'no-pole' wording correction."),
        "probe": {"phi0": PHI0, "V_MTJ": V_MTJ, "I_SOT": I_SOT, "R_MTJ": R_MTJ,
                  "ESOT": ESOT, "ESTT": ESTT, "R_SOT_FL_DL": R_SOT_FL_DL,
                  "h_ex_x_added": 4.0e4, "guard_eps": GUARD_EPS},
        "redteam_fixes_applied": [
            "do NOT compare dphi_cay/dt (coordinate artifact); show Cartesian |dm/dt|_cay regularity",
            "state field() has no 1/sin division while switching_vector still does arccos/arctan2 round-trip",
            "'no pole' wording is only PARTLY true -> 'kernel field assembly avoids 1/sin division'",
            "guard actually fires at theta<=1e-9 (strict |sin|<1e-8), not 1e-8",
        ],
        "wording_correction": (
            "正文「无 1/sinθ 极点」应改为「内核场组装 field() 规避了 1/sinθ 除法"
            "（笛卡尔 dm/dt 处处有限），而 switching_vector 仍执行 arccos/arctan2 "
            "球坐标往返，仅作为喂 field 的角参数、不进入 1/sinθ 分母」。"),
        "rows": rows,
        "decisive_checks": checks,
        "config": {"t_step_default": cc.t_step, "theta_SH": cc.theta_SH},
    }
    json_path = OUT_DIR / "probe_pole_singularity.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # ── Console report ──────────────────────────────────────────────────
    print("=== D2: spherical 1/sin(theta) pole probe ===")
    print(f"  GUARD-OFF log-log slope (theta<=1e-3) = {slope:.4f}  "
          f"(within -1.00+-0.05: {checks['guard_off_slope_within_-1.00pm0.05']})")
    print(f"  |dphi/dt|*sin(theta) CV (theta<=1e-3) = {cv_le_1e3:.3e}  "
          f"(const from 1e-3 scale: {checks['product_constant_from_theta_le_1e-3_scale']})")
    print(f"  |dphi/dt|*sin(theta) CV (theta<=1e-4) = {cv_deep:.3e}  "
          f"(deep-pole constant: {checks['product_constant_deep_pole_theta_le_1e-4']})")
    print(f"  |dphi/dt|*sin(theta) rel-dev at theta=0.1 = {reldev_1em1:.3e}  "
          f"(NOT constant: {checks['product_NOT_constant_at_theta_0.1']})")
    print(f"  guard first fires at theta = {first_fire_theta:.0e}  "
          f"(<=1e-9 not 1e-8: {checks['guard_fires_at_theta_le_1e-9_not_1e-8']})")
    print(f"  Cartesian |dm/dt| max/min over sweep = {cart_span:.4f}  "
          f"(finite/regular: {checks['cartesian_dmdt_finite_regular']})")
    print(f"PNG : {png_path}")
    print(f"JSON: {json_path}")


if __name__ == "__main__":
    main()
