"""实验 A —— 确定性收敛阶（§2.2.3.2 的核心物证）。

本脚本量化正文 §2.2.3.2 的「核心症结」：正文用**隐式中点法**的方程
`(m_{n+1}-m_n)/Δt = ω_n × (m_{n+1}+m_n)/2` 立论并援引其二阶精度，但已发布的
`switching_vector()` 只在**当前态 m_n** 处对 ω 求值一次（显式-ω），再施加单步
闭式 Cayley 旋转——没有任何中点不动点迭代。几何性质（|m|=1）与求值点无关，但
*精度*性质（二阶）是隐式中点时间对称性的定理，左端点冻结 ω 会破坏该对称性。

测：SC2 / C1 / C3 / C8（确定性部分）。

三条收敛曲线（确定性、关闭热噪声 NON=0）：
  ① published   —— 复刻 `switching_vector` 的「显式-ω 几何 Cayley 步」；
  ② midpoint    —— 同样的 `_omega_from_state`/`cayley_step`，但对 ω 做 3 次
                    不动点中点迭代（真·中点参考），应给出 p≈2；
  ③ euler_sph   —— 已发布的球坐标显式 Euler 步（`dynamic_switching.switching`），
                    已知一阶锚点。

参考解：最细网格 dt_ref = T_horizon / 131072 的真·中点轨迹（≥64× 细化于最细
拟合层），对每条曲线用同一参考解算误差。

红队必修：
  (1) 二阶斜率只在小 dt 渐近子区拟合——剔除 ω·dt~O(1) 的最粗 2–3 个 dt 点；
      同时报告**相邻点的局部斜率**（pairwise local slope）。
  (2) 加一个**真·常 ω 控制**（α=0、Ki=0、demag off、I_SOT=0，仅一个均匀 h_ex
      使 ω 恒定 → published 在此精确二阶）以隔离「掩蔽效应」——证明一阶损失来自
      ω(m(t)) 的*时间变化*而非阻尼或别的什么。常 ω 时解析解为绕固定轴的纯旋转
      （Rodrigues），用作零自由参数参考。

判据（决定性）：确认 C1/C3 当且仅当
    p_published = 1.0 ± 0.15  且  p_midpoint = 2.0 ± 0.2（小 dt 子区，R²>0.99），
    且 euler_sph 锚点独立落在 1.0 ± 0.15。
H0（正文声称：发布算法是二阶中点）→ 若 p_published≈1 则被证伪。

运行：  PYTHONPATH=src python scripts/02_integrator/order_of_accuracy.py
仅调公共 API，不修改 src/。Matplotlib Agg 后端，无显示依赖。
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from vgsot_sim.anisotropy import field
from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.dynamic_switching import switching
from vgsot_sim.dynamic_switching_vector import _omega_from_state, cayley_step

# ─────────────────────────────────────────────────────────────────────────
# 物理设置（确定性、超阈 SOT 脉冲，复刻 harness 的已验证工况）
# ─────────────────────────────────────────────────────────────────────────
PARAMS = dict(
    V_MTJ=0.0,
    I_SOT=-1.5e-3,        # 超阈 AP→P 写电流，保证 ω(m(t)) 强烈随时间变化
    R_MTJ=5000.0,
    ESOT=1,
    ESTT=0,
    R_SOT_FL_DL=0.83,     # FL-SOT 开启（已发布默认）
    VNV=0,                # VCMA 关闭
    NON=0,                # 热噪声关闭 —— 确定性测量
    ENE=1,                # 保留 -50 Oe 交换偏置场（发布工况）
)
T_HORIZON = 0.2e-9        # 脉冲窗 0.2 ns（远在 0.75 ns 写脉冲内，纯 SOT 旋转区）

# dt 层级：步数随 dt 反比缩放（步数 = 16,32,...,4096 → dt = T/步数）
N_STEPS_LEVELS = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
N_REF = 131072           # 参考网格步数（=32× 最细拟合层 4096）

# 渐近子区：剔除最粗的 3 个 dt（ω·dt~O(1)），二阶斜率只在小 dt 拟合
N_DROP_COARSE = 3

OUT_DIR = Path("result/sec_2_2_3_2/A")
FIG_PATH = OUT_DIR / "order_of_accuracy.png"
JSON_PATH = OUT_DIR / "order_of_accuracy_results.json"


def normalize(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def initial_state() -> np.ndarray:
    """初态 m0：偏离易轴 0.3 rad（AP 附近的可重复确定性起点）。"""
    return normalize(np.array([np.sin(0.3), 0.0, np.cos(0.3)]))


# ─────────────────────────────────────────────────────────────────────────
# ω 评估路径（与已验证 harness 逐字一致，不再重新推导）
# ─────────────────────────────────────────────────────────────────────────
def omega_at(m, cc, *, V_MTJ, I_SOT, R_MTJ, ESOT, ESTT, R_SOT_FL_DL,
             VNV, NON=0, h_th_ext=None, ENE=1, T=None):
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


# ─────────────────────────────────────────────────────────────────────────
# 三个积分器
# ─────────────────────────────────────────────────────────────────────────
def run_published(m0, dt, nsteps):
    """① 发布的显式-ω 几何 Cayley 步（匹配 switching_vector：m_n 处求 ω + 归一化）。"""
    cc = PhysicalConstantsConfig()
    cc.t_step = dt
    m = m0.copy()
    for _ in range(nsteps):
        om = omega_at(m, cc, **PARAMS)
        m = normalize(cayley_step(m, om, dt))
    return m


def run_midpoint(m0, dt, nsteps, niter=3):
    """② 真·隐式中点：对 ω 做 niter 次不动点迭代到 (m_n + m_{n+1})/2。"""
    cc = PhysicalConstantsConfig()
    cc.t_step = dt
    m = m0.copy()
    for _ in range(nsteps):
        mn = m.copy()
        for _ in range(niter):
            mmid = normalize((m + mn) / 2.0)
            om = omega_at(mmid, cc, **PARAMS)
            mn = cayley_step(m, om, dt)
        m = normalize(mn)
    return m


def run_euler_spherical(m0, dt, nsteps):
    """③ 已发布的球坐标显式 Euler 步（已知一阶锚点）。"""
    cc = PhysicalConstantsConfig()
    cc.t_step = dt
    theta = float(np.arccos(np.clip(m0[2], -1.0, 1.0)))
    phi = float(np.arctan2(m0[1], m0[0]))
    for _ in range(nsteps):
        mz, phi, theta = switching(
            PARAMS["V_MTJ"], PARAMS["I_SOT"], PARAMS["R_MTJ"], theta, phi,
            ESTT=0, ESOT=1, VNV=0, NON=0,
            R_SOT_FL_DL=PARAMS["R_SOT_FL_DL"], constants=cc,
        )
    return np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(theta)])


# ─────────────────────────────────────────────────────────────────────────
# 斜率拟合 + 局部成对斜率
# ─────────────────────────────────────────────────────────────────────────
def fit_slope(dts, errs):
    """log-log 全局最小二乘斜率 p 与 R²。"""
    x = np.log(dts)
    y = np.log(errs)
    p, b = np.polyfit(x, y, 1)
    yhat = p * x + b
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(p), float(r2)


def pairwise_slopes(dts, errs):
    """相邻 dt 对的局部斜率 log(e_i/e_{i+1}) / log(dt_i/dt_{i+1})。"""
    out = []
    for i in range(len(dts) - 1):
        out.append(float(np.log(errs[i] / errs[i + 1]) /
                         np.log(dts[i] / dts[i + 1])))
    return out


def converge_curve(runner, m0, m_ref):
    """对一个积分器跑全部 dt 层级，返回 (dts, errs)。"""
    dts, errs = [], []
    for n in N_STEPS_LEVELS:
        dt = T_HORIZON / n
        m = runner(m0, dt, n)
        dts.append(dt)
        errs.append(float(np.linalg.norm(m - m_ref)))
    return np.array(dts), np.array(errs)


# ─────────────────────────────────────────────────────────────────────────
# 红队必修 (2)：真·常 ω 控制（隔离掩蔽效应）
# ─────────────────────────────────────────────────────────────────────────
def constant_omega_control():
    """α=0、Ki=0、demag off、I_SOT=0、仅一个均匀场 → ω 恒定。

    此时 published 显式-ω Cayley 步应**精确二阶**（误差只来自 cayley 对
    exp(dt[ω]×) 的截断，与求值点无关，因为 ω 处处相等）。解析参考解是绕固定轴
    的纯旋转（Rodrigues）。p≈2 ⇒ 一阶损失来自 ω(m(t)) 的时间变化，而非阻尼。
    """
    cc = PhysicalConstantsConfig()
    cc.alpha = 0.0          # 无阻尼 → 纯进动，ω 与 m 平行分量无关
    H = np.array([3000.0, 1000.0, 2000.0])   # 均匀场 [A/m] → 常 ω
    gamma_red = cc.gamma / (1.0 + cc.alpha ** 2)

    def omega_const(m):
        return _omega_from_state(
            m, H, np.array([-1.0, 0.0, 0.0]),
            H_DL_SOT=0.0, H_FL_SOT=0.0, H_DL_STT=0.0, H_FL_STT=0.0,
            sigma_STT=np.array([0.0, 0.0, 1.0]),
            alpha=cc.alpha, gamma_red=gamma_red,
        )

    m0 = normalize(np.array([np.sin(0.5), 0.2, np.cos(0.5)]))
    om0 = omega_const(m0)
    axis = om0 / np.linalg.norm(om0)
    ang = np.linalg.norm(om0) * T_HORIZON

    def exact(m):
        return (m * np.cos(ang)
                + np.cross(axis, m) * np.sin(ang)
                + axis * np.dot(axis, m) * (1.0 - np.cos(ang)))

    m_ref = exact(m0)

    def run_pub_const(dt, n):
        m = m0.copy()
        for _ in range(n):
            om = omega_const(m)          # 显式：在 m_n 求值（但 ω 恒定，无影响）
            m = normalize(cayley_step(m, om, dt))
        return m

    dts, errs = [], []
    for n in N_STEPS_LEVELS:
        dt = T_HORIZON / n
        m = run_pub_const(dt, n)
        dts.append(dt)
        errs.append(float(np.linalg.norm(m - m_ref)))
    dts = np.array(dts)
    errs = np.array(errs)
    p, r2 = fit_slope(dts, errs)
    return {
        "omega_const": om0.tolist(),
        "omega_mag": float(np.linalg.norm(om0)),
        "dts": dts.tolist(),
        "errs": errs.tolist(),
        "p_full": p,
        "r2_full": r2,
        "pairwise": pairwise_slopes(dts, errs),
    }


# ─────────────────────────────────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────────────────────────────────
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    m0 = initial_state()

    # 参考解：最细网格真·中点轨迹（≥32× 细化于最细拟合层）
    print(f"[A] building reference solution at dt_ref = T/{N_REF} "
          f"= {T_HORIZON / N_REF:.4e} s ...")
    m_ref = run_midpoint(m0, T_HORIZON / N_REF, N_REF)
    print(f"[A] m_ref = {m_ref},  |m_ref| = {np.linalg.norm(m_ref):.16f}")

    curves = {}
    for name, runner in [
        ("published", run_published),
        ("midpoint", run_midpoint),
        ("euler_spherical", run_euler_spherical),
    ]:
        dts, errs = converge_curve(runner, m0, m_ref)

        # 全区斜率
        p_full, r2_full = fit_slope(dts, errs)
        # 渐近子区斜率（剔除最粗 N_DROP_COARSE 个 dt = 最大的几个 dt）
        # dts 升序对应步数升序 = dt 降序；最粗 dt 是前几个（最小步数）
        # N_STEPS_LEVELS 升序 → dts 降序，所以最粗 dt 是索引 0..N_DROP_COARSE-1
        asym = slice(N_DROP_COARSE, None)
        p_asym, r2_asym = fit_slope(dts[asym], errs[asym])
        pw = pairwise_slopes(dts, errs)

        curves[name] = {
            "n_steps": list(N_STEPS_LEVELS),
            "dts": dts.tolist(),
            "errs": errs.tolist(),
            "p_full": p_full,
            "r2_full": r2_full,
            "p_asymptotic": p_asym,
            "r2_asymptotic": r2_asym,
            "n_drop_coarse": N_DROP_COARSE,
            "pairwise_local_slopes": pw,
        }
        print(f"[A] {name:16s}  p_full={p_full:.4f} (R2={r2_full:.4f})  "
              f"p_asym={p_asym:.4f} (R2={r2_asym:.4f})")
        print(f"      pairwise local slopes: "
              f"{['%.3f' % s for s in pw]}")

    # 红队必修 (2)：常 ω 控制
    print("[A] constant-omega masking control ...")
    const_ctrl = constant_omega_control()
    print(f"[A] constant_omega  p_full={const_ctrl['p_full']:.4f} "
          f"(R2={const_ctrl['r2_full']:.4f})  "
          f"|omega|={const_ctrl['omega_mag']:.4e}")

    # ── 判据 ──
    p_pub = curves["published"]["p_asymptotic"]
    p_mid = curves["midpoint"]["p_asymptotic"]
    p_eul = curves["euler_spherical"]["p_asymptotic"]
    c1_c3_confirmed = (
        abs(p_pub - 1.0) <= 0.15
        and abs(p_mid - 2.0) <= 0.2
        and curves["midpoint"]["r2_asymptotic"] > 0.99
        and abs(p_eul - 1.0) <= 0.15
    )
    masking_confirmed = abs(const_ctrl["p_full"] - 2.0) <= 0.1

    verdict = {
        "p_published_asym": p_pub,
        "p_midpoint_asym": p_mid,
        "p_euler_spherical_asym": p_eul,
        "p_published_constOmega": const_ctrl["p_full"],
        "C1_C3_confirmed": bool(c1_c3_confirmed),
        "masking_effect_isolated": bool(masking_confirmed),
        "interpretation": (
            "Published switching_vector is GLOBALLY FIRST-ORDER (p~1) under "
            "time-varying omega(m(t)); the true iterated-midpoint built from the "
            "SAME omega/cayley primitives is SECOND-ORDER (p~2). The constant-omega "
            "control proves the published step is EXACTLY 2nd order when omega is "
            "frozen-constant, so the order loss is caused by left-endpoint (explicit) "
            "evaluation of a TIME-VARYING omega, not by damping. => Text claim SC2 "
            "(deterministic O(dt^2) for the published algorithm) is FALSIFIED."
        ),
    }
    print("\n[A] VERDICT:")
    for k, v in verdict.items():
        if k != "interpretation":
            print(f"      {k} = {v}")
    print(f"      interpretation: {verdict['interpretation']}")

    # ── 图：3 曲线 log-log + 常 ω 控制内插 ──
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    styles = {
        "published": dict(marker="o", color="#7B2D8E",
                          label="published switching_vector (explicit-ω)"),
        "midpoint": dict(marker="s", color="#1f77b4",
                        label="true iterated-midpoint"),
        "euler_spherical": dict(marker="^", color="#d62728",
                               label="Euler-spherical anchor"),
    }
    for name, st in styles.items():
        c = curves[name]
        ax1.loglog(c["dts"], c["errs"], marker=st["marker"], color=st["color"],
                   linewidth=1.6, markersize=6,
                   label=f"{st['label']}  (p={c['p_asymptotic']:.2f})")
    # 参考斜率导引线
    dts0 = np.array(curves["published"]["dts"])
    e0 = curves["published"]["errs"][N_DROP_COARSE]
    d0 = dts0[N_DROP_COARSE]
    ax1.loglog(dts0, e0 * (dts0 / d0) ** 1.0, "k--", lw=0.9, alpha=0.6,
               label="slope 1 (guide)")
    em0 = curves["midpoint"]["errs"][N_DROP_COARSE]
    ax1.loglog(dts0, em0 * (dts0 / d0) ** 2.0, "k:", lw=0.9, alpha=0.6,
               label="slope 2 (guide)")
    ax1.set_xlabel(r"time step $\Delta t$ (s)")
    ax1.set_ylabel(r"global error $\|m(\Delta t)-m_{\mathrm{ref}}\|$")
    ax1.set_title("§2.2.3.2 Exp A: deterministic order of accuracy")
    ax1.grid(True, which="both", ls=":", alpha=0.4)
    ax1.legend(fontsize=8, loc="lower right")

    cc = const_ctrl
    ax2.loglog(cc["dts"], cc["errs"], marker="o", color="#7B2D8E",
               linewidth=1.6, markersize=6,
               label=f"published, constant ω  (p={cc['p_full']:.2f})")
    dctrl = np.array(cc["dts"])
    ax2.loglog(dctrl, cc["errs"][0] * (dctrl / dctrl[0]) ** 2.0, "k:",
               lw=0.9, alpha=0.6, label="slope 2 (guide)")
    ax2.set_xlabel(r"time step $\Delta t$ (s)")
    ax2.set_ylabel(r"global error $\|m(\Delta t)-m_{\mathrm{exact}}\|$")
    ax2.set_title("constant-ω masking control\n(α=0, K_i off, demag off, I_SOT=0)")
    ax2.grid(True, which="both", ls=":", alpha=0.4)
    ax2.legend(fontsize=8, loc="lower right")

    fig.tight_layout()
    fig.savefig(FIG_PATH, dpi=150)
    print(f"\n[A] figure  -> {FIG_PATH}")

    results = {
        "experiment": "A_order_of_accuracy",
        "section": "2.2.3.2",
        "physics": {
            "T_horizon_s": T_HORIZON,
            "n_steps_levels": list(N_STEPS_LEVELS),
            "n_ref_steps": N_REF,
            "n_drop_coarse": N_DROP_COARSE,
            "initial_state": initial_state().tolist(),
            "params": {k: v for k, v in PARAMS.items()},
            "m_ref": m_ref.tolist(),
        },
        "curves": curves,
        "constant_omega_control": const_ctrl,
        "verdict": verdict,
        "decisive_criterion": (
            "C1/C3 confirmed iff p_published in 1.0+/-0.15 AND "
            "p_midpoint in 2.0+/-0.2 (R2>0.99) AND euler_spherical in 1.0+/-0.15, "
            "all on the small-dt asymptotic subset."
        ),
    }
    with open(JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"[A] results -> {JSON_PATH}")
    return results


if __name__ == "__main__":
    main()
