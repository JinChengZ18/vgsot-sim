"""
参数标定 — 让 vgsot-sim 的 SER 阈值匹配 Device A P→AP @ 5 ns 实验数据。

实验靶值（来自章节 §2.3.4，τ_0 = 1 ns 先验，Device A P→AP 同批次）：
    a = 0.793 V, b = 0.175 V (V_th = a − b·ln(t_w/ns))
    V_th(0.75 ns)  ≈ 844 mV  →  I_SOT ≈ 1.09 mA
    V_th(5 ns)     ≈ 511 mV  →  I_SOT ≈ 659 µA

当前仿真：SER 50% 阈值 ~140 µA（V_SOT ≈ 109 mV），约比实验低 4.7×。

候选调参（按效应排序）：
    1. θ_SH:  减小 → 提高 I_c0
    2. K_i:   增大 → 同时增大 Δ 与 I_c0
    3. M_s:   增大 → 通过 H_k 影响 I_c0
    4. α:     增大 → 线性增大 I_c0（CoFeB 典型 0.005–0.05，目前 0.05 已偏高）
    5. R_W:   增大 → 单位 I_SOT 对应更高 V_SOT

本脚本对几组候选组合各跑短 MC（25 trials/point），输出 V_th 与
β_s 估计，便于挑选最接近实验的参数组合。
"""
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

from vgsot_sim.configs import PhysicalConstantsConfig, SerSotNoVcmaThermalConfig
from vgsot_sim.ser_cases import ser_sot_no_vcma_thermal

PULSE_NS, TOTAL_NS = 0.75, 4.0   # ← match experimental detailed-P_sw protocol
N_TRIALS = 30                     # cheap probe per candidate

I_LIST_uA = np.array([200, 400, 700, 1000, 1300, 1600, 2000, 2500, 3000])
I_LIST    = -I_LIST_uA * 1e-6     # negative for switching

def sigmoid_baseline(I_uA, V_th_uA, k_inv, baseline):
    """SER(|I_SOT|) ≈ baseline + (1 − baseline)/(1 + exp((|I| − V_th)/k)).

    Accounts for the back-hopping plateau (residual SER) at very high
    currents. As I → ∞ the SER → baseline, not 0.
    """
    return baseline + (1.0 - baseline) / (1.0 + np.exp((I_uA - V_th_uA) / k_inv))


def half_drop_threshold(I_uA, ser):
    """Interpolate the |I_SOT| where SER crosses (1 + plateau)/2.

    Robust to the residual back-hopping plateau when the simple-Sigmoid
    fit gets unstable: pick the highest-I plateau value (≈ tail mean)
    and find where SER first falls below the (1+plateau)/2 mid-level.
    """
    plateau = float(np.mean(ser[-3:]))
    mid = (1.0 + plateau) / 2.0
    # crossing from above mid → below mid as I increases
    for i in range(len(ser) - 1):
        if ser[i] >= mid and ser[i + 1] < mid:
            x0, x1 = I_uA[i], I_uA[i + 1]
            y0, y1 = ser[i], ser[i + 1]
            if y0 == y1:
                return float(x0)
            return float(x0 + (mid - y0) * (x1 - x0) / (y1 - y0))
    return float("nan")


def run_one_candidate(label, cc_overrides):
    cc = PhysicalConstantsConfig()
    for k, v in cc_overrides.items():
        setattr(cc, k, v)
    sim_end = int(round(TOTAL_NS * 1e-9 / cc.t_step))
    mid1    = int(round(PULSE_NS * 1e-9 / cc.t_step))
    cfg = SerSotNoVcmaThermalConfig(
        i_sot_list=tuple(I_LIST), trials=N_TRIALS,
        sim_start_step=1, sim_mid1_step=mid1, sim_end_step=sim_end,
        pap=1, non=1, vnv=0, v_mtj=0.0, r_sot_fl_dl=0.83,
        target_mz=1.0, failure_tol=0.2, constants=cc,
    )
    res = ser_sot_no_vcma_thermal(cfg, show_progress=False,
                                    enable_self_heating=False, seed=2026)
    # Robust V_th via mid-drop interpolation (handles back-hopping plateau).
    I_uA = -res.x * 1e6
    I_th_uA = half_drop_threshold(I_uA, res.ser)
    # Try a baseline-Sigmoid fit too for the slope (β_s) where possible.
    try:
        plateau = float(np.mean(res.ser[-3:]))
        p0 = (I_th_uA if np.isfinite(I_th_uA) else float(I_uA[len(I_uA)//2]),
              50.0, max(plateau, 0.02))
        popt, _ = curve_fit(
            lambda i, vth, k, b: sigmoid_baseline(i, vth, k, b),
            I_uA, res.ser, p0=p0,
            bounds=([10, 1, 0.0], [10000, 1000, 0.6]), maxfev=4000)
        _, k_inv, _ = popt
        beta_V_inv = (1.0 / k_inv) / cc.R_W * 1e6
    except Exception:
        k_inv, beta_V_inv = float("nan"), float("nan")
    V_th_mV = I_th_uA * cc.R_W * 1e-3              # µA × Ω → mV
    print(f"  {label:40s}: V_th(0.75ns) ≈ {V_th_mV:6.0f} mV  "
          f"(I_th = {I_th_uA:6.0f} µA),  β_s ≈ {beta_V_inv:6.1f} V⁻¹")
    print(f"     full SER:  ", "  ".join(f"{i:.0f}:{s:.2f}" for i, s in zip(I_uA, res.ser)))
    return I_th_uA, V_th_mV, beta_V_inv


# Scan brackets the FL-SOT-corrected (Cayley) calibration point θ_SH≈0.066.
# (The pre-fix spherical-Euler stepper required ≈0.04 for the same V_th; the
#  default integrator is now Cayley, so this scan uses it via the stepper default.)
CANDIDATES = [
    ("θ_SH=0.040 (pre-fix value)",      dict(theta_SH=0.040)),
    ("θ_SH=0.050",                      dict(theta_SH=0.050)),
    ("θ_SH=0.060",                      dict(theta_SH=0.060)),
    ("θ_SH=0.063",                      dict(theta_SH=0.063)),
    ("θ_SH=0.066 (calibrated)",         dict(theta_SH=0.066)),
    ("θ_SH=0.069",                      dict(theta_SH=0.069)),
    ("θ_SH=0.075",                      dict(theta_SH=0.075)),
]

EXP_V_TH_mV  = 894.0   # Device A P→AP @ 0.75 ns (chapter §2.3.3 detailed P_sw)
EXP_I_TH_uA  = EXP_V_TH_mV / 776 * 1000     # ≈ 1152
EXP_BETA_INV = 44.6                          # V⁻¹

print("=" * 90)
print(f"  Calibration scan — target Device A P→AP @ {PULSE_NS:.2f} ns:")
print(f"    V_th = {EXP_V_TH_mV:.0f} mV,  I_th = {EXP_I_TH_uA:.0f} µA,  β_s = {EXP_BETA_INV:.1f} V⁻¹")
print("=" * 90)
results = []
for label, overrides in CANDIDATES:
    I_th, V_th, beta = run_one_candidate(label, overrides)
    results.append((label, overrides, I_th, V_th, beta))
print()
print("Closest matches:")
ranked = sorted(results, key=lambda r: abs((r[3] or 0) - EXP_V_TH_mV))
for rec in ranked[:3]:
    label, ov, I_th, V_th, beta = rec
    err_pct = (V_th - EXP_V_TH_mV) / EXP_V_TH_mV * 100
    print(f"  {label:40s} V_th = {V_th:5.0f} mV ({err_pct:+5.1f}%)")
