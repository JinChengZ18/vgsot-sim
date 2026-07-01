"""实验 F 回归测试 — 图表与参数溯源（命题 C4 / C13）。

锁定 §2.2.3.2 的两条溯源结论，使图脚本注释/正文中残留的「默认 Euler-球坐标 /
θ_SH=0.04」陈述无法在不被测试捕获的情况下重新成为「真」：

  C13: PhysicalConstantsConfig().theta_SH == 0.066 （而非 0.04）。
  C4 : run_piecewise_direct_excitation 的 integrator 默认值 == 'cayley'，
       且该默认路径实际调用的是 Cayley 向量步 switching_vector（n_euler==0）。

测试只用公共 API + 标准库 monkeypatch；不依赖任何重型蒙特卡罗系综。
"""
import inspect

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
import vgsot_sim.time_series_cases as tsc
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation


def test_theta_sh_is_0066_not_004():
    """C13: runtime θ_SH is the recalibrated 0.066, not the stale 0.04."""
    cc = PhysicalConstantsConfig()
    assert cc.theta_SH == 0.066, (
        f"theta_SH={cc.theta_SH!r}; the figure-caption '0.04' value is stale "
        "and must not be the runtime default."
    )


def test_driver_default_integrator_is_cayley():
    """C4: the figure driver's default integrator is 'cayley'."""
    default = inspect.signature(
        run_piecewise_direct_excitation
    ).parameters["integrator"].default
    assert default == "cayley", (
        f"default integrator={default!r}; figures that omit integrator= run "
        "this default, so it must be 'cayley' (not 'euler_spherical')."
    )


def test_default_path_routes_to_cayley_vector_step():
    """C4: with the default integrator, the driver calls switching_vector
    (Cayley) and never switching (spherical-Euler).

    Monkeypatch the two step functions *as bound in time_series_cases* (where
    they are actually called) with counters, run a short pulse, and assert the
    routing.
    """
    counts = {"cayley": 0, "euler": 0}
    orig_vec, orig_sph = tsc.switching_vector, tsc.switching

    def vec(*a, **k):
        counts["cayley"] += 1
        return orig_vec(*a, **k)

    def sph(*a, **k):
        counts["euler"] += 1
        return orig_sph(*a, **k)

    tsc.switching_vector, tsc.switching = vec, sph
    try:
        cc = PhysicalConstantsConfig()
        mid1 = 20
        sim_end = 40
        run_piecewise_direct_excitation(
            sim_start_step=1, sim_mid1_step=mid1, sim_mid2_step=sim_end,
            sim_end_step=sim_end, pap=1,
            v_mtj_stage1=0.0, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
            i_sot_stage1=-1300e-6, i_sot_stage2=0.0, i_sot_stage3=0.0,
            estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1,
            estt_stage3=0, esot_stage3=1,
            vnv=0, non=1, r_sot_fl_dl=0.83,
            constants=cc, show_progress=False, rng=np.random.default_rng(0),
            # integrator omitted → default 'cayley'
        )
    finally:
        tsc.switching_vector, tsc.switching = orig_vec, orig_sph

    assert counts["cayley"] == sim_end and counts["euler"] == 0, (
        f"default-integrator routing wrong: cayley={counts['cayley']}, "
        f"euler={counts['euler']} (expected {sim_end}, 0)."
    )


if __name__ == "__main__":
    test_theta_sh_is_0066_not_004()
    test_driver_default_integrator_is_cayley()
    test_default_path_routes_to_cayley_vector_step()
    print("FIGURE PROVENANCE TESTS PASSED")
