"""
实验 F — PART A：图表积分器溯源审计（命题 C4）。

目的
----
正文 §2.2.3.2 与图脚本注释之间存在历史矛盾：图脚本 `plot_ser_mc.py` 的
`--integrator` 帮助文本仍写「euler_spherical (default, matches the θ_SH=0.04
experimental calibration)」，`plot_single_trajectory.py` 的注释亦提到「θ_SH=0.04
… spherical-Euler」。这批注释是**陈旧的**：发布代码 `run_piecewise_direct_excitation`
的 `integrator` 默认值实为 `'cayley'`，而图 2.9 的脚本恰恰不显式传 `integrator=`，
因此走的是 Cayley 向量步，而非球坐标-Euler。

本脚本用一个 ~10 行的上下文管理器，对 `time_series_cases` 命名空间里实际被调用的
两个步进函数 (`switching_vector` = Cayley 向量步；`switching` = 球坐标-Euler 步)
打桩计数，无头跑两个图脚本，记录 (n_cayley, n_euler)；并用 `inspect.signature`
读回 `run_piecewise_direct_excitation` 的 `integrator` 默认值。

决定性判据 (C4)
----------------
图 2.9 (`plot_single_trajectory.py`)：n_euler == 0 且 n_cayley > 0，且
默认 integrator == 'cayley'  →  「默认 Euler-球坐标」的图脚本注释陈旧、被证伪。

只调用公共 API + 标准库 monkeypatch；不修改 src/。matplotlib 用 Agg 后端。
运行：  PYTHONPATH=src python scripts/09_simulation_figures/integrator_audit.py
"""
from __future__ import annotations

import contextlib
import inspect
import io
import json
import runpy
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless — figure scripts call plt.savefig/close

# Force UTF-8 on our own stdout so the Chinese log lines survive a GBK console.
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import vgsot_sim.time_series_cases as tsc

REPO_ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = Path(__file__).resolve().parent
OUT_DIR = REPO_ROOT / "result" / "sec_2_2_3_2" / "F"
OUT_DIR.mkdir(parents=True, exist_ok=True)


@contextlib.contextmanager
def count_integrator_calls():
    """Monkeypatch the two LLG steppers *as resolved inside time_series_cases*
    with call-counting wrappers, restoring the originals on exit.

    `run_piecewise_*` call the names `switching_vector` / `switching` that live
    in the `time_series_cases` module globals (imported at module load), so the
    counters must wrap those bound names — not the source-module attributes.
    `switching_vector` is the Cayley vector step; `switching` is the legacy
    spherical-Euler step.
    """
    counts = {"n_cayley": 0, "n_euler": 0}
    orig_vec, orig_sph = tsc.switching_vector, tsc.switching

    def vec_wrapper(*a, **k):
        counts["n_cayley"] += 1
        return orig_vec(*a, **k)

    def sph_wrapper(*a, **k):
        counts["n_euler"] += 1
        return orig_sph(*a, **k)

    tsc.switching_vector, tsc.switching = vec_wrapper, sph_wrapper
    try:
        yield counts
    finally:
        tsc.switching_vector, tsc.switching = orig_vec, orig_sph


def run_figure_script(script_path: Path, argv_extra=None):
    """Execute a figure script headless under counting, returning its counts.

    Uses runpy.run_path so the script body runs exactly as `python <script>`
    would (top-level statements drive the simulation), with sys.argv set so any
    argparse defaults inside the script are honoured (we pass NO --integrator,
    so each script picks its own default — that is the whole point of the audit).
    """
    saved_argv = sys.argv[:]
    sys.argv = [str(script_path)] + (argv_extra or [])
    try:
        with count_integrator_calls() as counts:
            # Swallow the script's stdout chatter into a UTF-8 StringIO sink
            # (a text devnull would inherit the Windows GBK console codec and
            # choke on the µ / Δ glyphs the figure scripts print).
            with contextlib.redirect_stdout(io.StringIO()):
                runpy.run_path(str(script_path), run_name="__main__")
        return dict(counts)
    finally:
        sys.argv = saved_argv


def introspect_default_integrator():
    sig = inspect.signature(tsc.run_piecewise_direct_excitation)
    return sig.parameters["integrator"].default


def fig29_passes_integrator_kw(script_path: Path) -> bool:
    """True if plot_single_trajectory.py passes an explicit `integrator=` kw
    to run_piecewise_direct_excitation. It must NOT — so it inherits the
    'cayley' driver default. (Static check complementing the runtime counters.)
    """
    src = script_path.read_text(encoding="utf-8")
    return "integrator=" in src


def main():
    fig29 = FIG_DIR / "plot_single_trajectory.py"   # 图 2.9
    fig211 = FIG_DIR / "plot_ser_mc.py"              # 图 2.11

    default_integrator = introspect_default_integrator()

    print("=" * 78)
    print("  实验 F — PART A：图表积分器溯源 (C4)")
    print("=" * 78)
    print(f"  run_piecewise_direct_excitation  integrator 默认值 = "
          f"{default_integrator!r}")
    print("-" * 78)

    audit = {}

    # --- 图 2.9 (decisive target) ---------------------------------------
    print(f"  跑 {fig29.name} (图 2.9, 无 --integrator) ...")
    c29 = run_figure_script(fig29)
    print(f"    n_cayley = {c29['n_cayley']}   n_euler = {c29['n_euler']}")
    fig29_overrides = fig29_passes_integrator_kw(fig29)
    audit["fig_2.9_plot_single_trajectory"] = {
        "script": str(fig29.relative_to(REPO_ROOT)).replace("\\", "/"),
        "n_cayley": c29["n_cayley"],
        "n_euler": c29["n_euler"],
        "integrator_used": ("cayley" if c29["n_euler"] == 0 and c29["n_cayley"] > 0
                            else ("euler_spherical" if c29["n_cayley"] == 0
                                  and c29["n_euler"] > 0 else "mixed/none")),
        "script_passes_integrator_kw": bool(fig29_overrides),
        "passes_driver_default": not fig29_overrides,  # omits integrator= → default
    }

    # --- 图 2.11 (context: its own default is euler_spherical) ------------
    print(f"  跑 {fig211.name} (图 2.11, 无 --integrator → 脚本自带默认) ...")
    c211 = run_figure_script(fig211)
    print(f"    n_cayley = {c211['n_cayley']}   n_euler = {c211['n_euler']}")
    audit["fig_2.11_plot_ser_mc"] = {
        "script": str(fig211.relative_to(REPO_ROOT)).replace("\\", "/"),
        "n_cayley": c211["n_cayley"],
        "n_euler": c211["n_euler"],
        "integrator_used": ("cayley" if c211["n_euler"] == 0 and c211["n_cayley"] > 0
                            else ("euler_spherical" if c211["n_cayley"] == 0
                                  and c211["n_euler"] > 0 else "mixed/none")),
        "note": ("plot_ser_mc.py argparse default is 'euler_spherical' and it "
                 "forwards integrator=INTEGRATOR explicitly; it is NOT fig 2.9."),
    }

    # --- decisive verdict for C4 -----------------------------------------
    fig29_is_cayley = (c29["n_euler"] == 0 and c29["n_cayley"] > 0)
    default_is_cayley = (default_integrator == "cayley")
    c4_confirmed = fig29_is_cayley and default_is_cayley

    result = {
        "experiment": "F",
        "part": "A_provenance",
        "claim": "C4",
        "claim_text_zh": ("图 2.9 究竟用哪个积分器；图脚本注释「默认 "
                          "Euler-球坐标 / θ_SH=0.04」是否陈旧"),
        "run_piecewise_direct_excitation_default_integrator": default_integrator,
        "audit": audit,
        "decisive_criterion": ("fig 2.9  n_euler==0 且 n_cayley>0  且  "
                               "default integrator == 'cayley'"),
        "fig_2.9_uses_cayley": bool(fig29_is_cayley),
        "default_integrator_is_cayley": bool(default_is_cayley),
        "C4_stale_comment_confirmed": bool(c4_confirmed),
        "verdict": ("STALE-COMMENT CONFIRMED: fig 2.9 runs Cayley, driver "
                    "default is 'cayley'; the '默认 Euler-球坐标' figure-script "
                    "comment is stale."
                    if c4_confirmed else
                    "NOT confirmed — see counts above."),
    }

    print("-" * 78)
    print(f"  C4 stale-comment confirmed : {c4_confirmed}")
    print("=" * 78)

    out_json = OUT_DIR / "integrator_provenance.json"
    out_json.write_text(json.dumps(result, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(f"  wrote {out_json}")
    return result


if __name__ == "__main__":
    main()
