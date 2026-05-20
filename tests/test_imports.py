"""Smoke test for the public-API import path.

Verifies the package imports cleanly and that the core entry points
return finite results when supplied with the default `PhysicalConstantsConfig`.
"""
import json
import subprocess
import sys
import textwrap


def test_imports():
    import vgsot_sim                                            # noqa: F401
    from vgsot_sim.configs import PhysicalConstantsConfig
    from vgsot_sim.initialize import init, compute_Rp
    from vgsot_sim.tmr import tmr

    cc = PhysicalConstantsConfig()
    r, theta, mz, phi = init(1, cc)
    assert r > 0 and -1.5 < mz <= 1.0
    r_p = compute_Rp(cc)
    assert r_p > 0
    r_at_mz = tmr(0.0, mz, cc)
    assert r_at_mz > 0


def test_import_preserves_matplotlib_rcparams():
    """Importing the package must not override a script's plotting style."""
    code = textwrap.dedent(
        """
        import json
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.rcParams["font.family"] = ["sans-serif"]
        plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans"]
        plt.rcParams["font.size"] = 9

        import vgsot_sim  # noqa: F401

        print(json.dumps({
            "font.family": list(plt.rcParams["font.family"]),
            "font.sans-serif": list(plt.rcParams["font.sans-serif"][:2]),
            "font.size": plt.rcParams["font.size"],
        }))
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        text=True,
        capture_output=True,
    )
    rc = json.loads(result.stdout)
    assert rc["font.family"] == ["sans-serif"]
    assert rc["font.sans-serif"] == ["Arial", "Liberation Sans"]
    assert rc["font.size"] == 9.0


if __name__ == "__main__":
    test_imports()
    print("IMPORT SUCCESS!")
