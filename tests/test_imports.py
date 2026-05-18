"""Smoke test for the public-API import path.

Verifies the package imports cleanly and that the core entry points
return finite results when supplied with the default `PhysicalConstantsConfig`.
"""
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


if __name__ == "__main__":
    test_imports()
    print("IMPORT SUCCESS!")
