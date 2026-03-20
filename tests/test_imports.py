def test_imports():
    import vgsot_sim
    from vgsot_sim.initialize import init
    from vgsot_sim.tmr import tmr
    r, theta, mz, phi = init(1)
    assert r > 0
    r2 = tmr(0.0, mz)
    assert r2 > 0
    print('IMPORT SUCESS!')



if __name__ == "__main__":
    test_imports()
