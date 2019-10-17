import unittest
import numpy as np
import numpy.testing as npt
from firesat import orbit
import firesat.solar as solar
import firesat.timefn as timefn
import firesat.utils as utils
import firesat.system as system
import datetime
from firesat.constants import AU_KM, DEG2RAD, RAD2DEG

class Test_Orbit(unittest.TestCase):

    def shortDescription(self):
        return None

    def test_orbit_lowfidelity(self):
        np.random.seed(1234)
        sat_params = system.setup()
        n = 10000
        x = utils.mvn(['H', 'phi'], n)
        qoi = orbit(x, sat_params, fidelity=0)
        qoi_means_true = np.array([
            4045.81892, 37913.0974, 3191.54502, 1.30984041e-5
        ])
        qoi_means = np.mean(qoi, axis=1)
        for i, q in enumerate(qoi_means):
            npt.assert_approx_equal(q, qoi_means_true[i])

    def test_orbit_medfidelity(self):
        np.random.seed(1234)
        sat_params = system.setup()
        n = 10
        x = utils.mvn(['H', 'phi'], n)
        qoi = orbit(x, sat_params, fidelity=1)
        qoi_means = np.mean(qoi, axis=1)
        qoi_means_true = np.array([
            3999.73915, 39189.6937, 6291.00484, 1.23764331e-5
        ])
        for i, q in enumerate(qoi_means):
            npt.assert_approx_equal(q, qoi_means_true[i])

np.random.seed(1234)
sat_params = system.setup()
n = 10
x = utils.mvn(['H', 'phi'], n)

def prof_orbit():
    orbit(x, sat_params)

if __name__ == "__main__":

    import cProfile

    np.random.seed(1234)
    sat_params = system.setup()
    n = 100
    x = utils.mvn(['H', 'phi'], n)
    # print(x)
    pf = cProfile.Profile()
    pf.enable()
    orbit(x, sat_params, fidelity=1)
    pf.disable()
    pf.dump_stats('orbit_prof_cython.pf')

    # T = Test_Orbit()
    # T.test_orbit_medfidelity()

    # suite = unittest.TestSuite()
    # loader = unittest.TestLoader()
    # tests = loader.loadTestsFromTestCase(Test_Orbit)
    # suite.addTests(tests)
    # unittest.TextTestRunner(verbosity=2).run(suite)
