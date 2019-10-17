import sys, os
p = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(p)

import unittest
import numpy as np
import numpy.testing as npt
import firesat.solar as solar
import firesat.timefn as timefn
import datetime
from firesat.constants import AU_KM, DEG2RAD, RAD2DEG

class Test_Solar(unittest.TestCase):

    def shortDescription(self):
        return None

    def test_sun_pos(self):
        """
        Vallado, Eg. 5-1, p. 280
        """
        jdt = 2453827.5  # April 2, 2006, 00:00 UTC
        jdt = np.asarray(jdt)
        r = solar.sun_pos(jdt)
        r_true = np.array([146186212, 28788976, 12481064], dtype=np.float)
        npt.assert_allclose(r, r_true, rtol=1e-4)

    def test_sun_pos2(self):
        """
        Vallado, Eg. 11-6, p. 913
        """
        dt = datetime.datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
        jdt = timefn.julian_date(dt)
        jdt = np.asarray(jdt)
        r = solar.sun_pos(jdt) / AU_KM
        r_true = np.array([0.9765, 0.1960, 0.0850])
        npt.assert_allclose(r, r_true, rtol=1e-3)

    def test_sun_pos3(self):
        """
        Vallado, Eg. 5-1, p. 280
        Vallado, Eg. 11-6, p. 913
        """
        # April 2, 1997, 01:08:0.00 UTC, April 2, 2006, 00:00 UTC
        jdt = np.array([2450540.54722222, 2453827.5])
        r = solar.sun_pos(jdt) / AU_KM
        r_true0 = np.array([0.9765, 0.1960, 0.0850])
        r_true1 = np.array([0.97719447, 0.19244242, 0.08343076])
        npt.assert_allclose(r[0], r_true0, rtol=1e-3)
        npt.assert_allclose(r[1], r_true1, rtol=1e-4)

    def test_sun_sat_angle(self):
        """
        Vallado, Eg. 11-6, p.913
        """
        dt = datetime.datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
        jdt = timefn.julian_date(dt)
        jdt = np.asarray(jdt)
        rsun = solar.sun_pos(jdt)
        rsat = np.array([-2811.2769, 3486.2632, 5069.5763])
        rsun = np.atleast_2d(rsun)
        rsat = np.atleast_2d(rsat)
        sunangle = solar.sun_sat_angle(rsat, rsun) * RAD2DEG
        npt.assert_almost_equal(sunangle, 76.0407, decimal=3)

    def test_sun_sat_angle2(self):
        """
        Vallado, Eg. 11-6, p.913
        """
        rsat = np.array([-2811.2769, 3486.2632, 5069.5763])
        rsun = np.array([0.9765, 0.1960, 0.0850]) * AU_KM
        sunangle = solar.sun_sat_angle(rsat, rsun) * RAD2DEG
        npt.assert_almost_equal(sunangle, 76.0407, decimal=3)


    def test_sun_sat_orthogonal_distance(self):
        """
        Vallado, Eg. 11-6, p.913
        """
        r = np.array([[-2811.2769, 3486.2632, 5069.5763]])  # sat, ECI coordinates
        zeta = 76.0407  # deg
        dist = solar.sun_sat_orthogonal_distance(r, zeta * DEG2RAD)
        npt.assert_almost_equal(dist, 6564.6870, decimal=4)

    def test_is_sat_illuminated(self):
        """
        Vallado, Eg. 11-6, p.913
        """
        dt = datetime.datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
        jdt = np.array([timefn.julian_date(dt)])
        rsat = np.array([[-2811.2769, 3486.2632, 5069.5763]])  # ECI coords
        rsun = solar.sun_pos(jdt)
        vis = solar.is_sat_illuminated(rsat, rsun)
        assert vis

    def test_is_sat_illuminated2(self):
        jdt = np.array([2450540.54722222, 2453827.5])
        rsun = solar.sun_pos(jdt)  # (3, 2) array
        rsat = np.array([
            [-2811.2769, 3486.2632, 5069.5763],
            [-2811.2769, 3486.2632, 5069.5763]
        ])  # ECI coords
        vis = solar.is_sat_illuminated(rsat, rsun)
        for i in range(len(vis)):
            assert vis[i]


if __name__ == "__main__":

    suite = unittest.TestSuite()
    loader = unittest.TestLoader()
    tests = loader.loadTestsFromTestCase(Test_Solar)
    suite.addTests(tests)
    unittest.TextTestRunner(verbosity=2).run(suite)
