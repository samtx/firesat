import sys, os

p = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(p)

import unittest
import numpy as np
import numpy.testing as npt
from firesat import timefn
import datetime


class Test_Timefn(unittest.TestCase):
    def shortDescription(self):
        return None

    def test_julian_date(self):
        """
        Vallado, eg.3-4
        """
        yr, mo, dy = 1996, 10.0, 26.0
        hr, mn, sec = 14.0, 20.0, 0.0
        jd = timefn.julian_date(yr, mo, dy, hr, mn, sec)
        jdT = 2450383.09722222
        npt.assert_almost_equal(jd, jdT, decimal=8)

    def test_julian_date2(self):
        """
        Vallado, eg. 11-5
        """
        yr, mo, dy = 1997, 4, 2
        hr, mn, sec = 1, 8, 0
        jd = timefn.julian_date(yr, mo, dy, hr, mn, sec)
        jdT = 2450540.5472
        npt.assert_almost_equal(jd, jdT, decimal=4)

    def test_julian_date_datetime(self):
        """
        Vallado, eg.3-4
        """
        yr, mo, dy = 1996, 10, 26
        hr, mn, sec = 14, 20, 0
        dt = datetime.datetime(yr, mo, dy, hr, mn, sec)
        jd = timefn.julian_date(dt)
        jdT = 2450383.09722222
        npt.assert_almost_equal(jd, jdT, decimal=8)

    def test_julian_date_datetime2(self):
        """
        Vallado, Eg. 5-1, p. 280
        """
        dt = datetime.datetime(2006, 4, 2)  # April 2, 2006, 00:00 UTC
        jdt = timefn.julian_date(dt)
        npt.assert_almost_equal(jdt, 2453827.5, decimal=12)

    # def test_julian_date_vectorized(self):
    #     """Use an array of datetimes to find the Julian Date"""
    #     dt_ary = np.arange(
    #         "2019-09-14T00:00:00", "2019-10-07T00:00:00", 200, dtype="datetime64"
    #     )
    #     jd_vectorized = np.vectorize(timefn.julian_date)
    #     jd_ary = jd_vectorized(dt_ary)
    #     # print(jd_ary)

    def test_jd_from_skyfield(self):
        """From skyfield.tests.test_earth_satellites.py"""
        ms = int(1e6) + 386000 - 439961
        dt = datetime.datetime(2004, 4, 6, 7, 51, 27, ms)
        jd = timefn.julian_date(dt)
        npt.assert_almost_equal(jd, 2453101.8274067827, decimal=12)

    def test_jd_from_skyfield2(self):
        """From skyfield.tests.test_earth_satellites.py"""
        ms = int(1e6) + 386000 - 439961
        dt = datetime.datetime(2004, 4, 6, 7, 51, 27, ms)
        jd = timefn.julian_date(dt)
        jd_desired = 2453101.8274067827
        npt.assert_almost_equal(jd, 2453101.8274067827, decimal=12)

    def test_jd_from_skyfield3(self):
        """From skyfield.tests.test_earth_satellites.py"""
        sec = 28.386 - 0.439961
        yr, mo, dy = 2004, 4, 6
        hr, mn = 7, 51
        jd = timefn.julian_date(yr, mo, dy, hr, mn, sec)
        jd_desired = 2453101.8274067827
        npt.assert_almost_equal(jd, jd_desired, decimal=12)

    def test_jday_to_invjday(self):
        yr, mo, dy = 1996, 10, 26
        hr, mn, sec = 14, 20, 1.123456
        jd = timefn.julian_date(yr, mo, dy, hr, mn, sec)
        out = timefn.invjday(jd)
        npt.assert_almost_equal(out[0], yr, decimal=2)
        npt.assert_almost_equal(out[1], mo, decimal=2)
        npt.assert_almost_equal(out[2], dy, decimal=2)
        npt.assert_almost_equal(out[3], hr, decimal=2)
        npt.assert_almost_equal(out[4], mn, decimal=2)
        npt.assert_almost_equal(out[5], sec, decimal=5)

    def test_jdy_floor_divide_scalar(self):
        # year, month, day = 1996, 10, 16
        # year, month, day = 2006, 4, 2
        year, month, day = 2000, 2, 29

        janfeb = month < 3

        part_A1 = 1461 * (year + 4800 - janfeb) // 4
        part_A2 = 367 * (month - 2 + janfeb * 12) // 12
        part_A3 = -(3 * ((year + 4900 - janfeb) // 100) // 4)
        out_A = day + part_A1 + part_A2 + part_A3 - 32075

        part_B1 = np.floor_divide(1461 * (year + 4800 - janfeb), 4)
        part_B2 = np.floor_divide(367 * (month - 2 + janfeb * 12), 12)
        part_B3 = -np.floor_divide(3 * np.floor_divide(year + 4900 - janfeb, 100), 4)
        out_B = day + part_B1 + part_B2 + part_B3 - 32075
        # fmtstr = ' {:10}  {:6}  {:6}'
        # print(fmtstr.format(part_A1, part_A2, part_A3))
        # print(fmtstr.format(part_B1, part_B2, part_B3))
        # print(month - 2 + janfeb * 12)
        npt.assert_almost_equal(out_A, out_B, verbose=True)

    def test_invjday_vectorized(self):
        jd = np.array(
            [
                2451723.53495062,
                2451723.78495062,
                2451724.03495062,
                2451724.28495062,
                2451724.53495062,
                2451724.78495062,
                2451725.03495062,
                2451725.28495062,
                2451725.53495062,
                2451725.78495062,
                2451726.03495062,
                2451726.28495062,
            ]
        )
        dt = timefn.invjday(jd)
        year, mon, day, hr, minute, sec = dt
        t_true = np.array(
            [
                [2000, 6, 28, 0, 50, 19.733571],
                [2000, 6, 28, 6, 50, 19.733571],
                [2000, 6, 28, 12, 50, 19.733571],
                [2000, 6, 28, 18, 50, 19.733571],
                [2000, 6, 29, 0, 50, 19.733571],
                [2000, 6, 29, 6, 50, 19.733571],
                [2000, 6, 29, 12, 50, 19.733571],
                [2000, 6, 29, 18, 50, 19.733571],
                [2000, 6, 30, 0, 50, 19.733571],
                [2000, 6, 30, 6, 50, 19.733571],
                [2000, 6, 30, 12, 50, 19.733571],
                [2000, 6, 30, 18, 50, 19.733571],
            ]
        )
        for i in range(len(jd)):
            npt.assert_almost_equal(year[i], t_true[i][0], decimal=1)
            npt.assert_almost_equal(mon[i], t_true[i][1], decimal=1)
            npt.assert_almost_equal(day[i], t_true[i][2], decimal=1)
            npt.assert_almost_equal(hr[i], t_true[i][3], decimal=1)
            npt.assert_almost_equal(minute[i], t_true[i][4], decimal=1)
            npt.assert_almost_equal(sec[i], t_true[i][5], decimal=6)
        #     print(f'{tsince[i]:6.1f}  {jd[i]:15.8f}  {year[i]:4.0f}  {mon[i]:2.0f}  {day[i]:2.0f}  {hr[i]:2.0f}:{minute[i]:2.0f}:{sec[i]:10.6f}')

    def test_jdt_tsince(self):
        jdt_start = timefn.julian_date(2000, 6, 28, 0, 50, 19.733571)
        tsince = np.linspace(0, 3960, 12)
        jd = timefn.jdt_tsince(jdt_start, tsince)
        jd_true = np.array(
            [
                2451723.53495062,
                2451723.78495062,
                2451724.03495062,
                2451724.28495062,
                2451724.53495062,
                2451724.78495062,
                2451725.03495062,
                2451725.28495062,
                2451725.53495062,
                2451725.78495062,
                2451726.03495062,
                2451726.28495062,
            ]
        )
        for i in range(len(jd)):
            npt.assert_almost_equal(jd[i], jd_true[i], decimal=8)


if __name__ == "__main__":

    # T = Test_Timefn()
    # T.test_julian_date_vectorized()

    test_cases = [Test_Timefn]
    suite = unittest.TestSuite()
    loader = unittest.TestLoader()
    for test_class in iter(test_cases):
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)

    unittest.TextTestRunner(verbosity=2).run(suite)
