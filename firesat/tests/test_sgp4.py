# Test spg4 extension

import numpy as np
import numpy.testing as npt
import unittest
import firesat.sgp4 as sgp4
from firesat.constants import DEG2RAD

class Test_SGP4(unittest.TestCase):

    def shortDescription(self):
        return None

    def test_vallado_example(self):
        """Revisting Spacetrack Report #3, 2006, TEME Example
        TLE:
        1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753
        2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667
        """
        satrec = sgp4.Satellite()
        satrec.whichconst = whichconst = sgp4.wgs72
        afspc_mode = False;  # 'improved' mode
        satrec.satnum = satnum = 5;
        satrec.epoch = epoch = 18441.784950620029;
        satrec.ecco = ecco = 0.1859667;
        satrec.argpo = argpo = 331.7664 * DEG2RAD;
        satrec.inclo = inclo = 34.2682 * DEG2RAD;
        satrec.mo = mo = 19.3264 * DEG2RAD;
        satrec.no = no = 10.82419157 / sgp4.min2rev;
        satrec.nodeo = nodeo = 348.7242 * DEG2RAD;
        satrec.ndot = ndot = 0.00000023;
        satrec.nddot = nddot = 0.0;
        satrec.bstar = bstar = 0.28098e-4;
        sgp4.sgp4init(whichconst, afspc_mode, satnum, epoch, bstar, ecco, argpo, inclo, mo, no, nodeo, satrec)

        t = np.linspace(0, 4320, 13)   # minutes after epoch
        tsize = len(t)
        r = np.zeros((tsize, 3))
        v = np.zeros((tsize, 3))

        for i in range(tsize):
            ti = t[i]
            ri, vi = sgp4.sgp4(satrec, ti)
            r[i] = ri
            v[i] = vi

        true_rv = np.array([
            [ 7022.46529266, -1400.08296755,     0.03995155,  1.893841015,  6.405893759,  4.534807250],
            [-7154.03120202, -3783.17682504, -3536.19412294,  4.741887409, -4.151817765, -2.093935425],
            [-7134.59340119,  6531.68641334,  3260.27186483, -4.113793027, -2.911922039, -2.557327851],
            [ 5568.53901181,  4492.06992591,  3863.87641983, -4.209106476,  5.159719888,  2.744852980],
            [ -938.55923943, -6268.18748831, -4294.02924751,  7.536105209, -0.427127707,  0.989878080],
            [-9680.56121728,  2802.47771354,   124.10688038, -0.905874102, -4.659467970, -3.227347517],
            [  190.19796988,  7746.96653614,  5110.00675412, -6.112325142,  1.527008184, -0.139152358],
            [ 5579.55640116, -3995.61396789, -1518.82108966,  4.767927483,  5.123185301,  4.276837355],
            [-8650.73082219, -1914.93811525, -3007.03603443,  3.067165127, -4.828384068, -2.515322836],
            [-5429.79204164,  7574.36493792,  3747.39305236, -4.999442110, -1.800561422, -2.229392830],
            [ 6759.04583722,  2001.58198220,  2783.55192533, -2.180993947,  6.402085603,  3.644723952],
            [-3791.44531559, -5712.95617894, -4533.48630714,  6.668817493, -2.516382327, -0.082384354],
            [-9060.47373569,  4658.70952502,   813.68673153, -2.232832783, -4.110453490, -3.157345433]
        ])
        for i in range(r.shape[0]):
            # print_rv(t[i], r[i], v[i])
            npt.assert_allclose(r[i], true_rv[i,[0,1,2]], rtol=0, atol=1e-8, verbose=True)
            npt.assert_allclose(v[i], true_rv[i,[3,4,5]], rtol=0, atol=1e-8, verbose=True)


    def test_vallado_example_2(self):
        """Revisting Spacetrack Report #3, 2006, XM-3 example
        TLE:
        1 28626U 05008A  06176.46683397 -.00000205 00000-0  10000-3 0 2190
        2 28626  0.0019 286.9433 0000335  13.7918 55.6504  1.00270176 4891
        """
        satrec = sgp4.Satellite()
        satrec.whichconst = whichconst = sgp4.wgs72
        afspc_mode = False;  # 'improved' mode
        satrec.satnum = satnum = 28626;
        satrec.epoch = epoch = 20630.466833970044;
        satrec.ecco = ecco = 0.0000335;
        satrec.argpo = argpo = 13.7918 * DEG2RAD;
        satrec.inclo = inclo = 0.0019 * DEG2RAD;
        satrec.mo = mo = 55.6504 * DEG2RAD;
        satrec.no = no = 1.00270176 / sgp4.min2rev;
        satrec.nodeo = nodeo = 286.9433 * DEG2RAD;
        satrec.ndot = ndot = -.00000205;
        satrec.nddot = nddot = 0.0;
        satrec.bstar = bstar = 0.10000-3;
        sgp4.sgp4init(whichconst, afspc_mode, satnum, epoch, bstar, ecco, argpo, inclo, mo, no, nodeo, satrec)

        t = np.linspace(0, 1440, 13)  # minutes after epoch
        tsize = len(t)
        r = np.zeros((tsize, 3))
        v = np.zeros((tsize, 3))

        for i in range(tsize):
            ti = t[i]
            ri, vi = sgp4.sgp4(satrec, ti)
            r[i] = ri
            v[i] = vi

        true_rv = np.array([
            [ 42080.71852213,  -2646.86387436,  0.81851294,  0.193105177,  3.068688251,  0.000438449],
            [ 37740.00085593,  18802.76872802,  3.45512584, -1.371035206,  2.752105932,  0.000336883],
            [ 23232.82515008,  35187.33981802,  4.98927428, -2.565776620,  1.694193132,  0.000163365],
            [  2467.44290178,  42093.60909959,  5.15062987, -3.069341800,  0.179976276, -0.000031739],
            [-18962.59052991,  37661.66243819,  4.04433258, -2.746151982, -1.382675777, -0.000197633],
            [-35285.00095313,  23085.44402778,  2.08711880, -1.683277908, -2.572893625, -0.000296282],
            [-42103.20138132,   2291.06228893, -0.13274964, -0.166974816, -3.070104560, -0.000311007],
            [-37580.31858370, -19120.40485693, -2.02755702,  1.394367848, -2.740341612, -0.000248591],
            [-22934.20761876, -35381.23870806, -3.16495932,  2.580167539, -1.672360951, -0.000134907],
            [ -2109.90332389, -42110.71508198, -3.36507889,  3.070935369, -0.153808390, -0.000005855],
            [ 19282.77774728, -37495.59250598, -2.71861462,  2.734400524,  1.406220933,  0.000103486],
            [ 35480.60990600, -22779.03375285, -1.52841859,  1.661210676,  2.587414593,  0.000168300],
            [ 42119.96263499,  -1925.77567263, -0.19827433,  0.140521206,  3.071541613,  0.000179561]
        ])

        # print('xxx')
        for i in range(tsize):
            # print_rv(t[i], r[i], v[i])
            npt.assert_allclose(r[i], true_rv[i,[0,1,2]], rtol=0, atol=1e-4, verbose=True)
            npt.assert_allclose(v[i], true_rv[i,[3,4,5]], rtol=0, atol=1e-4, verbose=True)


def print_rv(t, r, v):
    print(f't={t:8.2f}  r={r[0]:14.6f}  {r[1]:14.6f}  {r[2]:14.6f}  |r|={np.linalg.norm(r):8.5f}   ',end='')
    print(f'v={v[0]:10.6f}  {v[1]:10.6f}  {v[2]:10.6f} |v|={np.linalg.norm(v):.6}')


def main():
    tend = 1440*4
    t = np.linspace(0, tend, tend+1)
    tsize = t.size

    DEG2RAD = 0.017453292519943295474371680598    # pi/180
    MIN2REV = 229.183118052329291458590887486935  # 1440/(2*pi)
    radpersec_to_revperday = MIN2REV * 60  # 86400/(2*pi)
    mu = 3.986*(10**14)
    H = 1.7e7   # altitude [m]
    Re = 6378140  # Earth radius [m]

    # propagate using sgp4
    satrec = sgp4.Satellite()
    satrec.no = np.sqrt(mu/((Re+H)**3)) * 60
    sgp4.sgp4init(satrec.whichconst, False, satrec.satnum, satrec.epoch, satrec.bstar, satrec.ecco,
                    satrec.argpo, satrec.inclo, satrec.mo, satrec.no, satrec.nodeo, satrec)

    r = np.zeros((tsize, 3))
    v = np.zeros((tsize, 3))

    for i in range(tsize):
        ti = t[i]
        ri, vi = sgp4.sgp4(satrec, ti)
        r[i] = ri
        v[i] = vi

    for i in range(10):
        print_rv(t[i], r[i], v[i])

    # compute orbit time
    # find indecies where y-axis value goes from negative to positive, then linearly interpolate for times
    # this function finds the sign change, and where r2[i+1]>r2[i]
    idx = np.nonzero((r[:-1,1]*r[1:,1]<0) & (r[:-1,1]<r[1:,1]))[0]
    t0 = np.diff(t)[idx]/np.diff(r[:,1])[idx] * (-r[idx,1]) + t[idx]
    # print(t0)
    # print('crossing times')
    # for k, i in enumerate(idx):
    #     print_rv(t[i], r[i], v[i])
    #     print_rv(t[i+1], r[i+1], v[i+1])
    #     print(f't orbit at = {t0[k]:.3f}')

    dt_orbit = np.mean(np.diff(t0)) * 60
    v_lofi = np.sqrt(mu/(Re + H))
    dt_orbit_lofi = 2*np.pi*(Re + H)/v_lofi
    print(f'dt_orbit hifi = {dt_orbit:.3f} sec')
    print(f'dt_orbit lofi = {dt_orbit_lofi:.3f} sec')
    dtdiff = abs(dt_orbit - dt_orbit_lofi)
    print(f'abs diff = {dtdiff:.6f},  rel diff ={dtdiff/dt_orbit:0.6e}')

    # compute average velocity
    vnorm = np.linalg.norm(v, axis=1)
    vavg = np.mean(vnorm) * 1000
    vdiff = abs(vavg - v_lofi)
    print(f'velocity hifi = {vavg:0.6f}')
    print(f'velocity lofi = {v_lofi:0.6f}')
    print(f'abs diff = {vdiff:.6f},  rel diff ={vdiff/vavg:0.6e}')


if __name__ == "__main__":
    # main()

    loader = unittest.TestLoader()
    tests = loader.loadTestsFromTestCase(Test_SGP4)
    suite = unittest.TestSuite()
    suite.addTests(tests)
    unittest.TextTestRunner(verbosity=2).run(suite)
