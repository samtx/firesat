# Test atmosphere model

import unittest
import numpy as np
import numpy.testing as npt
import firesat.atmosphere as atmos

class Test_Atmosphere(unittest.TestCase):

    def shortDescription(self):
        return None

    def test_exponential_atmosphere_scalar(self):
        h = 712345  # [m] altitude
        h /= 1000   # convert to [km]
        rho = atmos.exponential_density_model(h)
        npt.assert_approx_equal(rho, 3.144284600e-14 )

    def test_exponential_atmosphere_vector(self):
        # [m] altitude
        # normally distributed, mean=700000[m], std=100000[m]
        h = np.array([
            747143.52,
            580902.43,
            843270.70,
            668734.81,
            627941.13,
            788716.29,
            785958.84,
            636347.65,
            701569.64,
            475731.50,
        ], dtype=float)
        h /= 1000.0   # convert to [km]
        rho = atmos.exponential_density_model(h)
        rho_true = np.array([
            2.12362259e-14,
            1.96120416e-13,
            8.26825457e-15,
            5.58486669e-14,
            9.85462508e-14,
            1.32877737e-14,
            1.37075025e-14,
            8.76630548e-14,
            3.55058573e-14,
            1.03827812e-12,
        ])
        for i in range(len(h)):
            npt.assert_approx_equal(rho[i], rho_true[i])


if __name__ == "__main__":

    suite = unittest.TestSuite()
    loader = unittest.TestLoader()
    tests = loader.loadTestsFromTestCase(Test_Atmosphere)
    suite.addTests(tests)
    unittest.TextTestRunner(verbosity=2).run(suite)
