import unittest
import numpy as np
import firesat.system as system
import firesat.utils as utils
from firesat import power

class Test_Firesat(unittest.TestCase):

    def shortDescription(self):
        return None

    def setUp(self):
        np.random.seed(42)

    def test_firesat_system(self):
        # Test Fire Satellite System
        n = int(1e5)
        usehifi = False
        sat_params = system.setup()
        input_vars = sat_params['rand_inputs']
        x = utils.mvn(input_vars, n)
        qoi = system.run(x, var_info=sat_params)
        self.assertTrue(True)

    def test_firesat_power_lowfidelity(self):
        # Test power model, low fidelity
        n = int(1e3)
        fidelity = 0
        sat_params = system.setup()
        input_vars = ['Po', 'F_s']
        x = utils.mvn(input_vars, n)
        cpl_vars = ['PACS', 'dt_orbit', 'dt_eclipse']
        y = utils.uniform(cpl_vars, n)
        qoi = power(x, y, var_info=sat_params, fidelity=fidelity)
        self.assertTrue(True)


if __name__ == "__main__":

    suite = unittest.TestSuite()
    loader = unittest.TestLoader()
    tests = loader.loadTestsFromTestCase(Test_Firesat)
    suite.addTests(tests)
    unittest.TextTestRunner(verbosity=2).run(suite)
