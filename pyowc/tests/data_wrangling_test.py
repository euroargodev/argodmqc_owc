import os
import unittest
import numpy as np
from scipy.io import loadmat

from pyowc.data.wrangling import interp_climatology, map_data_grid
from . import TESTS_CONFIG


# pylint: disable=too-many-instance-attributes
class InterpClimatology(unittest.TestCase):
    """
    Test cases for interp_climatology function
    """

    def setUp(self):
        """
        Set up some values for testing, pulling from matlab matrices.
        The matlab matrices were used in the matlab version of this code, so answers should match
        :return: Nothing
        """
        # load in the data for testing
        test = loadmat(os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "interp_climatology", "testfile.mat"]))
        results = loadmat(os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "interp_climatology", "results.mat"]))

        # set the test variables from the loaded .mat file
        self.grid_sal = test['S']
        self.grid_theta = test['Theta']
        self.grid_pres = test['P']
        self.float_sal = test['S_f']
        self.float_theta = test['Theta_f']
        self.float_pres = test['P_f']

        self.expected_interp_pres = results['P_h']
        self.expected_interp_sal = results['S_h']

    def test_throws_error(self):
        """
        Test that an we get NaNs for bad data
        :return: Nothing
        """
        print("Testing that interp_climatology gives NaNs if all data is bad")

        bad_grid_sal = np.full((5, 5), np.inf)

        sal, pres = interp_climatology(bad_grid_sal, self.grid_theta, self.grid_pres,
                                       self.float_sal, self.float_theta, self.float_pres)

        self.assertTrue(np.all(np.isnan(sal)), True)
        self.assertTrue(np.all(np.isnan(pres)), True)

    def test_returns_correct_shape(self):
        """
        Test that the returned matrix is the shape we expect it to be
        :return: Nothing
        """
        print("Testing that interp_climatology returns a matrix of the correct shape")

        sal, pres = interp_climatology(self.grid_sal, self.grid_theta, self.grid_pres,
                                       self.float_sal, self.float_theta, self.float_pres)

        self.assertTrue(sal.shape == self.expected_interp_sal.shape,
                        "salinity matrix shape is incorrect")
        self.assertTrue(pres.shape == self.expected_interp_pres.shape,
                        "pressure matrix shape is incorrect")

    def test_returns_same_shape(self):
        """
        Test that the salinity and pressures matrices are equal shapes
        :return: Nothing
        """
        print("Testing that interp_climatology returns a matrices of the same shape")

        sal, pres = interp_climatology(self.grid_sal, self.grid_theta, self.grid_pres,
                                       self.float_sal, self.float_theta, self.float_pres)

        self.assertTrue(sal.shape == pres.shape,
                        "salinity and pressure matrices shaped differently")

    def test_returns_correct_values(self):
        """
        Test that the output of this function matches the matlab version value by value.
        We cannot compare NaNs directly (as it will always give false, so only do a direct
        comparison if they are both not NaN
        :return: Nothing
        """
        print("Testing that interp_climatology returns matrices with correct values")

        sal, pres = interp_climatology(self.grid_sal, self.grid_theta, self.grid_pres,
                                       self.float_sal, self.float_theta, self.float_pres)

        for i in range(0, sal.shape[0]):
            for j in range(0, sal.shape[1]):
                if not (np.isnan(sal[i, j]) and np.isnan(self.expected_interp_sal[i, j])):
                    sal_diff = np.abs(sal[i, j] - self.expected_interp_sal[i, j])
                    pres_diff = np.abs(pres[i, j] - self.expected_interp_pres[i, j])

                    # if the difference is less than 10^-12 then consider them equal
                    if sal_diff < 1e-12 and pres_diff < 1e-12:
                        continue

                    self.assertTrue(sal[i, j] == self.expected_interp_sal[i, j],
                                    ("Values at ", i, " and ", j, " do not match for salinity"))
                    self.assertTrue(pres[i, j] == self.expected_interp_pres[i, j],
                                    ("Values at ", i, " and ", j, " do not match for pressure"))


# pylint: disable=too-many-instance-attributes
class MapDataGrid(unittest.TestCase):
    """
    Test cases for map_data_grid function
    """

    def setUp(self):
        self.sal = np.array([34.5, 34.3, 34])
        self.grid = np.array([-59.1, 57.2, 2018, 5108])
        self.data = np.array([[-58, 53.2, 1974, 5224],
                              [-56.5, 63.1, 1987, 4630],
                              [-52.1, 62, 1994, 4657]])
        self.lat = 4
        self.long = 8
        self.age = 20
        self.signal_var = 0.0337
        self.noise_var = 0.0045
        self.phi = 0.5
        self.map_pv_use = 0

    def test_return_types(self):
        """
        Test that the return types are what we expect
        :return: Nothing
        """
        print("Testing that map_data_grid returns the correct types")
        ans = map_data_grid(self.sal, self.grid, self.data,
                            self.lat, self.long, self.age,
                            self.signal_var, self.noise_var, self.phi, self.map_pv_use)

        self.assertTrue(isinstance(ans, tuple), "should return a tuple")
        self.assertTrue(isinstance(ans[0], np.float), "1st return should be a float")
        self.assertTrue(isinstance(ans[1], np.float), "2nd return should a float")
        self.assertTrue(isinstance(ans[2], np.ndarray), "3rd  return should a numpy array")
        self.assertTrue(isinstance(ans[3], np.ndarray), "4th return should be a numpy array")

    def test_return_sizes(self):
        """
        Test that the returned values are the expected size
        :return: Nothing
        """
        print("Testing that map_data_grid returns the expected sizes for each return")

        ans = map_data_grid(self.sal, self.grid, self.data,
                            self.lat, self.long, self.age,
                            self.signal_var, self.noise_var, self.phi, self.map_pv_use)

        expected = self.data.__len__()

        self.assertEqual(ans.__len__(), 4, "should contain 4 answers")
        self.assertEqual(ans[2].__len__(), expected,
                         "Should contain same amount of values as unique points in grid")
        self.assertEqual(ans[3].__len__(), expected,
                         "Should contain same amount of values as unique points in grid")

    def test_return_values(self):
        """
        Test that we get the expected values
        :return: Nothing
        """
        print("Testing that map_data_grid returns the expected values")

        ans = map_data_grid(self.sal, self.grid, self.data,
                            self.lat, self.long, self.age,
                            self.signal_var, self.noise_var, self.phi, self.map_pv_use)
        expected_grid = 34.294578006104572
        expected_grid_error = 0.222578348841383
        expected_data = np.array([34.476629765035661,
                                  34.273698730996550,
                                  34.049671503967787])
        expected_data_error = np.array([0.064804321291622,
                                        0.062091002213780,
                                        0.062089602227937])

        self.assertAlmostEqual(ans[0], expected_grid, 15, "grid mapped field is not as expected")
        self.assertAlmostEqual(ans[1], expected_grid_error, 15, "grid error is not as expected")

        for i in range(0, ans[2].__len__()):
            self.assertAlmostEqual(ans[2][i], expected_data[i], 15,
                                   "grid mapped field is not as expected")
            self.assertAlmostEqual(ans[3][i], expected_data_error[i], 15,
                                   "grid error is not as expected")


if __name__ == '__main__':
    unittest.main()
