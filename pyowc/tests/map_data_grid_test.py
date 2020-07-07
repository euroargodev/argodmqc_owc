"""
-----Map Data to Grid Test File-----

Written by: Edward Small
When: 24/09/2019

Contains unit tests to check the functionality of the `map_data_grid` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from pyowc.data.wrangling import map_data_grid

#pylint: disable=too-many-instance-attributes
class MapDataGridTestCase(unittest.TestCase):
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
