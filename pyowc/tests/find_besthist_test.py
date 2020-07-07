"""
-----Find Best Hist Test File-----

Written by: Edward Small
When: 31/10/2019

Contains unit tests to check the functionality of the `find_besthist` function
Difficult to test fully because the function selects 1/3 of its points randomly

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python

"""

import unittest
import random
import numpy as np
# from .find_besthist import find_besthist
from pyowc import core

class FindBestHistTestCase(unittest.TestCase):
    """
    Test cases for find_besthist function
    """

    # pylint: disable=too-many-instance-attributes
    def setUp(self):
        """
        Set up for the tests. Creates 1000 historical random data points that
        are the float data points +/- 5. This could fail randomly if, by pure
        chance, none of the generated data is in the ellipse, so use pseudo random
        numbers just in case
        :return: nothing
        """
        random.seed(1)
        self.grid_lat = np.random.rand(1000) * 5 * random.choice([-1, 1]) + -59.1868
        self.grid_long = np.random.rand(1000) * 5 * random.choice([-1, 1]) + 57.1794
        self.grid_dates = np.random.rand(1000) * 5 * random.choice([-1, 1]) + 5.1083 * 10 ** 3
        self.grid_z_values = np.random.rand(1000) * 5 * random.choice([-1, 1]) + 2.018 * 10 ** 3
        self.lat = -59.1868
        self.long = 57.1794
        self.z_value = 2.018 * 10 ** 3
        self.date = 5.1083 * 10 ** 3
        self.lat_large = 4
        self.lat_small = 2
        self.long_large = 8
        self.long_small = 4
        self.phi_large = 0.5
        self.phi_small = 0.1
        self.age_large = 20
        self.age_small = 10
        self.map_pv_usage = 1
        self.max_casts = 300

    def test_returns_array(self):
        """
        Check that find_besthist gives back an array
        :return: Nothing
        """
        print("Testing that find_besthist gives back a numpy array")

        index = core.finders.find_besthist(self.grid_lat, self.grid_long, self.grid_dates, self.grid_z_values,
                              self.lat, self.long, self.date, self.z_value,
                              self.lat_large, self.lat_small, self.long_large, self.long_small,
                              self.phi_large, self.phi_small, self.age_large, self.age_small,
                              self.map_pv_usage, self.max_casts)

        self.assertTrue(isinstance(index, np.ndarray), "find_besthist not returning array")

    def test_returns_empty(self):
        """
        Check that find_besthist gives back an empty array if no data fits inside
        the ellipse
        :return: Nothing
        """
        print("Testing that find_besthist gives back an empty array if no data inside ellipse")

        grid_lat = np.array([1, 1])
        grid_long = np.array([10, 10])
        grid_z_values = np.array([10, 10])
        grid_dates = np.array([0, 0])

        index = core.finders.find_besthist(grid_lat, grid_long, grid_dates, grid_z_values,
                              self.lat, self.long, self.date, self.z_value,
                              self.lat_large, self.lat_small, self.long_large, self.long_small,
                              self.phi_large, self.phi_small, self.age_large, self.age_small,
                              self.map_pv_usage, self.max_casts)

        self.assertTrue(index.__len__() == 0, "No data should have been selected")

    def test_returns_array_of_right_size(self):
        """
        Check that find_besthist gives back an array that is the size we asked
        :return: Nothing
        """
        print("Testing that find_besthist gives back the right sized array")

        index = core.finders.find_besthist(self.grid_lat, self.grid_long, self.grid_dates, self.grid_z_values,
                              self.lat, self.long, self.date, self.z_value,
                              self.lat_large, self.lat_small, self.long_large, self.long_small,
                              self.phi_large, self.phi_small, self.age_large, self.age_small,
                              self.map_pv_usage, self.max_casts)

        self.assertTrue(index.__len__() == self.max_casts, "index is incorrect size")

    def test_returns_expected_values(self):
        """
        Check that find_besthist gives us the values we expect
        Since 1/3 of the data could be selected randomly, we will just check
        that it removes data that isn't inside the ellipse
        :return: Nothing
        """
        print("Testing that find_besthist gives back an array containing expected values")
        grid_lat = np.array([1, 1, -60, -58])
        grid_long = np.array([10, 10, 59, 57])
        grid_dates = np.array([1, 1, 5.108 * 10 ** 3, 5.109 * 10 ** 3])
        grid_z_values = np.array([1, 1, 2.02 * 10 ** 3, 2.015 * 10 ** 3])
        expected = np.array([2, 3])

        index = core.finders.find_besthist(grid_lat, grid_long, grid_dates, grid_z_values,
                              self.lat, self.long, self.date, self.z_value,
                              self.lat_large, self.lat_small, self.long_large, self.long_small,
                              self.phi_large, self.phi_small, self.age_large, self.age_small,
                              self.map_pv_usage, self.max_casts)

        for i in range(0, index.__len__()):
            self.assertTrue(index[i] == expected[i], "output is incorrect (wrong index selected)")

if __name__ == '__main__':
    unittest.main()
