"""
-----Find Besthist Test File-----

Written by: Edward Small
When: 31/10/2019

Contains unit tests to check the functionality of the `find_besthist` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python

"""

import unittest
import random
import numpy as np
from ow_calibration.find_besthist.find_besthist import find_besthist


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
        random.seed(500)
        self.grid_lat = np.random.rand(1000) * 5 * random.choice([-1, 1]) + -59.1868
        self.grid_long = np.random.rand(1000) * 5 * random.choice([-1, 1]) + 57.1794
        self.grid_z_values = np.random.rand(1000) * 5 * random.choice([-1, 1]) + 2.018 * 10 ** 3
        self.grid_dates = np.random.rand(1000) * 5 * random.choice([-1, 1]) + 5.1083 * 10 ** 3
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

    def test_index_returns_array(self):
        """
        Check that find_besthist gives back an array
        :return: Nothing
        """
        print("Testing that find_besthist gives back a numpy array")

        index = find_besthist(self.grid_lat, self.grid_long, self.grid_z_values, self.grid_dates,
                              self.lat, self.long, self.date, self.z_value,
                              self.lat_large, self.lat_small, self.long_large, self.long_small,
                              self.phi_large, self.phi_small, self.age_large, self.age_small,
                              self.map_pv_usage, self.max_casts)

        self.assertTrue(isinstance(index, np.ndarray), "find_besthist not returning array")

    def test_index_returns_empty(self):
        """
        Check that find_besthist gives back an empty array if no data fits inside
        the ellipse
        :return: Nothing
        """
        print("Testing that find_besthist gives back a numpy array")

        grid_lat = np.array([1, 1])
        grid_long = np.array([10, 10])
        grid_z_values = np.array([10, 10])
        grid_dates = np.array([0, 0])

        index = find_besthist(grid_lat, grid_long, grid_dates, grid_z_values,
                              self.lat, self.long, self.date, self.z_value,
                              self.lat_large, self.lat_small, self.long_large, self.long_small,
                              self.phi_large, self.phi_small, self.age_large, self.age_small,
                              self.map_pv_usage, self.max_casts)

        self.assertTrue(index.__len__() == 0, "No data should have been selected")


if __name__ == '__main__':
    unittest.main()
