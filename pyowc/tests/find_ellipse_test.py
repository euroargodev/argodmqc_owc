"""
-----Find Ellipse Test File-----

Written by: Edward Small
When: 07/11/2019

Contains unit tests to check the functionality of the `find_ellipse` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python

"""

import unittest
import numpy as np
from pyowc import core

class MyTestCase(unittest.TestCase):
    """
    Test cases for "find ellipse" function
    """

    def test_find_ellipse_returns_float(self):
        """
        Check that find_ellipse returns a float if only given floats
        :return: Nothing
        """
        print("Testing that find_ellipse will return a float")

        ellipse = core.finders.find_ellipse(1, 1, 1, 1, 1, 1, 1, 1)

        self.assertTrue(isinstance(ellipse, float), "find ellipse did not return float")

    def test_find_ellipse_returns_array(self):
        """
        Check that find_ellipse returns an array if vectorised
        :return: Nothing
        """
        print("Testing that find_ellipse vectorised returns an array")

        find_ellipse_vec = np.vectorize(core.finders.find_ellipse)
        ellipse = find_ellipse_vec([1, 1], 1, 1, [1, 1], 1, 1, 1)

        self.assertTrue(isinstance(ellipse, np.ndarray), "find_ellipse did not return array")

    def test_find_ellipse_returns_expected_float(self):
        """
        Check that find_ellipse returns the correct float
        :return: Nothing
        """
        print("Testing that find_ellipse returns the right answer")

        ellipse = core.finders.find_ellipse(53.195, 57.1794, 8,
                               -57.996, -59.1868, 4,
                               0.5, -2.3547 * 10 ** -8, -2.452 * 10 ** -8)

        self.assertEqual(round(ellipse, 4), 0.2017, "find_ellipse returned incorrect float")

    def test_find_ellipse_returns_expected_float_without_pv(self):
        """
        Check that find_ellipse returns the correct float without potential vorticity
        :return: Nothing
        """
        print("Testing that find_ellipse returns the right answer without potential vorticity")

        ellipse = core.finders.find_ellipse(53.195, 57.1794, 8,
                               -57.996, -59.1868, 4,
                               0.5)

        self.assertEqual(round(ellipse, 4), 0.1934, "find_ellipse returned incorrect float "
                                                    "without potential vortcity")

    def test_find_ellipse_returns_expected_array(self):
        """
        Check that find_ellipse returns the correct array
        :return: Nothing
        """
        print("Testing that find_ellipse returns the right array answer")

        ellipse_expected = [0.2018, 0.3286, 0.4428]
        find_ellipse_vec = np.vectorize(core.finders.find_ellipse)
        hist_long = [53.195, 51.954, 53.107]
        hist_lat = [-57.996, -56.375, -54.496]
        hist_pv = [-0.2354 * 10 ** -7,
                   -0.2325 * 10 ** -7,
                   -0.267 * 10 ** -7]
        ellipse = find_ellipse_vec(hist_long, 57.1794, 8,
                                   hist_lat, -59.1868, 4,
                                   0.5, hist_pv, -2.452 * 10 ** -8)

        for i in range(0, hist_long.__len__()):
            self.assertEqual(round(ellipse[i], 4), ellipse_expected[i],
                             "incorrect result in array")

    def test_find_ellipse_returns_expected_array_without_pv(self):
        """
        Check that find_ellipse returns the correct array
        :return: Nothing
        """
        print("Testing that find_ellipse returns the right array answer")

        ellipse_expected = [0.1934, 0.3199, 0.4261]
        find_ellipse_vec = np.vectorize(core.finders.find_ellipse)
        hist_long = [53.195, 51.954, 53.107]
        hist_lat = [-57.996, -56.375, -54.496]
        ellipse = find_ellipse_vec(hist_long, 57.1794, 8,
                                   hist_lat, -59.1868, 4,
                                   0.5)

        for i in range(0, hist_long.__len__()):
            self.assertEqual(round(ellipse[i], 4), ellipse_expected[i],
                             "incorrect result in array")


if __name__ == '__main__':
    unittest.main()
