"""
-----Change Dates Test File-----

Written by: Edward Small
When: 08/11/2019

Contains unit tests to check functionality of the `change_dates` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

from ow_calibration.change_dates.change_dates import change_dates
import unittest
import math
import numpy as np


class MyTestCase(unittest.TestCase):
    """
    Test cases for change_dates function
    """
    def setUp(self):
        """
        Set up variables for testing
        :return: Nothing
        """

        self.good_date1 = 19740201223748
        self.good_date2 = 19740202171947
        self.good_date3 = 19740326211012
        self.bad_date1 = 0
        self.bad_date2 = 12341234
        self.bad_date3 = -5

    def test_returns_array(self):
        """
        Check that the function returns an array
        :return: Nothing
        """
        print("Testing that change_dates returns an array")

        date = change_dates([self.good_date1, self.good_date2, self.good_date3])

        self.assertTrue(isinstance(date, np.ndarray), "return type is not array")

    def test_returns_correct_array(self):
        """
        Check that the function reutnrs the expected values
        :return: Nothing
        """
        print("Testing that change_dates returns expected values")

        expected = [1974.087513318113, 1974.089648021309, 1974.232553272451]

        date = change_dates([self.good_date1, self.good_date2, self.good_date3])
        for i in range(0, date.__len__()):
            self.assertEqual(round(date[i], 6), round(expected[i], 6))


if __name__ == '__main__':
    unittest.main()

change_dates([19740201223748, 19740202171947, 19740326211012, 0])