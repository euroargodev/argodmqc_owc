"""
-----Change Dates Test File-----

Written by: Edward Small
When: 08/11/2019

Contains unit tests to check functionality of the `change_dates` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from ow_calibration.change_dates.change_dates import change_dates


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
        Check that the function returns the expected values
        :return: Nothing
        """
        print("Testing that change_dates returns expected values")

        expected = [1974.087513318113, 1974.089648021309, 1974.232553272451]

        date = change_dates([self.good_date1, self.good_date2, self.good_date3])
        for i in range(0, date.__len__()):
            self.assertEqual(round(date[i], 6), round(expected[i], 6))

    def test_returns_date_for_missing_minutes(self):
        """
        Check that, even if the calendar date has no minutes, we can
        still get a decimalised date
        :return: Nothing
        """
        print("Testing that change_dates returns date if minutes missing")

        not_expected = [0, 0]

        date = change_dates([197402012237, 197402021719])

        for i in range(0, date.__len__()):
            self.assertNotEqual(round(date[i], 6), round(not_expected[i], 6))

    def test_returns_date_for_missing_hours(self):
        """
        Check that, even if the calendar date has no hours, we can
        still get a decimalised date
        :return: Nothing
        """
        print("Testing that change_dates returns date if hours missing")

        not_expected = [0, 0]

        date = change_dates([1974020122, 1974020217])

        for i in range(0, date.__len__()):
            self.assertNotEqual(round(date[i], 6), round(not_expected[i], 6))

    def test_returns_zeroes_for_negative_date(self):
        """
        Check that giving n negative dates will force the function to
        return array of 0's of length n
        :return: Nothing
        """
        print("Testing that change_dates returns 0's for negatives")

        expected = [0, 0, 0]

        date = change_dates([self.bad_date3, -7, -10])

        for i in range(0, date.__len__()):
            self.assertEqual(round(date[i], 6), round(expected[i], 6))

    def test_returns_zeroes_for_short_date(self):
        """
        Check that giving dates that are too short returns 0s
        :return: Nothing
        """
        print("Testing that change_dates returns 0's for short dates")

        expected = [0, 0, 0]

        date = change_dates([444, 123456, 101])

        for i in range(0, date.__len__()):
            self.assertEqual(round(date[i], 6), round(expected[i], 6))

    def test_returns_zeroes_for_character_date(self):
        """
        Check that giving dates that are not numbers reutnrs 0's
        :return: Nothing
        """
        print("Testing that change_dates returns 0's for character dates")

        expected = [0, 0, 0]

        date = change_dates(["<very>", "{bad}", "(date) 'given'"])

        for i in range(0, date.__len__()):
            self.assertEqual(round(date[i], 6), round(expected[i], 6))


if __name__ == '__main__':
    unittest.main()
