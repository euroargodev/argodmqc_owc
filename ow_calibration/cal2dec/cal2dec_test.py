"""
-----cal2dec Test File-----

Written by: Edward Small
When: 26/09/2019

Contains unit tests to check functionality of the `cal2dec` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
from ow_calibration.cal2dec.cal2dec import cal2dec


class Cal2decTestCase(unittest.TestCase):
    """
    Test cases for cal2dec function
    """

    def test_returns_float(self):
        """
        Check return type is a float
        :return: Nothing
        """
        print("Testing return type is a float...")
        date = cal2dec(0, 1)
        self.assertTrue(isinstance(date, float), "cal2dec should return a float")

    def test_throws_error_if_month_too_large(self):
        """
        Check that an error is thrown if the month exceeds 12
        :return: Nothing
        """
        print("Testing exception is thrown if month is out of scope...")
        with self.assertRaises(Exception) as month_out_of_scope:
            cal2dec(13, 1, 0, 0)

        self.assertTrue('Month is out of scope' in str(month_out_of_scope.exception))

    def test_returns_0_for_first_day(self):
        """
        Check that the first possible day is 0
        :return: Nothing
        """
        print("Testing that day 1 returns 0...")
        self.assertEqual(cal2dec(0, 1, 0, 0), 0, "Should return 0 on first day")

    def test_returns_365_for_last_day(self):
        """
        Check that the last hour in the last day of the last month return 365
        :return: Nothing
        """
        print("Testing that the 24th hour on the last day returns 365")
        self.assertEqual(cal2dec(11, 31, 24, 0), 365, "last value should be 365")

    def test_no_day_larger_than_366(self):
        """
        Check for error if day is larger than 366 (leap year)
        :return: Nothing
        """
        print("Testing that an error is thrown if a date exceeds 365")
        with self.assertRaises(Exception) as date_out_of_scope:
            cal2dec(11, 31, 300, 300)

        self.assertTrue('Day is out of scope of the year' in str(date_out_of_scope.exception))


if __name__ == '__main__':
    unittest.main()
