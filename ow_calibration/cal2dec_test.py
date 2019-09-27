"""
-----cal2dec Test File-----

Written by: Edward Small
When: 26/09/2019

Contains unit tests to check functionality of the `cal2dec` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
from .cal2dec import cal2dec


class MyTestCase(unittest.TestCase):

    def test_returns_float(self):
        date = cal2dec(0, 1)
        self.assertTrue(type(cal2dec(0, 1, 0, 0)) is float, "cal2dec should return a float")

    def test_throws_error_if_month_too_large(self):
        with self.assertRaises(Exception) as month_out_of_scope:
            cal2dec(13, 1, 0, 0)

        self.assertTrue('Month is out of scope' in str(month_out_of_scope.exception))

    def test_returns_0_for_first_day(self):
        self.assertEqual(cal2dec(0, 1, 0, 0), 0, "Should return 0 on first day")

    def test_returns_365_for_last_day(self):
        self.assertEqual(cal2dec(11, 31, 24, 0), 365, "last value should be 365")


if __name__ == '__main__':
    unittest.main()
