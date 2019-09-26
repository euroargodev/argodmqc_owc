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
        date = cal2dec(0, 0, 0, 0)
        self.assertEqual(type(date), float, "return type isn't float")


if __name__ == '__main__':
    unittest.main()
