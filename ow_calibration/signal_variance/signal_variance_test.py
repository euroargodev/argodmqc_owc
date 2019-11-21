"""
-----LSignal Variance Test File-----

Written by: Edward Small
When: 21/11/2019

Contains unit tests to check the functionality of the `signal_variance` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""


import unittest
from .signal_variance import signal_variance


class SignalVarianceTestCase(unittest.TestCase):
    """
    Test cases for signal_variance function
    """

    def test_returns_float(self):
        """
        Check that we return a float if given some data
        :return: Nothing
        """
        print("Testing that signal_variance returns a float")

        var = signal_variance([1, 2, 3, 4, 5])
        self.assertTrue(isinstance(var, float), "signal variance is not a float")


if __name__ == '__main__':
    unittest.main()
