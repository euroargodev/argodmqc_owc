"""
-----Noise Variance Test File-----

Written by: Edward Small
When: 24/09/2019

Contains unit tests to check the functionality of the `noise_variance` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from pyowc.core.stats import noise_variance

class NoiseVarianceTestCase(unittest.TestCase):
    """
    Test cases for noise_variance function
    """

    def setUp(self):
        self.sal = np.array([34.4988, 34.3267, 34.0346])
        self.lat = np.array([-57.9960, -56.4902, -52.0800])
        self.long = np.array([53.1950, 63.1170, 61.9770])

    def test_returns_float(self):
        """
        Check that noise_variance returns a float
        :return: Nothing
        """
        print("Testing that noise_variance returns a float")

        noise_var1 = noise_variance(self.sal, self.lat, self.long)
        noise_var2 = noise_variance(np.array([0, 2]), np.array([1, -1]), np.array([-1, 1]))

        self.assertTrue(isinstance(noise_var1, float), "noise variance is not a float")
        self.assertTrue(isinstance(noise_var2, float), "noise variance is not a float")

    def test_returns_0_if_no_unique_points(self):
        """
        Check that noise_variance returns 0 if it cannot find any unique points
        :return: Nothing
        """
        print("Testing that noise_variance returns 0 for no unique points")

        noise_var = noise_variance(np.array([0, 0]), np.array([-1, -1]), np.array([1, 1]))
        self.assertEqual(noise_var, 0, "Variance is not 0 for equal points")

    def test_returns_expected(self):
        """
        Check that noise_variance returns the expected answer
        :return: Nothing
        """
        print("Testing that noise_variance returns the expected value")

        expected = 0.033377205000001
        noise_var = noise_variance(self.sal, self.lat, self.long)

        self.assertAlmostEqual(noise_var, expected, 16, "Did not receive expected answer")


if __name__ == '__main__':
    unittest.main()
