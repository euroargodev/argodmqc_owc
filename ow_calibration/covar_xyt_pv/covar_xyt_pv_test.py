"""
-----Covariance xyt pv test file Test File-----

Written by: Edward Small
When: 24/09/2019

Contains unit tests to check the functionality of the `covar_xyt_pv` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from .covar_xyt_pv import covar_xyt_pv

class Covarxytpv(unittest.TestCase):
    """
    Test cases for covar_xyt_pv function
    """

    def setUp(self):
        self.points1 = np.array([[-0.057996, 0.053195, 1.9740875, 5.229838],
              [-0.0564902, 0.0631170, 1.9870367, 4.6300392],
              [-0.05208, 0.0619770, 1.9941118, 4.6536932]]) * 10 ** 3
        self.points2 = self.points1
        self.lat = 4
        self.long = 8
        self.age = 20
        self.phi = 0.5
        self.p_v = 0

    def test_returns_array(self):
        covar = covar_xyt_pv(self.points1, self.points2, self.lat, self.long,
                             self.age, self.phi, self.p_v)

        self.assertTrue(isinstance(covar, np.ndarray), "Covariance isn't numpy array")



if __name__ == '__main__':
    unittest.main()
