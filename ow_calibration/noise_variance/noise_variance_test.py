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
from .noise_variance import noise_variance


class NoiseVarianceTestCase(unittest.TestCase):
    def SetUp(self):
        self.sal = np.array([34.4988, 34.3267, 34.0346])
        self.lat = np.array([-57.9960, -56.4902, -52.0800])
        self.long = np.array([53.1950, 63.1170, 61.9770])


if __name__ == '__main__':
    unittest.main()
