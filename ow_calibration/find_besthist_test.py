"""
-----Find Besthist Test File-----

Written by: Edward Small
When: 31/10/2019

Contains unit tests to check the functionality of the `find_besthist` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python

"""

import unittest
import numpy as np
import random
from ow_calibration.find_besthist.find_besthist import find_besthist


class FindBestHistTestCase(unittest.TestCase):
    """
    Test cases for find_besthist function
    """
    def do_something(self):
        lat_test = np.random.rand(1000) * 10 * random.choice([-1, 1]) + -57.996
        long_test = np.random.rand(1000) * 10 * random.choice([-1, 1]) + 57.1794
        z_test = np.random.rand(1000) * 10 * random.choice([-1, 1]) + 2.018 * 10 ** 3
        age_test = np.random.rand(1000) * 10 * random.choice([-1, 1]) + 5.1083 * 10 ** 3

        hello = find_besthist(lat_test, long_test, z_test, age_test,
                      -59.1868, 57.1794, 2.018 * 10 ** 3, 5.1083 * 10 ** 3,
                      4, 2, 8, 4, 0.5, 0.1, 20, 10, 1, 300)
        print(hello)


if __name__ == '__main__':
    unittest.main()

"""# random.seed(1)
lat_test = np.random.rand(1000) * 10 * random.choice([-1, 1]) + -57.996
long_test = np.random.rand(1000) * 10 * random.choice([-1, 1]) + 57.1794
z_test = np.random.rand(1000) * 10 * random.choice([-1, 1]) + 2.018 * 10 ** 3
age_test = np.random.rand(1000) * 10 * random.choice([-1, 1]) + 5.1083 * 10 ** 3

find_besthist(lat_test, long_test, z_test, age_test,
              -59.1868, 57.1794, 2.018 * 10 ** 3, 5.1083 * 10 ** 3,
              4, 2, 8, 4, 0.5, 0.1, 20, 10, 1, 300)

find_besthist([-57.996, -56.375], [53.195, 51.954],
              [1.9741 * 10 ** 3, 1.9741 * 10 ** 3], [5.23 * 10 ** 3, 5.2231 * 10 ** 3],
              -59.1868, 57.1794, 2.018 * 10 ** 3, 5.1083 * 10 ** 3,
              4, 2, 8, 4, 0.5, 0.1, 10, 20, 0, 300)"""