"""
-----Build Covariance Matrix Test File-----

Written by: Edward Small
When: 12/03/2020

Contains unit tests to check functionality of the `build_cov` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
import scipy.io as scipy
from ow_calibration.build_cov.build_cov import build_cov


class MyTestCase(unittest.TestCase):
    """
    Test cases for build_cov function
    """

    def setUp(self):
        """
        Set up for test
        :return: Nothing
        """

        self.config = {"MAPSCALE_LONGITUDE_SMALL": 4,
                       "MAPSCALE_LATITUDE_SMALL": 2,
                       "MAPSCALE_PHI_SMALL": 0.1,
                       "MAP_USE_PV": 0}
        self.ptmp = np.array([[0.7058, 0.7039, 0.8285],
                              [0.6713, 0.6664, 0.7432],
                              [0.8257, 0.8224, 0.7804],
                              [0.7452, 0.7411, 1.1980],
                              [0.7836, 0.7802, 1.1504],
                              [1.2008, 1.2010, 1.2497],
                              [1.1496, 1.1481, 1.3036],
                              [1.2520, 1.2553, 1.0921],
                              [1.3039, 1.3046, np.nan],
                              [1.0947, 1.0962, np.nan]])
        self.coord_float = np.array([[0.0572, -0.0592, 5.1083],
                                     [0.0578, -0.0591, 5.0993],
                                     [0.0586, -0.0585, 5.0861]]) * 1.0e+03

    def test_returns_numpy_array(self):
        """
        Check that build_cov returns a numpy array
        :return: Nothing
        """

        print("Testing that build_cov returns a numpy array")

        test = build_cov(self.ptmp, self.coord_float, self.config)

        self.assertEqual(type(test), np.ndarray, "build_cov did not return a numpy array")

    def test_returns_correct_size(self):
        """
        Check that build_cov returns a matrix that is the correct size
        :return: Nothing
        """

        print("Testing that build_cov returns correct shape matrix")

        test = build_cov(self.ptmp, self.coord_float, self.config)

        self.assertEqual(test.shape,
                         (self.coord_float.shape[0] * self.ptmp.shape[0],
                          self.coord_float.shape[0] * self.ptmp.shape[0]),
                         "build_cov returned a matrix of incorrect size")

    def test_returns_correct_elements(self):
        """
        Check that build_cov returns a matrix that is the correct size
        :return: Nothing
        """

        print("Testing that build_cov returns correct shape matrix")

        expected = scipy.loadmat("data/test_data/build_cov/cov.mat")['test_cov_1']
        expected_size = expected.shape

        test = build_cov(self.ptmp, self.coord_float, self.config)

        for i in range(0, expected_size[0]):
            for j in range(0, expected_size[1]):
                self.assertAlmostEqual(test[i, j], expected[i, j], 15,
                                       "covariance matrix is incorrect")


if __name__ == '__main__':
    unittest.main()
