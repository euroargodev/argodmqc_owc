"""
-----Covariance x, y, potential vorticity matrix Test File-----

Written by: Edward Small
When: 12/03/2020

Contains unit tests to check functionality of the `covarxy_pv` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from ow_calibration.build_cov.covarxy_pv.covarxy_pv import covarxy_pv


class MyTestCase(unittest.TestCase):
    """
    Test cases for covarxy_pv function
    """

    def setUp(self):
        """
        Set up for test
        :return: Nothing
        """

        self.input_coords = np.array([0.0572, -0.0592, 5.1083]) * 1.0e+03
        self.coords = np.array([[0.0572, -0.0592, 5.1083],
                                [0.0578, -0.0591, 5.0993],
                                [0.0586, -0.0585, 5.0861]]) * 1.0e+03
        self.long = 4
        self.lat = 2
        self.phi = 0.1
        self.no_pv = 0
        self.yes_pv = 1

    def test_return_shape(self):
        """
        Check that the returned covariance matrix is the correct size
        :return: Nothing
        """
        print("Testing that covarxy_pv returns 1 dimensional matrix")

        cov = covarxy_pv(self.input_coords, self.coords, self.long, self.lat, self.phi, self.no_pv)

        self.assertEqual(cov.shape, (3,), "covarxy_pv matrix is incorrect shape")

    def test_returns_correct(self):
        """
        Check that the returned covriance matrix contains the correct values
        :return: Nothing
        """

        print("Testing that covarxy_pv returns correct values")

        ans = np.array([1, 0.9134, 0.5941])

        cov = covarxy_pv(self.input_coords, self.coords, self.long, self.lat, self.phi, self.no_pv)

        for i in range(0, ans.size):
            self.assertAlmostEqual(ans[i], cov[i],
                                   4, "elements in covarxy_pv matrix are not correct")

    def test_returns_correct_pv(self):
        """
        Check that the returned covriance matrix contains the correct values with
        potential vorticity
        :return: Nothing
        """

        print("Testing that covarxy_pv returns correct values with potential vorticity")

        ans = np.array([1, 0.9101, 0.5828])

        cov = covarxy_pv(self.input_coords, self.coords, self.long, self.lat,
                         self.phi, self.yes_pv)

        for i in range(0, ans.size):
            self.assertAlmostEqual(ans[i], cov[i],
                                   4, "elements in covarxy_pv matrix are not correct")

    def test_returns_ones(self):
        """
        Check that we get a matrix of almost ones if data is very close (according to scale)
        :return: nothing
        """

        print("Testing that covarxy_pv returns almost ones if data is close")

        cov = covarxy_pv(self.input_coords, self.coords, 99999999, 99999999, self.phi, self.no_pv)
        for i in cov:
            self.assertAlmostEqual(i, 1, 15, "elements in covarxy_pv matrix are not correct")


if __name__ == '__main__':
    unittest.main()
