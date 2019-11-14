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
        """
        set up variables to use for testing
        :return: Nothing
        """
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
        """
        Check that the function returns an array if given an array
        :return: Nothing
        """
        print("Testing that covar_xyt_pv returns an array")
        covar = covar_xyt_pv(self.points1, self.points2, self.lat, self.long,
                             self.age, self.phi, self.p_v)

        self.assertTrue(isinstance(covar, np.ndarray), "Covariance isn't numpy array")

    def test_throws_error(self):
        """
        Check that the function throws an exception if "pointsx" has no length
        :return: Nothing
        """
        print("Testing that covar_xyt_pv throws exception if points is not a matrix")

        with self.assertRaises(Exception) as no_length:
            covar_xyt_pv(1, "hello", self.lat, self.long,
                         self.age, self.phi, self.p_v)

        self.assertTrue("A set of points has no length associated with it"
                        in str(no_length.exception))

    def test_returns_ones_for_same_value(self):
        """
        Check that entering the same value for points 1 and 2 gives a matrix of ones
        :return: Nothing
        """
        print("Testing that we return 1's if input points are identical")

        covar = covar_xyt_pv(np.array([[1, 2, 3, 4], [1, 2, 3, 4]]),
                             np.array([[1, 2, 3, 4], [1, 2, 3, 4]]),
                             self.lat, self.long, self.age, self.phi, self.p_v)
        expected = np.array([[1, 1], [1, 1]])

        for i in range(0, covar.__len__()):
            for j in range(0, covar[i].__len__()):
                self.assertEqual(covar[i][j], expected[i][j], "covariance isn't just 1's")

    def test_returns_zeroes_for_extremely_different_value(self):
        """
        Check that entering extreme value for points 1 and 2 gives a matrix of zeroes
        :return: Nothing
        """
        print("Testing that we return 0's if input points are extremely different")

        covar = covar_xyt_pv(np.array([[1, 2, 3, 4], [1000, 2000, 3000, 4000]]),
                             np.array([[-1000, -2000, -3000, 4000],
                                       [-1*9999, 2*9999, -3*9999, -4*999]]),
                             self.lat, self.long, self.age, self.phi, self.p_v)
        expected = np.array([[0, 0], [0, 0]])

        for i in range(0, covar.__len__()):
            for j in range(0, covar[i].__len__()):
                self.assertAlmostEqual(covar[i][j]*1000000, expected[i][j]*1000000,
                                       "covariance isn't just 1's")

    def test_return_matrix_shape_correct(self):
        """
        Check that, if given an m*4 points 1 matrix and a n n*5 points 2 matrix
        the function returns an m*n matrix
        :return: Nothing
        """
        print("Testing that the returned covariance matrix is the correct shape")

        covar = covar_xyt_pv(np.array([[1, 2, 3, 4], [1, 2, 3, 4]]),
                             np.array([[1, 2, 3, 4], [1, 2, 3, 4]]),
                             self.lat, self.long, self.age, self.phi, self.p_v)

        expected = (2, 2)

        for i in range(0, covar.shape.__len__()):
            self.assertEqual(covar.shape[i], expected[i], "covariance matrix is wrong shape")


if __name__ == '__main__':
    unittest.main()
