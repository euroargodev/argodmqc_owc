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
from pyowc import core

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
        covar = core.stats.covar_xyt_pv(self.points1, self.points2, self.lat, self.long,
                             self.age, self.phi, self.p_v)

        self.assertTrue(isinstance(covar, np.ndarray), "Covariance isn't numpy array")

    def test_returns_ones_for_same_value(self):
        """
        Check that entering the same value for points 1 and 2 gives a matrix of ones
        :return: Nothing
        """
        print("Testing that covar_xyt_pv return 1's if input points are identical")

        covar = core.stats.covar_xyt_pv(np.array([[1, 2, 3, 4], [1, 2, 3, 4]]),
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
        print("Testing that covar_xyt_pv return 0's if input points are extremely different")

        covar = core.stats.covar_xyt_pv(np.array([[1, 2, 3, 4], [1000, 2000, 3000, 4000]]),
                             np.array([[-1000, -2000, -3000, 4000],
                                       [-1 * 9999, 2 * 9999, -3 * 9999, -4 * 999]]),
                             self.lat, self.long, self.age, self.phi, self.p_v)
        expected = np.array([[0, 0], [0, 0]])

        for i in range(0, covar.__len__()):
            for j in range(0, covar[i].__len__()):
                self.assertAlmostEqual(covar[i][j] * 1000000, expected[i][j] * 1000000,
                                       "covariance isn't just 1's")

    def test_return_matrix_shape_correct_multidimensional(self):
        """
        Check that, if given an m*4 points 1 matrix and a n n*5 points 2 matrix
        the function returns an m*n matrix (m > 1)
        :return: Nothing
        """
        print("Testing that covar_xyt_pv returns matrix of the correct shape for multidemensional")

        covar = core.stats.covar_xyt_pv(self.points1, self.points2, self.lat, self.long,
                             self.age, self.phi, self.p_v)

        expected = (3, 3)

        for i in range(0, covar.shape.__len__()):
            self.assertEqual(covar.shape[i], expected[i], "covariance matrix is wrong shape")

    def test_return_matrix_shape_correct_one_demensional(self):
        """
        Check that, if given an 1*4 points 1 matrix and a n n*5 points 2 matrix
        the function returns an 1*n matrix (m > 1)
        :return: Nothing
        """
        print("Testing that covar_xyt_pv returns matrix of the correct shape one dimensional")

        covar = core.stats.covar_xyt_pv(np.array([-0.057996, 0.053195, 1.9740875, 5.229838]) * 10 ** 3,
                             self.points2, self.lat, self.long, self.age, self.phi, self.p_v)
        expected_shape = (1, 3)

        self.assertEqual(covar.shape, expected_shape, "covariance matrix is wrong shape")

    def test_returns_expected_answers(self):
        """
        Check that we get the answers we expect
        :return: Nothing
        """
        print("Testing that covar_xyt_pv returns the expected result")

        expected = np.array([[1, 0.001350185414046, 0.001712995003897],
                             [0.001350185414046, 1, 0.600332430927527],
                             [0.001712995003897, 0.600332430927527, 1]])
        covar = core.stats.covar_xyt_pv(self.points1, self.points2, self.lat, self.long,
                             self.age, self.phi, self.p_v)

        for i in range(0, covar.__len__()):
            for j in range(0, covar[i].__len__()):
                self.assertAlmostEqual(covar[i][j], expected[i][j], 15,
                                       "covariances is not as expected")

    def test_allows_1_dimensional_data(self):
        """
        Check that using 1 dimensional data sets does not through an error
        :return: Nothing
        """
        print("Testing that covar_xyt_pv can use one dimensional data")

        covar = core.stats.covar_xyt_pv(np.array([-0.057996, 0.053195, 1.9740875, 5.229838]) * 10 ** 3,
                             self.points2, self.lat, self.long, self.age, self.phi, self.p_v)

        expected_ans = [1, 0.122561894349782, 0.012339774448825]

        for i in range(covar.__len__()):
            self.assertAlmostEqual(covar[0][i], expected_ans[i], 15,
                                   "unexpected answer with one dimensional data")


if __name__ == '__main__':
    unittest.main()
