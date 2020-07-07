"""
-----Potential Vorticity Test File-----

Written by: Edward Small
When: 07/11/2019

Contains unit tests to check the functionality of the `potential_vorticity` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python

"""

import unittest
import numpy as np
from pyowc.core.utils import potential_vorticity

class PotentialVorticityTestCase(unittest.TestCase):
    """
    Test cases for potential_vorticity function
    """
    def test_bpv_returns_float(self):
        """
        Checks that potential_vorticity returns a float if
        inputs are just floats.
        :return: Nothing
        """
        print("Testing that potential_vorticity returns float...")

        pv_result = potential_vorticity(1, 1)
        self.assertTrue(isinstance(pv_result, float), "potential vorticity is not float")

    def test_bpv_vectorised_returns_array(self):
        """
        Check that potential_vorticity returns an array if the function
        is vectorised and given lists as inputs
        :return: Nothing
        """
        print("Testing that potential_vorticity returns array...")

        lat = [1, 2, 3, 4]
        z_value = [4, 3, 2, 1]
        barotropic_potential_vorticity_vec = np.vectorize(potential_vorticity)
        pv_result = barotropic_potential_vorticity_vec(lat, z_value)
        self.assertTrue(isinstance(pv_result, np.ndarray), "potential vorticity vec is not list")

    def test_bpv_never_returns_0(self):
        """
        Check that potential_vorticity doesn't return 0 if the equation equals
        0. Should return a very small number.
        :return: Nothing
        """
        print("Testing that potential_vorticity returns non-zero...")

        pv_result = potential_vorticity(0, 1)
        self.assertNotEqual(pv_result, 0, "potential vorticity should never equal 0")
        self.assertLess(pv_result, 0.001, "potential vorticity should be very small")

    def test_bpv_returns_expected_float(self):
        """
        Check that potential_vorticity returns the expected value after calculation
        :return: Nothing
        """
        print("Testing that potential_vorticity returns correct value")

        lat = -59.1868
        z_value = 5108.3
        pv_result = potential_vorticity(lat, z_value)
        self.assertAlmostEqual(pv_result * 1000, -2.452 * 10 ** -8 * 1000)

    def test_bpv_returns_expected_array(self):
        """
        Check that potential_vorticity returns the expected array after calculation
        :return: Nothing
        """
        print("Testing that potential_vorticity returns correct array")

        lat = [-57.996, -56.375, -54.496]
        z_value = [5230, 5223.1, 4447.1]
        barotropic_potential_vorticity_vec = np.vectorize(potential_vorticity)

        pv_result = barotropic_potential_vorticity_vec(lat, z_value)
        pv_expected = [-0.2365 * 10 ** -7, -0.2325 * 10 ** -7, -0.267 * 10 ** -7]

        for i in range(0, lat.__len__()):
            self.assertAlmostEqual(pv_result[i] * 1000, pv_expected[i] * 1000)


if __name__ == '__main__':
    unittest.main()
