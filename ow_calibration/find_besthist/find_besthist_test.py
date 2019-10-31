"""
-----Load Configuration Test File-----

Written by: Edward Small
When: 31/10/2019

Contains unit tests to check the functionality of the `find_besthist` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python

"""

import unittest
import numpy as np
from .find_besthist import barotropic_potential_vorticity, spatial_correlation


class FindBestHistTestCase(unittest.TestCase):
    """
    Test cases for find_besthist function, and associated functions
    """

    def test_bpv_returns_float(self):
        """
        Checks that barotropic_potential_vorticity returns a float if
        inputs are just floats.
        :return: Nothing
        """
        print("Testing that barotropic_potential_vorticity returns float...")

        pv_result = barotropic_potential_vorticity(1, 1)
        self.assertTrue(isinstance(pv_result, float), "potential vorticity is not float")

    def test_bpv_vectorised_returns_array(self):
        """
        Check that barotropic_potential_vorticity returns an array if the function
        is vectorised and given lists as inputs
        :return: Nothing
        """
        print("Testing that barotropic_potential_vorticity returns array...")

        lat = [1, 2, 3, 4]
        z_value = [4, 3, 2, 1]
        barotropic_potential_vorticity_vec = np.vectorize(barotropic_potential_vorticity)
        pv_result = barotropic_potential_vorticity_vec(lat, z_value)
        self.assertTrue(isinstance(pv_result, np.ndarray), "potential vorticity vec is not list")

    def test_bpv_never_returns_0(self):
        """
        Check that barotropic_potential_vorticity doesn't return 0 if the equation equals
        0. Should return a very small number.
        :return: Nothing
        """
        print("Testing that barotropic_potential_vorticity returns non-zero...")

        pv_result = barotropic_potential_vorticity(0, 1)
        self.assertNotEqual(pv_result, 0, "potential vorticity should never equal 0")
        self.assertLess(pv_result, 0.001, "potential vorticity should be very small")

    def test_bpv_returns_expected_float(self):
        """
        Check that barotropic_potential_vorticity returns the expected value after calculation
        :return: Nothing
        """
        print("Testing that barotropic_potential_vorticity returns correct value")

        lat = -59.1868
        z_value = 5108.3
        pv_result = barotropic_potential_vorticity(lat, z_value)
        self.assertAlmostEqual(pv_result * 1000, -2.452 * 10 ** -8 * 1000)

    def test_bpv_returns_expected_array(self):
        """
        Check that barotropic_potential_vorticity returns the expected array after calculation
        :return: Nothing
        """
        print("Testing that barotropic_potential_vorticity returns correct array")

        lat = [-57.996, -56.375, -54.496]
        z_value = [5230, 5223.1, 4447.1]
        barotropic_potential_vorticity_vec = np.vectorize(barotropic_potential_vorticity)

        pv_result = barotropic_potential_vorticity_vec(lat, z_value)
        pv_expected = [-0.2365 * 10 ** -7, -0.2325 * 10 ** -7, -0.267 * 10 ** -7]

        for i in range(0, lat.__len__()):
            self.assertAlmostEqual(pv_result[i] * 1000, pv_expected[i] * 1000)

    def test_spatial_correlation_returns_float(self):
        """
        Check that spatial_correlation returns a float if given a float
        :return: Nothing
        """
        print("Testing that spatial_correlation returns a float")

        sc_result = spatial_correlation(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)

        self.assertTrue(isinstance(sc_result, float), "spatial correlation is not a float")

    def test_spatial_correlation_returns_array(self):
        """
        Check that spatial_correlation returns an array if given an array
        :return: Nothing
        """
        print("Testing that spatial_correlation returns an array")

        # pylint: disable=too-many-locals
        hist_long = [53.195, 51.954, 53.107]
        float_long = 57.1794
        longitude_large = 8
        hist_lat = [-57.996, -56.375, -54.496]
        float_lat = -59.1868
        latitude_large = 4
        hist_dates = [1.9741 * 10 ** 3, 1.9471 * 10 ** 3, 1.9472 * 10 ** 3]
        float_date = 2.018 * 10 ** 3
        age_large = 20
        phi_large = 0.5
        pv_hist = [-0.0236 * 10 ** -6, -0.0233 * 10 ** -6, -0.0267 * 10 ** -6]
        pv_float = -2.452
        spatial_correlation_vec = np.vectorize(spatial_correlation)
        sc_result = spatial_correlation_vec(hist_long, float_long, longitude_large,
                                            hist_lat, float_lat, latitude_large,
                                            hist_dates, float_date, age_large,
                                            phi_large, pv_hist, pv_float)

        self.assertTrue(isinstance(sc_result, np.ndarray), "spatial correlation is not array")


if __name__ == '__main__':
    unittest.main()
