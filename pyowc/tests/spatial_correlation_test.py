"""
-----Spatial Correlation Test File-----

Written by: Edward Small
When: 01/11/2019

Contains unit tests to check the functionality of the `spatial_correlation` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python

"""

import unittest
import numpy as np
from pyowc.core.stats import spatial_correlation


class MyTestCase(unittest.TestCase):
    """
    Test cases for "spatial_correlation" function
    """

    def test_spatial_correlation_returns_float(self):
        """
        Check that spatial_correlation returns a float if given a float
        :return: Nothing
        """
        print("Testing that spatial_correlation returns a float")

        sc_result = spatial_correlation(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)

        self.assertTrue(isinstance(sc_result, float), "spatial correlation is not a float")

    #pylint: disable=too-many-locals
    def test_spatial_correlation_returns_array(self):
        """
        Check that spatial_correlation returns an array if given an array
        :return: Nothing
        """
        print("Testing that spatial_correlation returns an array")

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

    #pylint: disable=too-many-locals
    def test_spatial_correlation_returns_expected_float(self):
        """
        Check that spatial_correlation returns the correct answer if given a float
        with or without potential vorticity
        :return: Nothing
        """
        print("Testing that spatial_correlation returns correct float value")

        # Without potential vorticity
        hist_long = 54
        float_long = 55
        longitude_large = 4
        hist_lat = -58
        float_lat = -59
        latitude_large = 2
        hist_dates = 2 * 10 ** 3
        float_date = 2.01 * 10 ** 3
        age_large = 10
        phi_large = 0.1
        sc_result = spatial_correlation(hist_long, float_long, longitude_large,
                                        hist_lat, float_lat, latitude_large,
                                        hist_dates, float_date, age_large,
                                        phi_large)

        self.assertEqual(round(sc_result, 4), 1.3125)

        # With potential vorticity
        pv_hist = -0.3 * 10 ** -7
        pv_float = -0.35 * 10 ** -7
        sc_pv_result = spatial_correlation(hist_long, float_long, longitude_large,
                                           hist_lat, float_lat, latitude_large,
                                           hist_dates, float_date, age_large,
                                           phi_large, pv_hist, pv_float)

        self.assertEqual(round(sc_pv_result, 4), 2.489)

    #pylint: disable=too-many-locals
    def test_spatial_correlation_returns_expected_array(self):
        """
                Check that spatial_correlation returns an array if given an array
                :return: Nothing
                """
        print("Testing that spatial_correlation returns an array")

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
        pv_hist = [-0.2354*10**-7, -0.2325*10**-7, -0.267*10**-7]
        pv_float = -2.452 * 10 ** -8
        spatial_correlation_vec = np.vectorize(spatial_correlation)
        sc_expected = [5.158, 13.4935, 14.1804]
        sc_result = spatial_correlation_vec(hist_long, float_long, longitude_large,
                                            hist_lat, float_lat, latitude_large,
                                            hist_dates, float_date, age_large,
                                            phi_large, pv_hist, pv_float)

        for i in range(0, sc_expected.__len__()):
            self.assertEqual(
                round(sc_result[i], 4), sc_expected[i], "spatial_correlation arrays aren't equal"
            )


if __name__ == '__main__':
    unittest.main()
