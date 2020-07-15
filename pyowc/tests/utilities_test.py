
import os
import unittest
import math
import numpy as np

from pyowc import utilities as utils
from . import TESTS_CONFIG


class Cal2dec(unittest.TestCase):
    """
    Test cases for cal2dec function
    """

    def test_returns_float(self):
        """
        Check return type is a float
        :return: Nothing
        """
        print("Testing return type is a float...")
        date = utils.cal2dec(0, 1)
        self.assertTrue(isinstance(date, float), "cal2dec should return a float")

    #PUT RETURN HERE

    def test_throws_error_if_month_too_large(self):
        """
        Check that an error is thrown if the month exceeds 12
        :return: Nothing
        """
        print("Testing exception is thrown if month is out of scope...")
        with self.assertRaises(Exception) as month_out_of_scope:
            utils.cal2dec(13, 1, 0, 0)

        self.assertTrue('Month is out of scope' in str(month_out_of_scope.exception))

    def test_returns_0_for_first_day(self):
        """
        Check that the first possible day is 0
        :return: Nothing
        """
        print("Testing that day 1 returns 0...")
        self.assertEqual(utils.cal2dec(0, 1, 0, 0), 0, "Should return 0 on first day")

    def test_returns_365_for_last_day(self):
        """
        Check that the last hour in the last day of the last month return 365
        :return: Nothing
        """
        print("Testing that the 24th hour on the last day returns 365")
        self.assertEqual(utils.cal2dec(11, 31, 24, 0), 365, "last value should be 365")

    def test_no_day_larger_than_366(self):
        """
        Check for error if day is larger than 366 (leap year)
        :return: Nothing
        """
        print("Testing that an error is thrown if a date exceeds 365")
        with self.assertRaises(Exception) as date_out_of_scope:
            utils.cal2dec(11, 31, 300, 300)

        self.assertTrue('Day is out of scope of the year' in str(date_out_of_scope.exception))


class ChangeDates(unittest.TestCase):
    """
    Test cases for change_dates function
    """

    def setUp(self):
        """
        Set up variables for testing
        :return: Nothing
        """

        self.good_date1 = 19740201223748
        self.good_date2 = 19740202171947
        self.good_date3 = 19740326211012
        self.bad_date1 = 0
        self.bad_date2 = 12341234
        self.bad_date3 = -5

    def test_returns_array(self):
        """
        Check that the function returns an array
        :return: Nothing
        """
        print("Testing that change_dates returns an array")

        date = utils.change_dates([self.good_date1, self.good_date2, self.good_date3])

        self.assertTrue(isinstance(date, np.ndarray), "return type is not array")

    def test_returns_correct_array(self):
        """
        Check that the function returns the expected values
        :return: Nothing
        """
        print("Testing that change_dates returns expected values")

        expected = [1974.087513318113, 1974.089648021309, 1974.232553272451]

        date = utils.change_dates([self.good_date1, self.good_date2, self.good_date3])
        for i in range(0, date.__len__()):
            self.assertEqual(round(date[i], 6), round(expected[i], 6))

    def test_returns_date_for_missing_minutes(self):
        """
        Check that, even if the calendar date has no minutes, we can
        still get a decimalised date
        :return: Nothing
        """
        print("Testing that change_dates returns date if minutes missing")

        not_expected = [0, 0]

        date = utils.change_dates([197402012237, 197402021719])

        for i in range(0, date.__len__()):
            self.assertNotEqual(round(date[i], 6), round(not_expected[i], 6))

    def test_returns_date_for_missing_hours(self):
        """
        Check that, even if the calendar date has no hours, we can
        still get a decimalised date
        :return: Nothing
        """
        print("Testing that change_dates returns date if hours missing")

        not_expected = [0, 0]

        date = utils.change_dates([1974020122, 1974020217])

        for i in range(0, date.__len__()):
            self.assertNotEqual(round(date[i], 6), round(not_expected[i], 6))

    def test_returns_zeroes_for_negative_date(self):
        """
        Check that giving n negative dates will force the function to
        return array of 0's of length n
        :return: Nothing
        """
        print("Testing that change_dates returns 0's for negatives")

        expected = [0, 0, 0]

        date = utils.change_dates([self.bad_date3, -7, -10])

        for i in range(0, date.__len__()):
            self.assertEqual(round(date[i], 6), round(expected[i], 6))

    def test_returns_zeroes_for_short_date(self):
        """
        Check that giving dates that are too short returns 0s
        :return: Nothing
        """
        print("Testing that change_dates returns 0's for short dates")

        expected = [0, 0, 0]

        date = utils.change_dates([444, 123456, 101])

        for i in range(0, date.__len__()):
            self.assertEqual(round(date[i], 6), round(expected[i], 6))

    def test_returns_zeroes_for_character_date(self):
        """
        Check that giving dates that are not numbers reutnrs 0's
        :return: Nothing
        """
        print("Testing that change_dates returns 0's for character dates")

        expected = [0, 0, 0]

        date = utils.change_dates(["<very>", "{bad}", "(date) 'given'"])

        for i in range(0, date.__len__()):
            self.assertEqual(round(date[i], 6), round(expected[i], 6))


class PotentialVorticity(unittest.TestCase):
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

        pv_result = utils.potential_vorticity(1, 1)
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
        barotropic_potential_vorticity_vec = np.vectorize(utils.potential_vorticity)
        pv_result = barotropic_potential_vorticity_vec(lat, z_value)
        self.assertTrue(isinstance(pv_result, np.ndarray), "potential vorticity vec is not list")

    def test_bpv_never_returns_0(self):
        """
        Check that potential_vorticity doesn't return 0 if the equation equals
        0. Should return a very small number.
        :return: Nothing
        """
        print("Testing that potential_vorticity returns non-zero...")

        pv_result = utils.potential_vorticity(0, 1)
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
        pv_result = utils.potential_vorticity(lat, z_value)
        self.assertAlmostEqual(pv_result * 1000, -2.452 * 10 ** -8 * 1000)

    def test_bpv_returns_expected_array(self):
        """
        Check that potential_vorticity returns the expected array after calculation
        :return: Nothing
        """
        print("Testing that potential_vorticity returns correct array")

        lat = [-57.996, -56.375, -54.496]
        z_value = [5230, 5223.1, 4447.1]
        barotropic_potential_vorticity_vec = np.vectorize(utils.potential_vorticity)

        pv_result = barotropic_potential_vorticity_vec(lat, z_value)
        pv_expected = [-0.2365 * 10 ** -7, -0.2325 * 10 ** -7, -0.267 * 10 ** -7]

        for i in range(0, lat.__len__()):
            self.assertAlmostEqual(pv_result[i] * 1000, pv_expected[i] * 1000)


class Sorter(unittest.TestCase):
    """
    Sorter test case
    """

    def setUp(self):
        """
        Set up repeated test values
        :return: Nothing
        """
        self.sites = np.array([-1, 0])
        self.msites = np.arange(-1, 1.01, 0.01)

    def test_returns_array(self):
        """
        Check that the sorter returns a numpy array
        :return: Nothing
        """

        print("Testing that sorter returns numpy array")

        sorted_test = utils.sorter(self.sites, self.msites)

        self.assertEqual(type(sorted_test), np.ndarray,
                         "Return type for sorter is incorrect")

    def test_returns_length(self):
        """
        Check that the sorter returns a numpy array of the correct length
        :return: Nothing
        """

        print("Testing that sorter returns numpy array of correct length")

        sorted_test = utils.sorter(self.sites, self.msites)

        self.assertEqual(sorted_test.__len__(), self.msites.__len__(),
                         "Return for sorter is incorrect length")

    def test_returns_ones(self):
        """
        Check that the sorter returns a numpy array of the correct length
        :return: Nothing
        """

        print("Testing that sorter returns numpy array of ones "
              "if boundaries are bad")

        sorted_test = utils.sorter(np.array([-1, -1]), self.msites)

        for i in range(sorted_test.__len__()):
            self.assertEqual(sorted_test[i], 1, "Expected all ones")

    def test_returns_in_bounds(self):
        """
        Check that the sorter returns 0's for in bounds, 1's for out bounds
        :return: Nothing
        """

        print("Testing sorter boundaries")

        sorted_test = utils.sorter(self.sites, self.msites)

        for i in range(math.floor(sorted_test.__len__()/2)):
            self.assertEqual(sorted_test[i], 0, "Expected 0's")

        for i in range(math.floor(sorted_test.__len__()/2), sorted_test.__len__()):
            self.assertEqual(sorted_test[i], 1, "Expected 1's")


class SpatialCorrelation(unittest.TestCase):
    """
    Test cases for "spatial_correlation" function
    """

    def test_spatial_correlation_returns_float(self):
        """
        Check that spatial_correlation returns a float if given a float
        :return: Nothing
        """
        print("Testing that spatial_correlation returns a float")

        sc_result = utils.spatial_correlation(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)

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
        spatial_correlation_vec = np.vectorize(utils.spatial_correlation)
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
        sc_result = utils.spatial_correlation(hist_long, float_long, longitude_large,
                                        hist_lat, float_lat, latitude_large,
                                        hist_dates, float_date, age_large,
                                        phi_large)

        self.assertEqual(round(sc_result, 4), 1.3125)

        # With potential vorticity
        pv_hist = -0.3 * 10 ** -7
        pv_float = -0.35 * 10 ** -7
        sc_pv_result = utils.spatial_correlation(hist_long, float_long, longitude_large,
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
        spatial_correlation_vec = np.vectorize(utils.spatial_correlation)
        sc_expected = [5.158, 13.4935, 14.1804]
        sc_result = spatial_correlation_vec(hist_long, float_long, longitude_large,
                                            hist_lat, float_lat, latitude_large,
                                            hist_dates, float_date, age_large,
                                            phi_large, pv_hist, pv_float)

        for i in range(0, sc_expected.__len__()):
            self.assertEqual(
                round(sc_result[i], 4), sc_expected[i], "spatial_correlation arrays aren't equal"
            )


class WrapLongitude(unittest.TestCase):
    """
    Test cases for wrap_longitude function
    """

    def test_does_not_wrap(self):
        """
        Shouldn't wrap if longitude values are > 40 and < 320
        :return: Nothing
        """
        print("Testing that wrap_longitude does not wrap values unnecessarily")

        long = np.array([50, 100, 150, 68, 300])
        new_long = utils.wrap_longitude(long)

        for i in range(long.__len__()):
            self.assertTrue(long[i] == new_long[i],
                            "Shouldn't wrap values if they are > 40 and < 320")

    def test_wraps_values(self):
        """
        Should wrap values of the longitudes are around 0 or 360
        :return: Nothing
        """
        print("Testing that wrap_longitude wraps values of they are near 0 or 360")

        long = np.array([-10, 20, 30, 350, 340, 300, 359, 100, 35])
        expected_long = np.array([350, 380, 390, 350, 340, 300, 359, 100, 395])
        new_long = utils.wrap_longitude(long)

        for i in range(long.__len__()):
            self.assertTrue(new_long[i] == expected_long[i],
                            "Should wrap values if they are < 40 and > 320")

    def test_length_equal(self):
        """
        This function shouldn't throw away any spatial values. Output length should equal input
        :return: Nothing
        """
        print("Testing that wrap_longitude returns array of equal length to input")

        long = np.array([-100, 200, 35, 352, 34, 397, 359, 100, 35])
        new_long = utils.wrap_longitude(long)

        self.assertTrue(long.__len__() == new_long.__len__(),
                        "input and output should be equal in length")


if __name__ == '__main__':
    unittest.main()
