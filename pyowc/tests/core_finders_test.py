""" Tests for core.finders module functions """
import os
import unittest
import numpy as np
from scipy.io import loadmat
import random

from pyowc import core
from . import TESTS_CONFIG


# pylint: disable=fixme
# pylint: disable=too-many-instance-attributes
class Find10Thetas(unittest.TestCase):
    """
    Test cases for find_10thetas function
    """

    def setUp(self):
        # set up data to run find_10thetas
        float_source = TESTS_CONFIG['TEST_FLOAT_SOURCE']
        mapped_values_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'],
                                               "float_mapped_test", "map_" + float_source + ".mat"])
        source_values_path = os.path.sep.join([TESTS_CONFIG['FLOAT_SOURCE_DIRECTORY'],
                                               float_source + TESTS_CONFIG['FLOAT_SOURCE_POSTFIX']])

        mapped_values = loadmat(mapped_values_path)
        source_values = loadmat(source_values_path)

        self.sal = source_values['SAL']
        self.ptmp = source_values['PTMP']
        self.pres = source_values['PRES']

        self.la_ptmp = mapped_values['la_ptmp']

        # matlab answers (for file 3901960)

        self.tlevels = np.array([1.93961369305644,
                                 1.94046894952240,
                                 1.86741374662372,
                                 1.99059382007223,
                                 1.83113347266278,
                                 2.02694205058300,
                                 1.74856507801603,
                                 1.79263814207055,
                                 2.02025501581947,
                                 1.97075995179013]).reshape(-1, 1)


        self.plevels = np.array([255.099998474121,
                                 605.099998474121,
                                 755.099998474121,
                                 505.099998474121,
                                 805.099998474121,
                                 405.099998474121,
                                 905.099998474121,
                                 855.099998474121,
                                 355.099998474121,
                                 555.099998474121]).reshape(-1, 1)

        self.var_s_thetalevels = np.array([np.nan,
                                           0.0417641466600311,
                                           0.00827684711968173,
                                           0.00190605048358521,
                                           0.00690521338337522,
                                           0.00657688744326181,
                                           0.00598545254829351,
                                           0.00802480111250230,
                                           0.00463859078209821,
                                           0.00686978971670764,
                                           0.00231049559579914,
                                           0.00930018562434973,
                                           0.00791732371803068,
                                           0.00366687591997541,
                                           0.00581845577896215,
                                           0.00639230605471250,
                                           0.00631804077884318,
                                           0.00902606909144420,
                                           0.0137569367353719,
                                           0.0130841767947614,
                                           0.0174757883927057,
                                           0.0133959994423553,
                                           0.0145825527951333,
                                           0.0260367153139685,
                                           0.0286189522655563,
                                           0.0308997771862160,
                                           0.0371860758533865,
                                           0.0383576521771487,
                                           0.0474128181226278,
                                           0.0409629948289543,
                                           0.0423780718895870,
                                           0.0448321203445280,
                                           0.0447966291280737,
                                           0.0465707923888888,
                                           0.0494526560160367,
                                           0.0490136597858479,
                                           0.0603739836990882,
                                           0.0588213891031911,
                                           0.0618893978482367]).reshape(-1, 1)

        self.thetalevels = np.array([np.nan,
                                     1.11807385800664,
                                     1.74921107937534,
                                     1.93961369305644,
                                     1.99500304911626,
                                     2.02025501581947,
                                     2.02694205058300,
                                     2.00848644675449,
                                     1.99059382007223,
                                     1.97075995179013,
                                     1.94046894952240,
                                     1.92705529038417,
                                     1.89375692581917,
                                     1.86741374662372,
                                     1.83113347266278,
                                     1.79263814207055,
                                     1.74856507801603,
                                     1.71174509186099,
                                     1.68506762599664,
                                     1.67653435402351,
                                     1.62733097100766,
                                     1.57739642670448,
                                     1.52565233494691,
                                     1.47091411129001,
                                     1.41236755919116,
                                     1.36299326635668,
                                     1.30830169888757,
                                     1.25487243189936,
                                     1.20028176877869,
                                     1.14931046649699,
                                     1.09558302804420,
                                     1.04691797009209,
                                     1.00240760777759,
                                     0.954500074838010,
                                     0.910028348222131,
                                     0.866025315061063,
                                     0.823372829184810,
                                     0.781478879297776,
                                     0.741563524245540]).reshape(-1, 1)

    def test_theta_levels_shape(self):
        """
        Check that we get 10 levels
        :return: Nothing
        """
        print("Testing that we only get 10 theta levels")
        t_levels, p_levels, index, var_sal_theta, theta_levels = core.finders.find_10thetas(self.sal,
                                                                                            self.ptmp,
                                                                                            self.pres,
                                                                                            self.la_ptmp,
                                                                                            [], [], [], [], 0.5)

        self.assertEqual(t_levels.shape, self.tlevels.shape,
                         "Got incorrect number of theta levels")
        self.assertEqual(p_levels.shape, self.plevels.shape,
                         "Got incorrect number of pressure levels")
        self.assertEqual(index.shape, (10, self.sal.shape[1]),
                         "have incorrect number of indices")
        self.assertEqual(var_sal_theta.shape, self.var_s_thetalevels.shape,
                         "Got incorrect number of theta levels")
        self.assertEqual(theta_levels.shape, self.thetalevels.shape,
                         "Got incorrect number of theta levels")

    def test_theta_levels_values(self):
        """
        Check that we get 10 levels
        :return: Nothing
        """
        print("Testing that we only get 10 theta levels")
        t_levels, p_levels, index, var_sal_theta, theta_levels = core.finders.find_10thetas(self.sal,
                                                                                            self.ptmp,
                                                                                            self.pres,
                                                                                            self.la_ptmp,
                                                                                            [], [], [], [], 0.5)

        for i in range(t_levels.__len__()):
            self.assertAlmostEqual(t_levels[i, 0], self.tlevels[i, 0], 10,
                                   "Got incorrect theta level")
            self.assertAlmostEqual(p_levels[i, 0], self.plevels[i, 0], 10,
                                   "Got incorrect pressure level")

        for i in range(self.var_s_thetalevels.shape[0]):
            for j in range(self.var_s_thetalevels.shape[1]):
                if not np.isnan(var_sal_theta[i, j]):
                    self.assertAlmostEqual(var_sal_theta[i, j], self.var_s_thetalevels[i, j], 10,
                                           "salinity variance is incorrect")

        for i in range(self.thetalevels.shape[0]):
            for j in range(self.thetalevels.shape[1]):
                if not np.isnan(theta_levels[i, j]):
                    self.assertAlmostEqual(theta_levels[i, j], self.thetalevels[i, j], 10,
                                           "salinity variance is incorrect")

        self.assertGreater(index.__len__(), 0, "Should have some indices")

    def test_theta_levels_mid(self):
        """
                Check that we get 10 levels
                :return: Nothing
                """
        print("Testing that we can include the middle band of data")
        t_levels, p_levels, index, var_sal_theta, theta_levels = core.finders.find_10thetas(
            self.sal,
            self.ptmp,
            self.pres,
            self.la_ptmp,
            [10], [-10],
            [1750], [500],
            0.5)

        for i in t_levels:
            self.assertTrue(-10 < i < 10, "No theta levels should be in this range")

        for i in p_levels:
            self.assertTrue(500 < i < 1750, "No pressure levels should be in this range")

        self.assertEqual(index.__len__(), 10, "should get 10 levels for bounded data")

        self.assertGreater(var_sal_theta.__len__(), 0,
                           "should get 10 levels for bounded data")
        self.assertGreater(theta_levels.__len__(), 0,
                           "should get 10 levels for bounded data")

    def test_theta_levels_no_mid(self):
        """
        Check that we get 10 levels
        :return: Nothing
        """
        print("Testing that we can exclude the middle band of data")
        t_levels, p_levels, index, var_sal_theta, theta_levels = core.finders.find_10thetas(
            self.sal,
            self.ptmp,
            self.pres,
            self.la_ptmp,
            [-0.1], [0.1],
            [1000], [1100],
            0.5)

        self.assertEqual(index.__len__(), 10, "should get 10 levels for bounded data")

        for i in t_levels:
            self.assertFalse(-0.1 < i < 0.1, "No theta levels should be in this range")

        for i in p_levels:
            self.assertFalse(1000 < i < 1100, "No pressure levels should be in this range")

        # TODO (15/04/21) these test are already called within test_theta_levels_shape maybe these should be remove
        self.assertEqual(var_sal_theta.shape, self.var_s_thetalevels.shape,
                         "Got incorrect number of theta levels")
        self.assertEqual(theta_levels.shape, self.thetalevels.shape,
                         "Got incorrect number of theta levels")

    def test_theta_levels_bounds(self):
        """
                Check that we get 10 levels
                :return: Nothing
                """
        print("Testing that we can exclude data above/below points")
        t_levels, p_levels, index, var_sal_theta, theta_levels = core.finders.find_10thetas(
            self.sal,
            self.ptmp,
            self.pres,
            self.la_ptmp,
            [], [0.1],
            [], [500],
            0.5)

        for i in p_levels:
            self.assertGreater(i, 500, "pressure levels should exceed lower bound")
        for i in t_levels:
            self.assertGreater(i, 0.1, "pressure levels should not exceed upper bound")

        self.assertEqual(index.__len__(), 10, "should get 10 levels")

        self.assertTrue(var_sal_theta.__len__() > 0,
                        "should have variance of salinity on each level")
        self.assertTrue(theta_levels.__len__() > 0,
                        "should have theta levels")

        t_levels, p_levels, index, var_sal_theta, theta_levels = core.finders.find_10thetas(
            self.sal,
            self.ptmp,
            self.pres,
            self.la_ptmp,
            [2], [],
            [1500], [],
            0.5)

        for i in p_levels:
            self.assertLess(i, 1500, "pressure levels should not exceed upper bound")

        for i in t_levels:
            self.assertLess(i, 2, "pressure levels should not exceed upper bound")

        self.assertEqual(index.__len__(), 10, "should get 10 levels")

        self.assertTrue(var_sal_theta.__len__() > 0,
                        "should have variance of salinity on each level")
        self.assertTrue(theta_levels.__len__() > 0,
                        "should have theta levels")

    def test_less_50_bar(self):
        """
        Check that we still get 10 thetas if intervals smaller than 50 bar
        :return: Nothing
        """
        print("Testing we still get theta levels at <50 bar intervals")
        t_levels, p_levels, index, var_sal_theta, theta_levels = core.finders.find_10thetas(
            self.sal,
            self.ptmp,
            self.pres,
            self.la_ptmp,
            [], [],
            [700], [500],
            0.5)

        self.assertTrue(theta_levels.__len__() > 0, "should get levels")


class Find25Boxes(unittest.TestCase):
    """
    Test cases for find_25boxes function"
    """

    def setUp(self):
        """
        Setting up variables to be used in test cases
        Is run before each test case
        :return: Nothing
        """
        self.pa_wmo_boxes = loadmat(
            os.path.sep.join([TESTS_CONFIG['CONFIG_DIRECTORY'], TESTS_CONFIG['CONFIG_WMO_BOXES']]))
        self.pn_float_long = 57.1794
        self.pn_float_lat = -59.1868

    def test_can_load_matrix(self):
        """
        Check that we can open a matlab matrix file (.mat)
        :return: Nothing
        """
        print("Testing that a matlab matrix can be loaded in...")
        self.assertNotEqual(self.pa_wmo_boxes, '', "Failed to load matlab matrix")

    def test_nearest_neighbour_nan(self):
        """
        Check that the nearest neighbour algorithm returns grid point
        :return: Nothing
        """
        print("Testing nearest_neighbour returns correct grid point")
        x_axis = np.arange(0, 3, 1, int)
        y_axis = np.arange(0, 3, 1, int)
        data = np.full((1, 9), np.arange(3, 12)).reshape(3, 3)
        self.assertEqual(
            core.finders.nearest_neighbour(
                x_axis, y_axis, data, 0, 0
            ), 3, "Nearest neighbour interpolation is wrong"
        )

    def test_nearest_neighbour_returns_point_in_grid(self):
        """
        Check that, even if the nearest neighbour function is given a number outside
        the grid, it still returns a value that exists within the grid
        :return: Nothing
        """
        print("Testing nearest_neighbour returns data from the grid")
        x_axis = np.arange(0, 10000, 5, int)
        y_axis = np.arange(0, 10000, 5, int)
        data = np.full((1, 10000 * 10000), np.arange(1, 10000 * 10000 + 1)).reshape(10000, 10000)
        self.assertIn(
            core.finders.nearest_neighbour(
                x_axis, y_axis, data, 93820495.3958, -9283499184.439
            ), data, "Answer is not in data grid"
        )

    def test_returned_array_is_correct_size(self):
        """
        Check that we get back 25 boxes, each with 4 values
        :return: Nothing
        """
        print("Testing that the returned value is 25x4 matrix")
        pa_wmo_numbers = core.finders.find_25boxes(self.pn_float_long, self.pn_float_lat, self.pa_wmo_boxes)
        self.assertEqual(pa_wmo_numbers.shape, (25, 4), 'returned matrix shape is incorrect')

    def test_returns_array_of_nan_for_nan_input(self):
        """
        Check that even if we get empty inputs, we still receive boxes
        :return: Nothing
        """
        print("Testing nan input")
        pa_wmo_numbers = core.finders.find_25boxes(np.nan, np.nan, self.pa_wmo_boxes)
        self.assertEqual(
            np.allclose(
                pa_wmo_numbers, np.full((25, 4), np.nan), equal_nan=True
            ), True, "nan input doesn't lead to nan output"
        )

    def test_returned_value_exists_in_wmo_boxes(self):
        """
        Check that each of the 25 boxes we receive definitely exists in the original
        collection of wmo boxes
        :return:
        """
        print("Testing returned value exists in wmo_boxes.mat")
        pa_wmo_numbers = core.finders.find_25boxes(self.pn_float_long, self.pn_float_lat, self.pa_wmo_boxes)
        result = np.full((25, 1), False)

        count = 0
        for i in np.array(pa_wmo_numbers):
            for j in np.array(self.pa_wmo_boxes.get('la_wmo_boxes')):
                if i[0] == j[0] and i[1] == j[1] and i[2] == j[2] and i[3] == j[3]:
                    result[count] = True
            count += 1

        self.assertEqual(
            np.all([result, np.full((25, 1), True)]), True, "Some values don't match wmo_boxes.mat"
        )


class FindBestHist(unittest.TestCase):
    """
    Test cases for find_besthist function
    """

    # pylint: disable=too-many-instance-attributes
    def setUp(self):
        """
        Set up for the tests. Creates 1000 historical random data points that
        are the float data points +/- 5. This could fail randomly if, by pure
        chance, none of the generated data is in the ellipse, so use pseudo random
        numbers just in case
        :return: nothing
        """
        random.seed(1)
        self.grid_lat = np.random.rand(1000) * 5 * random.choice([-1, 1]) + -59.1868
        self.grid_long = np.random.rand(1000) * 5 * random.choice([-1, 1]) + 57.1794
        self.grid_dates = np.random.rand(1000) * 5 * random.choice([-1, 1]) + 5.1083 * 10 ** 3
        self.grid_z_values = np.random.rand(1000) * 5 * random.choice([-1, 1]) + 2.018 * 10 ** 3
        self.lat = -59.1868
        self.long = 57.1794
        self.z_value = 2.018 * 10 ** 3
        self.date = 5.1083 * 10 ** 3
        self.lat_large = 4
        self.lat_small = 2
        self.long_large = 8
        self.long_small = 4
        self.phi_large = 0.5
        self.phi_small = 0.1
        self.age_large = 20
        self.age_small = 10
        self.map_pv_usage = 1
        self.max_casts = 300

    def test_returns_array(self):
        """
        Check that find_besthist gives back an array
        :return: Nothing
        """
        print("Testing that find_besthist gives back a numpy array")

        index = core.finders.find_besthist(self.grid_lat, self.grid_long, self.grid_dates, self.grid_z_values,
                                           self.lat, self.long, self.date, self.z_value,
                                           self.lat_large, self.lat_small, self.long_large, self.long_small,
                                           self.phi_large, self.phi_small, self.age_large, self.age_small,
                                           self.map_pv_usage, self.max_casts)

        self.assertTrue(isinstance(index, np.ndarray), "find_besthist not returning array")

    def test_returns_empty(self):
        """
        Check that find_besthist gives back an empty array if no data fits inside
        the ellipse
        :return: Nothing
        """
        print("Testing that find_besthist gives back an empty array if no data inside ellipse")

        grid_lat = np.array([1, 1])
        grid_long = np.array([10, 10])
        grid_z_values = np.array([10, 10])
        grid_dates = np.array([0, 0])

        index = core.finders.find_besthist(grid_lat, grid_long, grid_dates, grid_z_values,
                                           self.lat, self.long, self.date, self.z_value,
                                           self.lat_large, self.lat_small, self.long_large, self.long_small,
                                           self.phi_large, self.phi_small, self.age_large, self.age_small,
                                           self.map_pv_usage, self.max_casts)

        self.assertTrue(index.__len__() == 0, "No data should have been selected")

    def test_returns_array_of_right_size(self):
        """
        Check that find_besthist gives back an array that is the size we asked
        :return: Nothing
        """
        print("Testing that find_besthist gives back the right sized array")

        index = core.finders.find_besthist(self.grid_lat, self.grid_long, self.grid_dates, self.grid_z_values,
                                           self.lat, self.long, self.date, self.z_value,
                                           self.lat_large, self.lat_small, self.long_large, self.long_small,
                                           self.phi_large, self.phi_small, self.age_large, self.age_small,
                                           self.map_pv_usage, self.max_casts)

        self.assertTrue(index.__len__() == self.max_casts, "index is incorrect size")

    def test_returns_expected_values(self):
        """
        Check that find_besthist gives us the values we expect
        Since 1/3 of the data could be selected randomly, we will just check
        that it removes data that isn't inside the ellipse
        :return: Nothing
        """
        print("Testing that find_besthist gives back an array containing expected values")
        grid_lat = np.array([1, 1, -60, -58])
        grid_long = np.array([10, 10, 59, 57])
        grid_dates = np.array([1, 1, 5.108 * 10 ** 3, 5.109 * 10 ** 3])
        grid_z_values = np.array([1, 1, 2.02 * 10 ** 3, 2.015 * 10 ** 3])
        expected = np.array([2, 3])

        index = core.finders.find_besthist(grid_lat, grid_long, grid_dates, grid_z_values,
                                           self.lat, self.long, self.date, self.z_value,
                                           self.lat_large, self.lat_small, self.long_large, self.long_small,
                                           self.phi_large, self.phi_small, self.age_large, self.age_small,
                                           self.map_pv_usage, self.max_casts)

        for i in range(0, index.__len__()):
            self.assertTrue(index[i] == expected[i], "output is incorrect (wrong index selected)")


class FindEllipse(unittest.TestCase):
    """
    Test cases for "find ellipse" function
    """

    def test_find_ellipse_returns_float(self):
        """
        Check that find_ellipse returns a float if only given floats
        :return: Nothing
        """
        print("Testing that find_ellipse will return a float")

        ellipse = core.finders.find_ellipse(1, 1, 1, 1, 1, 1, 1, 1)

        self.assertTrue(isinstance(ellipse, float), "find ellipse did not return float")

    def test_find_ellipse_returns_array(self):
        """
        Check that find_ellipse returns an array if vectorised
        :return: Nothing
        """
        print("Testing that find_ellipse vectorised returns an array")

        find_ellipse_vec = np.vectorize(core.finders.find_ellipse)
        ellipse = find_ellipse_vec([1, 1], 1, 1, [1, 1], 1, 1, 1)

        self.assertTrue(isinstance(ellipse, np.ndarray), "find_ellipse did not return array")

    def test_find_ellipse_returns_expected_float(self):
        """
        Check that find_ellipse returns the correct float
        :return: Nothing
        """
        print("Testing that find_ellipse returns the right answer")

        ellipse = core.finders.find_ellipse(53.195, 57.1794, 8,
                                            -57.996, -59.1868, 4,
                                            0.5, -2.3547 * 10 ** -8, -2.452 * 10 ** -8)

        self.assertEqual(round(ellipse, 4), 0.2017, "find_ellipse returned incorrect float")

    def test_find_ellipse_returns_expected_float_without_pv(self):
        """
        Check that find_ellipse returns the correct float without potential vorticity
        :return: Nothing
        """
        print("Testing that find_ellipse returns the right answer without potential vorticity")

        ellipse = core.finders.find_ellipse(53.195, 57.1794, 8,
                                            -57.996, -59.1868, 4,
                                            0.5)

        self.assertEqual(round(ellipse, 4), 0.1934, "find_ellipse returned incorrect float "
                                                    "without potential vortcity")

    def test_find_ellipse_returns_expected_array(self):
        """
        Check that find_ellipse returns the correct array
        :return: Nothing
        """
        print("Testing that find_ellipse returns the right array answer")

        ellipse_expected = [0.2018, 0.3286, 0.4428]
        find_ellipse_vec = np.vectorize(core.finders.find_ellipse)
        hist_long = [53.195, 51.954, 53.107]
        hist_lat = [-57.996, -56.375, -54.496]
        hist_pv = [-0.2354 * 10 ** -7,
                   -0.2325 * 10 ** -7,
                   -0.267 * 10 ** -7]
        ellipse = find_ellipse_vec(hist_long, 57.1794, 8,
                                   hist_lat, -59.1868, 4,
                                   0.5, hist_pv, -2.452 * 10 ** -8)

        for i in range(0, hist_long.__len__()):
            self.assertEqual(round(ellipse[i], 4), ellipse_expected[i],
                             "incorrect result in array")

    def test_find_ellipse_returns_expected_array_without_pv(self):
        """
        Check that find_ellipse returns the correct array
        :return: Nothing
        """
        print("Testing that find_ellipse returns the right array answer")

        ellipse_expected = [0.1934, 0.3199, 0.4261]
        find_ellipse_vec = np.vectorize(core.finders.find_ellipse)
        hist_long = [53.195, 51.954, 53.107]
        hist_lat = [-57.996, -56.375, -54.496]
        ellipse = find_ellipse_vec(hist_long, 57.1794, 8,
                                   hist_lat, -59.1868, 4,
                                   0.5)

        for i in range(0, hist_long.__len__()):
            self.assertEqual(round(ellipse[i], 4), ellipse_expected[i],
                             "incorrect result in array")
