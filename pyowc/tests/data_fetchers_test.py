""" Tests for data.fetchers module functions """

import os
import unittest
import numpy as np
import warnings

# from pyowc import data
from pyowc.data.fetchers import get_data, get_region_data, get_region_hist_locations
from . import TESTS_CONFIG


class GetData(unittest.TestCase):
    """
    Test cases for get_data function
    """

    def setUp(self):
        """
        Sets up some constant variables for testing
        :return: Nothing
        """
        self.float_name = TESTS_CONFIG['TEST_FLOAT_SOURCE']
        self.config = TESTS_CONFIG

    def test_returns_no_data(self):
        """
        If given some data to find, it returns the correct data
        :return: Nothing
        """
        print("Testing that get_data will return no data if box is all 0's")

        wmo_boxes = np.array([self.config['TEST_FLOAT_WMO_BOXES'][0], 0, 0, 0])

        this_data = get_data(wmo_boxes, 1, self.config, self.float_name)
        self.assertTrue(this_data.__len__() == 0, "should return no data")

    def test_returns_correct_data(self):
        """
        See if it gets the right type of data each time
        :return: Nothing
        """
        print("Testing that get_data returns the expected amount of data")

        wmo_boxes = np.array([self.config['TEST_FLOAT_WMO_BOXES'][0], 1, 1, 1])

        data_ctd = get_data(wmo_boxes, 1, self.config, self.float_name)
        self.assertTrue(data_ctd['long'].shape[0] == self.config['TEST_FLOAT_WMO_BOXES_Nhist'][0],
                        "should return some data when fetching ctd")

        data_bot = get_data(wmo_boxes, 2, self.config, self.float_name)
        self.assertTrue(data_bot['long'].shape[0] == self.config['TEST_FLOAT_WMO_BOXES_Nhist'][1],
                        "should return some data when fetching argo")

        data_argo = get_data(wmo_boxes, 3, self.config, self.float_name)
        self.assertTrue(data_argo['long'].__len__() == self.config['TEST_FLOAT_WMO_BOXES_Nhist'][2],
                        "should return some data when fetching argo")
    def test_removes_argo_float(self):
        """
        See if it removes the argo float currently being analysed
        :return: Nothing
        """
        print("Testing that get_data will remove the argo being analysed")

        wmo_boxes = np.array([self.config['TEST_FLOAT_WMO_BOXES'][0], 0, 0, 1])
        float_removed = self.config['TEST_FLOAT_IN_HIST']

        data_normal = get_data(wmo_boxes, 3, self.config, self.float_name)
        data_removed = get_data(wmo_boxes, 3, self.config, float_removed)

        self.assertTrue(data_normal['long'].__len__() - data_removed['long'].__len__() ==
                        self.config['TEST_FLOAT_N_TO_REMOVE'],
                        "Should have removed data associated with the float being processed")

    def test_returns_correct_shape(self):
        """
        See if the returned data has the expected shape
        :return: Nothing
        """
        print("Testing that get_data returns data with the expected shape")

        wmo_boxes = np.array([self.config['TEST_FLOAT_WMO_BOXES'][0], 1, 0, 0])
        this_data = get_data(wmo_boxes, 1, self.config, self.float_name)

        self.assertTrue(this_data['sal'].shape == this_data['ptmp'].shape == this_data['temp'].shape,
                        "Ocean characteristic data should be the same shape")
        self.assertTrue(this_data['long'].shape == this_data['lat'].shape == this_data['dates'].shape,
                        "Spatial/temporal data should be the same shape")
        self.assertTrue(this_data['sal'].shape[1] == this_data['long'].shape[0],
                        "Should be a profile for every location")


class GetRegionData(unittest.TestCase):
    """
    Test cases for get_region_data function
    """

    def setUp(self):
        self.float_name = TESTS_CONFIG['TEST_FLOAT_SOURCE']
        self.config = TESTS_CONFIG
        self.wmo_boxes = np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 1, 1, 1],
                                   [self.config['TEST_FLOAT_WMO_BOXES'][1], 1, 1, 1]])
        self.index = np.array([0, 4, 5, 11, 13, 15, 20, 21, 42])
        self.pres = np.array([3, 5, 15.1, 25.1, 36, 40, 45, 46, 500000])

    def test_returns_6(self):
        """
        check that the function returns 6 arrays
        :return: Nothing
        """
        print("Testing that get_region_data returns 6 arrays")

        test = get_region_data(self.wmo_boxes, self.float_name, self.config,
                               self.index, self.pres)

        self.assertTrue(test.__len__() == 6, "Should return 6 arrays")

    def test_return_shape(self):
        """
        check that the 6 arrays are the expected shape
        :return: Nothing
        """
        print("Testing that get_region_data return values are the correct shape")

        test = get_region_data(self.wmo_boxes, self.float_name, self.config,
                               self.index, self.pres)

        self.assertTrue(test['grid_sal'].shape == test['grid_ptmp'].shape == test['grid_pres'].shape,
                        "salinity, pressure, and potential temperature "
                        "arrays should be same shape")
        self.assertTrue(test['grid_lat'].shape == test['grid_long'].shape == test['grid_dates'].shape,
                        "longitude, latitude, and date arrays should be the same shape")
        self.assertTrue(test['grid_sal'].shape[1] == test['grid_ptmp'].shape[1] == test['grid_pres'].shape[1] ==
                        test['grid_lat'].shape[0] == test['grid_long'].shape[0] == test['grid_dates'].shape[0] ==
                        self.index.__len__(),
                        "Should get the same number of casts as we have asked for")

    def test_gets_different_data(self):
        """
        Check that we can fetch different combinations of data
        :return: Nothing
        """
        print("Testing that get_region_data can return different data types")

        test_ctd = get_region_data(np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 1, 0, 0]]),
                                                 self.float_name, self.config, self.index, self.pres)
        test_bot = get_region_data(np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 0, 1, 0]]),
                                                 self.float_name, self.config, self.index, self.pres)
        test_argo = get_region_data(np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 0, 0, 1]]),
                                                  self.float_name, self.config, self.index, self.pres)

        self.assertTrue(test_ctd['grid_sal'].shape[1] != test_argo['grid_sal'].shape[1],
                        "Should get a different data set, if we have specified it")
        self.assertTrue(test_bot['grid_sal'].shape[1] != test_argo['grid_sal'].shape[1],
                        "Should get a different data set, if we have specified it")

    def test_should_only_get_specific_indices(self):
        """
        If we only give n indices, we should get n data points back
        :return: Nothing
        """
        print("Testing that get_region returns correct amount of data")

        test_many = get_region_data(self.wmo_boxes, self.float_name, self.config, self.index, self.pres)

        self.assertTrue(test_many['grid_sal'].shape[1] == self.index.__len__())

        test_one = get_region_data(self.wmo_boxes, self.float_name, self.config,
                                   [50], self.pres)

        self.assertTrue(test_one['grid_sal'].shape[1] == 1)

    def test_raise_exception_bad_indices(self):
        """
        Check that, if only bad indices are given, an exception is raised
        :return: Nothing
        """
        print("Testing exception is raised if indices are bad")

        with self.assertRaises(Exception) as no_index:
            get_region_data(self.wmo_boxes, self.float_name, self.config,
                            [], self.pres)

        self.assertTrue('NO DATA FOUND' in str(no_index.exception))

        with self.assertRaises(Exception) as big_index:
            get_region_data(self.wmo_boxes, self.float_name, self.config,
                            [99999999999999999], self.pres)

        self.assertTrue('NO DATA FOUND' in str(big_index.exception))


class GetRegionHistLoc(unittest.TestCase):
    """
    Test cases for get_hist_region_locations function
    """

    def setUp(self):
        """
        Set up values that will be used for testing
        :return: Nothing
        """
        self.config = TESTS_CONFIG
        self.float_name = TESTS_CONFIG['TEST_FLOAT_IN_HIST']
        self.wmo_box = np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 1, 0, 0]])
        self.wmo_boxes = np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 1, 1, 1],
                                   [self.config['TEST_FLOAT_WMO_BOXES'][1], 1, 1, 1]])

    def test_returns_three(self):
        """
        Check that the function returns 3 unique values
        :return: Nothing
        """
        print("Testing that get_region_his_locations returns 3 unique values of equal length")

        lat, long, age = get_region_hist_locations(self.wmo_box, self.float_name, self.config)
        self.assertEqual(lat.__len__(), long.__len__(),
                         "Should have equal numbers of latitude and longitude")
        self.assertEqual(lat.__len__(), age.__len__(),
                         "should have equal number of spatial and temporal data")

        for i in range(0, lat.__len__()):
            self.assertNotEqual(lat[i], long[i], "latitude and longitude should be unique")
            self.assertNotEqual(lat[i], age[i], "spatial and temporal data should be unique")

    def test_current_float_removed(self):
        """
        Check that if the argo float currently being processed appears in the historical data that
        it is removed
        :return: Nothing
        """
        print("Testing that get_region_his_locations removes the current float")

        wmo_box_argo = np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 0, 0, 1]])

        # Data not including the float:
        lat_rem, long_rem, age_rem = get_region_hist_locations(wmo_box_argo,
                                                               self.float_name,
                                                               self.config)
        # All data:
        lat_no_rem, long_no_rem, age_no_rem = get_region_hist_locations(wmo_box_argo,
                                                                        'DUMMY_WMO',
                                                                        self.config)

        self.assertNotEqual(lat_rem.__len__(), lat_no_rem.__len__(),
                            "Did not remove current argo float")

        self.assertGreater(long_no_rem.__len__(), long_rem.__len__(),
                           "Should have more values for the data with no argo removal")

        self.assertEqual(age_no_rem.__len__() - age_rem.__len__(), self.config['TEST_FLOAT_N_TO_REMOVE'],
                         "should have removed %i values" % self.config['TEST_FLOAT_N_TO_REMOVE'])

    def test_no_data(self):
        """
        Check that if we receive no data then we return the expected values
        :return: Nothing
        """
        print("Testing that get_region_hist_locations returns expected values for no data")

        wmo_boxes_no_data = np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 0, 0, 0]])

        with self.assertRaises(ValueError) as no_data:
            get_region_hist_locations(wmo_boxes_no_data, 'none', self.config)

        self.assertTrue('get_region_hist_locations found no data for your specification. '
                        'Are your wmo_boxes files set up correctly?' in str(no_data.exception))

    def test_can_choose_data(self):
        """
        Check that, by applying different flags, we can fetch differing data types
        :return: Nothing
        """
        print("Testing that get_region_his_locations can fetch different data types")

        wmo_box_argo = np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 0, 0, 1]])
        wmo_box_bot = np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 0, 1, 0]])
        wmo_box_ctd = np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 1, 0, 0]])

        data_argo = get_region_hist_locations(wmo_box_argo, 'none', self.config)
        data_bot = get_region_hist_locations(wmo_box_bot, 'none', self.config)
        data_ctd = get_region_hist_locations(wmo_box_ctd, 'none', self.config)

        for i in range(0, data_argo.__len__()):
            self.assertNotEqual(data_argo[i].__len__(), data_bot[i].__len__(),
                                "Should have different data for argo and bottle")
            self.assertNotEqual(data_argo[i].__len__(), data_ctd[i].__len__(),
                                "Should have different data for argo and ctd")
            # no bottle database reference,
            self.assertNotEqual(data_bot[i].__len__(), data_ctd[i].__len__(),
                                "Should have different data for bottle and ctd")

    def test_can_combine_data(self):
        """
        Test that flags can be set to retrieve all the data
        :return: Nothing
        """
        print("Testing that get_region_hist_locations can fetch all the data")

        wmo_box_all = np.array([[self.config['TEST_FLOAT_WMO_BOXES'][0], 1, 1, 1]])
        all_data = get_region_hist_locations(wmo_box_all, 'none', self.config)

        for data in all_data:
            self.assertEqual(data.__len__(), self.config['TEST_FLOAT_N_DATA'],
                             "Should get the same amount of data as matlab")

    def test_can_combine_boxes(self):
        """
        Test that it works with multiple boxes
        :return: Nothing
        """
        print("Testing that get_region_hist_loc will fetch data from multiple WMO boxes")

        all_data = get_region_hist_locations(self.wmo_boxes, 'none', self.config)
        some_data = get_region_hist_locations([self.wmo_boxes[0]], 'none', self.config)
        half_data = get_region_hist_locations([self.wmo_boxes[0], [self.config['TEST_FLOAT_WMO_BOXES'][1], 0, 0, 0]],
                                              'none', self.config)

        for i in range(0, all_data.__len__()):
            self.assertGreater(all_data[i].__len__(), some_data[i].__len__(),
                               "Should have more data when using multiple boxes")
            self.assertEqual(some_data[i].__len__(), half_data[i].__len__(),
                             "Should get equal data with using one ox and using "
                             "two boxes, but getting no data from second box")


if __name__ == '__main__':
    unittest.main()
