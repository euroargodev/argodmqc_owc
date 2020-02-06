"""
-----Get Historical Region Location Data Test File-----

Written by: Edward Small
When: 05/12/2019

Contains unit tests to check the functionality of the `get_hist_region_locations` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from ow_calibration.get_region.get_region_hist_locations import get_region_hist_locations


class MyTestCase(unittest.TestCase):
    """
    Test cases for get_hist_region_locations function
    """

    def setUp(self):
        """
        Set up values that will be used for testing
        :return: Nothing
        """
        self.config = {'HISTORICAL_DIRECTORY': "data/climatology",
                       'HISTORICAL_CTD_PREFIX': "/historical_ctd/ctd_",
                       'HISTORICAL_BOTTLE_PREFIX': "/historical_bot/bot_",
                       'HISTORICAL_ARGO_PREFIX': "/historical_argo/argo_"}
        self.float_name = '1900193'
        self.wmo_box = np.array([[3505, 1, 0, 0]])
        self.wmo_boxes = np.array([[3505, 1, 1, 1], [3506, 1, 1, 1]])

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

        wmo_box_argo = np.array([[3505, 0, 0, 1]])

        lat_rem, long_rem, age_rem = get_region_hist_locations(wmo_box_argo,
                                                               self.float_name,
                                                               self.config)
        lat_no_rem, long_no_rem, age_no_rem = get_region_hist_locations(wmo_box_argo,
                                                                        'DOES NOT EXIST',
                                                                        self.config)

        self.assertNotEqual(lat_rem.__len__(), lat_no_rem.__len__(),
                            "Did not remove current argo float")

        self.assertGreater(long_no_rem.__len__(), long_rem.__len__(),
                           "Should have more values for the data with no argo removal")

        self.assertEqual(age_no_rem.__len__() - age_rem.__len__(), 30,
                         "should have removed 30 values")

    def test_no_data(self):
        """
        Check that if we receive no data then we return the expected values
        :return: Nothing
        """
        print("Testing that get_region_hist_locations returns expected values for no data")

        wmo_boxes_no_data = np.array([[3505, 0, 0, 0]])
        lat_no_data, long_no_data, age_no_data = get_region_hist_locations(wmo_boxes_no_data,
                                                                           'none',
                                                                           self.config)

        self.assertEqual(lat_no_data, long_no_data, "latitude and longitude should be equal (999)")
        self.assertEqual(age_no_data, 'NaN', "age should be NaN")

    def test_can_choose_data(self):
        """
        Check that, by applying different flags, we can fetch differing data types
        :return: Nothing
        """
        print("Testing that get_region_his_locations can fetch different data types")

        wmo_box_argo = np.array([[3505, 0, 0, 1]])
        wmo_box_bot = np.array([[3505, 0, 1, 0]])
        wmo_box_ctd = np.array([[3505, 1, 0, 0]])

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

        wmo_box_all = np.array([[3505, 1, 1, 1]])
        all_data = get_region_hist_locations(wmo_box_all, 'none', self.config)

        for data in all_data:
            self.assertEqual(data.__len__(), 830, "Should get the same amount of data as matlab")

    def test_can_combine_boxes(self):
        """
        Test that it works with multiple boxes
        :return: Nothing
        """
        print("Testing that get_region_hist_loc will fetch data from multiple WMO boxes")

        all_data = get_region_hist_locations(self.wmo_boxes, 'none', self.config)
        some_data = get_region_hist_locations([self.wmo_boxes[0]], 'none', self.config)
        half_data = get_region_hist_locations([self.wmo_boxes[0], [3506, 0, 0, 0]],
                                              'none', self.config)

        for i in range(0, all_data.__len__()):
            self.assertGreater(all_data[i].__len__(), some_data[i].__len__(),
                               "Should have more data when using multiple boxes")
            self.assertEqual(some_data[i].__len__(), half_data[i].__len__(),
                             "Should get equal data with using one ox and using "
                             "two boxes, but getting no data from second box")


if __name__ == '__main__':
    unittest.main()
