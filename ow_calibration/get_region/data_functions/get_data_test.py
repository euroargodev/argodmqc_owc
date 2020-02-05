"""
-----Get Data Test File-----

Written by: Edward Small
When: 05/12/2019

Contains unit tests to check the functionality of the `get_data` function.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from ow_calibration.get_region.data_functions.get_data import get_data


class GetDataTestCase(unittest.TestCase):
    """
    Test cases for get_data function
    """

    def setUp(self):
        """
        Sets up some constant variables for testing
        :return: Nothing
        """

        self.float_name = "3901960"
        # self.float_name_removed = "1900479"
        self.config = {'HISTORICAL_DIRECTORY': "data/climatology",
                       'HISTORICAL_CTD_PREFIX': "/historical_ctd/ctd_",
                       'HISTORICAL_BOTTLE_PREFIX': "/historical_bot/bot_",
                       'HISTORICAL_ARGO_PREFIX': "/historical_argo/argo_",
                       'MAP_P_DELTA': 250}

    def test_returns_no_data(self):
        """
        If given some data to find, it returns the correct data
        :return: Nothing
        """
        print("Testing that get_data will return no data if box is all 0's")

        wmo_boxes = np.array([3505, 0, 0, 0])

        data = get_data(wmo_boxes, 1, self.config, self.float_name)
        self.assertTrue(data.__len__() == 0, "should return no data")

    def test_returns_correct_data(self):
        """
        See if it gets the right type of data each time
        :return: Nothing
        """
        print("Testing that get_data returns the expected amount of data")

        wmo_boxes = np.array([3505, 1, 1, 1])

        data_ctd = get_data(wmo_boxes, 1, self.config, self.float_name)
        self.assertTrue(data_ctd['long'].shape[1] == 10, "should return some data when fetching ctd")

        data_bot = get_data(wmo_boxes, 2, self.config, self.float_name)
        self.assertTrue(data_bot['long'].shape[1] == 33, "should return some data when fetching argo")

        data_argo = get_data(wmo_boxes, 3, self.config, self.float_name)
        self.assertTrue(data_argo['long'][0].__len__() == 787, "should return some data when fetching argo")

    def test_removes_argo_float(self):
        """
        See if it removes the argo float currently being analysed
        :return: Nothing
        """
        print("Testing that get_data will remove the argo being analysed")

        wmo_boxes = np.array([3505, 0, 0, 1])
        float_removed = "1900479"

        data_normal = get_data(wmo_boxes, 3, self.config, self.float_name)
        data_removed = get_data(wmo_boxes, 3, self.config, float_removed)

        # 7 pieces of historical data are from float 1900479 in the selected box, whereas
        # 0 pieces of historical data are from float 3901960 in the selected box, so there
        # should be 7 less data in data_removed

        self.assertTrue(data_normal['long'][0].__len__() - data_removed['long'][0].__len__() == 7,
                        "Should have removed data associated with the float being processed")

    def test_returns_correct_shape(self):
        """
        See if the returned data has the expected shape
        :return: Nothing
        """
        print("Testing that get_data returns data with the expected shape")

        wmo_boxes = np.array([3505, 1, 0, 0])
        data = get_data(wmo_boxes, 1, self.config, self.float_name)

        self.assertTrue(data['sal'].shape == data['ptmp'].shape == data['temp'].shape,
                        "Ocean characteristic data should be the same shape")
        self.assertTrue(data['long'].shape == data['lat'].shape == data['dates'].shape,
                        "Spatial/temporal data should be the same shape")
        self.assertTrue(data['sal'].shape[1] == data['long'].shape[1],
                        "Should be a profile for every location")


if __name__ == '__main__':
    unittest.main()
