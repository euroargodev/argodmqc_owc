"""
-----Get Historical Data Test File-----

Written by: Edward Small
When: 05/12/2019

Contains unit tests to check the functionality of the `get_hist_region_locations` function.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""
import unittest
import numpy as np
from ow_calibration.get_region.get_region_data import get_region_data


class MyTestCase(unittest.TestCase):
    """
    Test cases for get_region_data function
    """

    def setUp(self):
        self.wmo_boxes = np.array([[3505, 1, 1, 1],
                                   [3506, 1, 1, 1]])
        self.float_name = "3901960"
        self.index = np.array([0, 4, 5, 11, 13, 15, 20, 21, 42])
        self.pres = np.array([3, 5, 15.1, 25.1, 36, 40, 45, 46, 500000])
        self.config = {'HISTORICAL_DIRECTORY': "data/climatology",
                       'HISTORICAL_CTD_PREFIX': "/historical_ctd/ctd_",
                       'HISTORICAL_BOTTLE_PREFIX': "/historical_bot/bot_",
                       'HISTORICAL_ARGO_PREFIX': "/historical_argo/argo_",
                       'MAP_P_DELTA': 250}

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

        self.assertTrue(test[0].shape == test[1].shape == test[2].shape,
                        "salinity, pressure, and potential temperature "
                        "arrays should be same shape")
        self.assertTrue(test[3].shape == test[4].shape == test[5].shape,
                        "longitude, latitude, and date arrays should be the same shape")
        self.assertTrue(test[0].shape[1] == test[1].shape[1] == test[2].shape[1] ==
                        test[3].shape[0] == test[4].shape[0] == test[5].shape[0] ==
                        self.index.__len__(),
                        "Should get the same number of casts as we have asked for")

    def test_gets_different_data(self):
        """
        Check that we can fetch different combinations of data
        :return: Nothing
        """
        print("Testing that get_region_data can return different data types")

        test_ctd = get_region_data(np.array([[3505, 1, 0, 0]]), self.float_name, self.config,
                                   self.index, self.pres)
        test_bot = get_region_data(np.array([[3505, 0, 1, 0]]), self.float_name, self.config,
                                   self.index, self.pres)
        test_argo = get_region_data(np.array([[3505, 0, 0, 1]]), self.float_name, self.config,
                                    self.index, self.pres)

        self.assertTrue(test_ctd[0].shape[1] != test_argo[0].shape[1],
                        "Should get a different data set, if we have specified it")
        self.assertTrue(test_bot[0].shape[1] != test_argo[0].shape[1],
                        "Should get a different data set, if we have specified it")

    def test_should_only_get_specific_indices(self):
        """
        If we only give n indices, we should get n data points back
        :return: Nothing
        """
        print("Testing that get_region returns correct amount of data")

        test_many = get_region_data(self.wmo_boxes, self.float_name, self.config,
                                   self.index, self.pres)

        self.assertTrue(test_many[0].shape[1] == self.index.__len__())

        test_one = get_region_data(self.wmo_boxes, self.float_name, self.config,
                                   [50], self.pres)

        self.assertTrue(test_one[0].shape[1] == 1)

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



if __name__ == '__main__':
    unittest.main()
