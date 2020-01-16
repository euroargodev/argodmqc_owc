"""
-----Get Historical Data Test File-----

Written by: Edward Small
When: 05/12/2019

Contains unit tests to check the functionality of the `get_hist_region_locations` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from ow_calibration.get_region.get_region_data import get_region_data


class MyTestCase(unittest.TestCase):
    def test_something(self):
        wmo = np.array([[3505, 1, 0, 0],
                        [3506, 1, 0, 0]])
        float_name = "3901960"
        index = np.array([1, 5, 6, 12, 14])
        pres = np.array([3, 5, 15.1, 25.1, 36])
        config = {'HISTORICAL_DIRECTORY': "data/climatology",
                           'HISTORICAL_CTD_PREFIX': "/historical_ctd/ctd_",
                           'HISTORICAL_BOTTLE_PREFIX': "/historical_bot/bot_",
                           'HISTORICAL_ARGO_PREFIX': "/historical_argo/argo_",
                           'MAP_P_DELTA': 250}

        get_region_data(wmo, float_name, config, index, pres)


if __name__ == '__main__':
    unittest.main()
