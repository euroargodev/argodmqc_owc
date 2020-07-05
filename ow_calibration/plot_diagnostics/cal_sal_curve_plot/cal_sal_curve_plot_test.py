"""
-----Theta Salinity Curve Test File-----

Written by: Edward Small
When: 11/05/2019

Contains unit tests to check the functionality of the `cal_sal_curve_plot` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import unittest
from unittest.mock import patch
import scipy.io as scipy
import numpy as np
from ow_calibration.find_10thetas.find_10thetas import find_10thetas
from ow_calibration.plot_diagnostics.cal_sal_curve_plot.cal_sal_curve_plot import cal_sal_curve_plot


# pylint: disable=bare-except
# pylint: disable=unused-argument
class MyTestCase(unittest.TestCase):
    """
    Test cases for cal_sal_curve_plot function
    """

    #@patch("ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot.plt.show")
    #def test_plot_runs(self, mockshow):
    def test_plot_runs(self):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that calibrated salinity curve plot throws no errors")

        # get the data we need

        cal_data = scipy.loadmat("data/test_data/float_calib_test/cal_3901960.mat")
        float_data = scipy.loadmat("data/float_source/3901960.mat")

        sal = float_data['SAL']
        cal_sal = cal_data['cal_SAL']
        cal_sal_err = cal_data['cal_SAL_err']
        sta_sal = cal_data['sta_SAL']
        sta_sal_err = cal_data['sta_SAL_err']
        sta_mean = cal_data['sta_mean']
        pcond_factor = np.array(cal_data['pcond_factor'])
        pcond_factor_err = np.array(cal_data['pcond_factor_err'])
        float_name = "3901960"
        profile_no = float_data['PROFILE_NO']

        cal_sal_curve_plot(sal, cal_sal, cal_sal_err, sta_sal, sta_sal_err, sta_mean,
                           pcond_factor, pcond_factor_err, profile_no, float_name)


if __name__ == '__main__':
    unittest.main()