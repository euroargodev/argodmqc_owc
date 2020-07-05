"""
-----Salinity variance on theta levels plot test file-----

Written by: Edward Small
When: 11/05/2019

Contains unit tests to check the functionality of the `sal_var_plot` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import unittest
from unittest.mock import patch
import scipy.io as scipy
import numpy as np
from ow_calibration.find_10thetas.find_10thetas import find_10thetas
from ow_calibration.plot_diagnostics.sal_var_plot.sal_var_plot import sal_var_plot


# pylint: disable=bare-except
# pylint: disable=unused-argument
class MyTestCase(unittest.TestCase):
    """
    Test cases for sal_var_plot function
    """

    #@patch("ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot.plt.show")
    #def test_plot_runs(self, mockshow):
    def test_plot_runs(self):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that salinity variance plot throws no errors")

        grid_data = scipy.loadmat("data/test_data/float_mapped_test/map_3901960.mat")
        float_data = scipy.loadmat("data/float_source/3901960.mat")
        cal_data = scipy.loadmat("data/test_data/float_calib_test/cal_3901960.mat")

        sal = np.array(float_data['SAL'])
        theta = np.array(float_data['PTMP'])
        pres = float_data['PRES']
        grid_sal = grid_data['la_mapped_sal']
        grid_ptmp = grid_data['la_ptmp']
        grid_errors = grid_data['la_mapsalerrors']
        cal_sal = cal_data['cal_SAL']
        cal_sal_errors = cal_data['cal_SAL_err']

        boundaries = [0, 0, 0, 0, 0.5]

        # Check various types run

        sal_var_plot(2, sal, pres, theta, grid_sal, grid_errors,
        grid_ptmp, cal_sal, cal_sal_errors, boundaries)


if __name__ == '__main__':
    unittest.main()
