"""
-----profile plot test file-----

Written by: Edward Small
When: 11/05/2019

Contains unit tests to check the functionality of the `t_s_profile_plot` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import unittest
from unittest.mock import patch
import scipy.io as scipy
import numpy as np
from ow_calibration.plot_diagnostics.t_s_profile_plot.t_s_profile_plot import t_s_profile_plot
from ow_calibration.find_10thetas.find_10thetas import find_10thetas


# pylint: disable=bare-except
# pylint: disable=unused-argument
# pylint: disable=too-many-locals
class MyTestCase(unittest.TestCase):
    """
    Test cases for t_s_plot function
    """

    @patch("ow_calibration.plot_diagnostics.t_s_profile_plot.t_s_profile_plot.plt.show")
    def test_plot_runs(self, mockshow):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """
        print("Check t_s_profile_plot runs")
        float_data = scipy.loadmat("data/float_source/3901960.mat")
        grid_data = scipy.loadmat("data/test_data/float_mapped_test/map_3901960.mat")

        sal = np.array(float_data['SAL'])
        ptmp = np.array(float_data['PTMP'])
        pres = float_data['PRES']
        grid_ptmp = grid_data['la_ptmp']

        thetas = find_10thetas(copy.deepcopy(sal),
                               copy.deepcopy(ptmp),
                               copy.deepcopy(pres),
                               copy.deepcopy(grid_ptmp),
                               0, 0,
                               0, 0,
                               0.5)

        sal_var = thetas[3]
        theta_levels = thetas[4]
        tlevels = thetas[0]
        plevels = thetas[1]

        self.assertEqual(t_s_profile_plot(sal, ptmp, pres, sal_var,
                                          theta_levels, tlevels, plevels, "3901960"), None)


if __name__ == '__main__':
    unittest.main()
