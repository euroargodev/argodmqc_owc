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
from scipy.io import loadmat
from pyowc.plot import plots

#pylint: disable=bare-except
#pylint: disable=unused-argument
#pylint: disable=too-many-locals
class MyTestCase(unittest.TestCase):
    """
    Test cases for sal_var_plot function
    """

    @patch("pyowc.plot.plots.plt.show")
    def test_plot_runs(self, mockshow):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that salinity variance plot throws no errors")

        grid_data = loadmat("../../data/test_data/float_mapped_test/map_3901960.mat")
        float_data = loadmat("../../data/float_source/3901960.mat")
        cal_data = loadmat("../../data/test_data/float_calib_test/cal_3901960.mat")

        sal = float_data['SAL']
        theta = float_data['PTMP']
        pres = float_data['PRES']
        grid_sal = grid_data['la_mapped_sal']
        grid_ptmp = grid_data['la_ptmp']
        grid_errors = grid_data['la_mapsalerrors']
        cal_sal = cal_data['cal_SAL']
        cal_sal_errors = cal_data['cal_SAL_err']

        no_boundaries = [0, 0, 0, 0, 0.5]
        low_bound = [250, 0, 250, 0, 0.5]
        up_bound = [0, 1, 0, 1, 0.5]

        profile_no = float_data['PROFILE_NO']

        # Check various types run

        self.assertEqual(plots.sal_var_plot(2, copy.deepcopy(sal), copy.deepcopy(pres),
                                      copy.deepcopy(theta), copy.deepcopy(grid_sal),
                                      copy.deepcopy(grid_errors), copy.deepcopy(grid_ptmp),
                                      copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_errors),
                                      no_boundaries,
                                      profile_no, "3902960"), None)
        self.assertEqual(plots.sal_var_plot(2, copy.deepcopy(sal), copy.deepcopy(pres),
                                      copy.deepcopy(theta), copy.deepcopy(grid_sal),
                                      copy.deepcopy(grid_errors), copy.deepcopy(grid_ptmp),
                                      copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_errors),
                                      low_bound,
                                      profile_no, "3902960"), None)
        self.assertEqual(plots.sal_var_plot(2, copy.deepcopy(sal), copy.deepcopy(pres),
                                      copy.deepcopy(theta), copy.deepcopy(grid_sal),
                                      copy.deepcopy(grid_errors), copy.deepcopy(grid_ptmp),
                                      copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_errors),
                                      up_bound,
                                      profile_no, "3902960"), None)


if __name__ == '__main__':
    unittest.main()
