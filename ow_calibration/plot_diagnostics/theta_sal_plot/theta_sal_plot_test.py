"""
-----Theta Salinity Curve Test File-----

Written by: Edward Small
When: 11/05/2019
Contains unit tests to check the functionality of the `trajectory_plot` function
To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import scipy.io as scipy
from unittest.mock import patch
from ow_calibration.plot_diagnostics.theta_sal_plot.theta_sal_plot import theta_sal_plot


# pylint: disable=bare-except
# pyline: disable=unused-argument
class MyTestCase(unittest.TestCase):
    """
    Test cases for theta_sal_plot function
    """

    #@patch("ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot.plt.show")
    # def test_plot_runs(self, mockshow):
    def test_plot_runs(self):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that theta salinity plot throws no errors")

        grid_data = scipy.loadmat("data/test_data/float_mapped_test/map_3901960.mat")
        float_data = scipy.loadmat("data/float_source/3901960.mat")

        sal = float_data['SAL'].transpose()
        theta = float_data['PTMP'].transpose()
        grid_sal = grid_data['la_mapped_sal']
        grid_ptmp = grid_data['la_ptmp']
        grid_errors = grid_data['la_mapsalerrors']

        # Check various types run

        theta_sal_plot(sal, theta, grid_sal, grid_ptmp, grid_errors)


if __name__ == '__main__':
    unittest.main()
