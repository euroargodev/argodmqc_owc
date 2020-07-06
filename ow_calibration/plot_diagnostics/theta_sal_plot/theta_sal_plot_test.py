"""
-----Theta Salinity Curve Test File-----

Written by: Edward Small
When: 11/05/2019
Contains unit tests to check the functionality of the `trajectory_plot` function
To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import unittest
from unittest.mock import patch
import scipy.io as scipy
import numpy as np
from ow_calibration.find_10thetas.find_10thetas import find_10thetas
from ow_calibration.plot_diagnostics.theta_sal_plot.theta_sal_plot import theta_sal_plot


# pylint: disable=bare-except
# pylint: disable=unused-argument
class MyTestCase(unittest.TestCase):
    """
    Test cases for theta_sal_plot function
    """

    @patch("ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot.plt.show")
    def test_plot_runs(self, mockshow):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that theta salinity plot throws no errors")

        grid_data = scipy.loadmat("data/test_data/float_mapped_test/map_3901960.mat")
        float_data = scipy.loadmat("data/float_source/3901960.mat")

        sal = np.array(float_data['SAL'])
        theta = np.array(float_data['PTMP'])
        pres = float_data['PRES']
        grid_sal = grid_data['la_mapped_sal']
        grid_ptmp = grid_data['la_ptmp']
        grid_errors = grid_data['la_mapsalerrors']

        thetas = find_10thetas(copy.deepcopy(sal),
                               copy.deepcopy(theta),
                               copy.deepcopy(pres),
                               copy.deepcopy(grid_ptmp),
                               0, 0,
                               0, 0,
                               0.5)

        index = thetas[2]

        # Check various types run

        self.assertEqual(theta_sal_plot(sal.transpose(), theta.transpose(),
                                        grid_sal, grid_ptmp, grid_errors, index), None)


if __name__ == '__main__':
    unittest.main()
