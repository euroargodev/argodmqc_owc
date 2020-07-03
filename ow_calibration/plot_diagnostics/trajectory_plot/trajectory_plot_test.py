"""
-----Trajectory Plot Test File-----

Written by: Edward Small
When: 11/05/2019

Contains unit tests to check the functionality of the `trajectory_plot` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import scipy.io as scipy
from unittest.mock import patch
from ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot import trajectory_plot


# pylint: disable=bare-except
# pyline: disable=unused-argument
class MyTestCase(unittest.TestCase):
    """
    Test cases for trajectory_plot function
    """

    @patch("ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot.plt.show")
    def test_plot_runs(self, mockshow):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that trajectory plot throws no errors")

        float_name = "TESTING NAME"

        grid_data = scipy.loadmat("data/float_mapped/map_3901960.mat")
        float_data = scipy.loadmat("data/float_source/3901960.mat")

        mapped_loc = grid_data['selected_hist']

        float_lat = float_data['LAT'].flatten()
        float_long = float_data['LONG'].flatten()

        # Check various types run

        self.assertEqual(trajectory_plot(float_name, mapped_loc, float_long, float_lat),
                         None)
        self.assertEqual(trajectory_plot(float_name, mapped_loc, float_long, float_lat, bathy=True),
                         None)
        self.assertEqual(trajectory_plot(float_name, mapped_loc, float_long, float_lat,
                                         bathy=True, style='shade'),
                         None)


if __name__ == '__main__':
    unittest.main()
