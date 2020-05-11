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
from ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot import trajectory_plot
from ow_calibration.plot_diagnostics.trajectory_plot.create_dataframe import create_dataframe


# pylint: disable=bare-except
class MyTestCase(unittest.TestCase):
    """
    Test cases for trajectory_plot function
    """

    def test_plot_runs(self):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that trajectory plot throws no errors")

        grid_data = scipy.loadmat("data/float_mapped/map_3901960.mat")
        float_data = scipy.loadmat("data/float_source/3901960.mat")

        grid, floats = create_dataframe(grid_data, float_data)


        try:
            trajectory_plot(1, 1, floats, grid, "3901960")

        except:
            self.fail("Trajectory plotting routine failed unexpectedly")


if __name__ == '__main__':
    unittest.main()
