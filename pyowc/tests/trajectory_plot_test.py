"""
-----Trajectory Plot Test File-----

Written by: Edward Small
When: 11/05/2019

Contains unit tests to check the functionality of the `trajectory_plot` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""


import unittest
from scipy.io import loadmat
from pyowc.plot import plots, utils


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

        config = {'CONFIG_DIRECTORY': "../../data/constants/"}

        grid_data = loadmat("../../data/float_mapped/map_3901960.mat")
        float_data = loadmat("../../data/float_source/3901960.mat")

        grid, floats = utils.create_dataframe(grid_data, float_data)

        try:
            plots.trajectory_plot(1, 1, floats, grid, "3901960", config)

        except:
            self.fail("Trajectory plotting routine failed unexpectedly")


if __name__ == '__main__':
    unittest.main()
