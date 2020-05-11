import unittest
import pandas as pd
import scipy.io as scipy
from ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot import trajectory_plot
from ow_calibration.plot_diagnostics.trajectory_plot.create_dataframe import create_dataframe


class MyTestCase(unittest.TestCase):

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
