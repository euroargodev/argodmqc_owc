import unittest
import pandas as pd
import scipy.io as scipy
from ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot import trajectory_plot


class MyTestCase(unittest.TestCase):
    def test_something(self):
        grid_data = scipy.loadmat("data/float_mapped/map_3901960.mat")
        grid_lat = grid_data['selected_hist'][:, 1].flatten()
        grid_long = grid_data['selected_hist'][:, 0].flatten()

        float_data = scipy.loadmat("data/float_source/3901960.mat")
        float_lat = float_data['LAT'].flatten()
        float_long = float_data['LONG'].flatten()
        float_no = float_data['PROFILE_NO'].flatten()

        grid = pd.DataFrame(
            {
                'Latitude': grid_lat,
                'Longitude': grid_long})

        floats = pd.DataFrame(
            {'number': float_no,
             'Latitude': float_lat,
             'Longitude': float_long})

        try:
            trajectory_plot(1, 1, floats, grid, "3901960")

        except:
            self.fail("Trajectory plotting routine failed unexpectedly")


if __name__ == '__main__':
    unittest.main()
