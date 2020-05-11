import pandas as pd
import scipy.io as scipy
import matplotlib.pyplot as plt
from ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot import trajectory_plot


def plot_diagnostics(float_dir, float_name, config):
    # load the data into dataframes
    grid_data_loc = config['FLOAT_MAPPED_DIRECTORY'] + config['FLOAT_MAPPED_PREFIX'] + \
                    float_name + config['FLOAT_MAPPED_POSTFIX']
    float_data_loc = config['FLOAT_SOURCE_DIRECTORY'] + float_dir + \
                     float_name + config['FLOAT_SOURCE_POSTFIX']

    grid_data = scipy.loadmat(grid_data_loc)
    float_data = scipy.loadmat(float_data_loc)

    # create trajectory plot ------------------------------
    grid_lat = grid_data['selected_hist'][:, 1].flatten()
    grid_long = grid_data['selected_hist'][:, 0].flatten()

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

    trajectory_plot(0, 0, floats, grid, float_name)

    plt.show()
