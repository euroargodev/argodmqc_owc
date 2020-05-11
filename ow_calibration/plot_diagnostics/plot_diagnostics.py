"""
-----plot diagnostics----

Written by: Breck Owens
When: 14/06/2011
 Converted to python by: Edward Small
When: 09/05/2020

Master function for loading and plotting all of the analysed data

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""


import scipy.io as scipy
import matplotlib.pyplot as plt
from ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot import trajectory_plot
from ow_calibration.plot_diagnostics.trajectory_plot.create_dataframe import create_dataframe


def plot_diagnostics(float_dir, float_name, config):
    # load the data into dataframes
    grid_data_loc = config['FLOAT_MAPPED_DIRECTORY'] + config['FLOAT_MAPPED_PREFIX'] + \
                    float_name + config['FLOAT_MAPPED_POSTFIX']
    float_data_loc = config['FLOAT_SOURCE_DIRECTORY'] + float_dir + \
                     float_name + config['FLOAT_SOURCE_POSTFIX']

    grid_data = scipy.loadmat(grid_data_loc)
    float_data = scipy.loadmat(float_data_loc)

    # create trajectory plot ------------------------------
    grid, floats = create_dataframe(grid_data, float_data)

    trajectory_plot(0, 0, floats, grid, float_name)

    plt.show()
