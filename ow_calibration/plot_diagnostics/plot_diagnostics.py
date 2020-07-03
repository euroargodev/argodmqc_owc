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
    """
    run the plotting procedures
    :param float_dir: location of float source
    :param float_name: name of the float source
    :param config: user configuration dictionary
    :return: nothing, but will save the plots as PDFs
    """

    # Get float data, mapped data, and calibrated data --------

    mapped_data = config['FLOAT_MAPPED_DIRECTORY'] + config['FLOAT_MAPPED_PREFIX'] + \
                    float_name + config['FLOAT_MAPPED_POSTFIX']
    float_source_data = config['FLOAT_SOURCE_DIRECTORY'] + float_dir + \
                     float_name + config['FLOAT_SOURCE_POSTFIX']

    mapped_data = scipy.loadmat(mapped_data)
    float_source_data = scipy.loadmat(float_source_data)


    mapped_loc = mapped_data['selected_hist']

    float_lat = float_source_data['LAT'].flatten()
    float_long = float_source_data['LONG'].flatten()
    float_no = float_source_data['PROFILE_NO'].flatten()

    # create trajectory plot ------------------------------

    trajectory_plot(mapped_loc, float_long, float_lat, bathy=True, cmap='block')
