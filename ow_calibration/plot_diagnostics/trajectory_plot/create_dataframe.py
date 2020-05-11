"""
-----Create dataframe-----

Written by: Edward Small
When: 11/05/2020

Refactored into it's own function. Creates a nice data frame to read
location data from, which is useful for plotting the float trajectory

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""


import pandas as pd


def create_dataframe(grid_data, float_data):
    """
    Return data frame of location data for trajectory plotting
    :param grid_data: data from the historical data set
    :param float_data: data from the profile
    :return: dataframes containing float and grid location data
    """
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

    return grid, floats
