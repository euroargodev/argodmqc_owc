""" Utilities Functions for plots
"""
import pandas as pd


def create_dataframe(grid_data, float_data):
    """ Return data frame of location data for trajectory plotting

        Parameters
        ----------
        grid_data: data from the historical data set
        float_data: data from the profile

        Returns
        -------
        Dataframes containing float and grid location data
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
