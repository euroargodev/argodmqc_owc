"""
-----Get Historical Region Location Data-----

Written by: Breck Owens
When: xx/12/2006
Converted to python by: Edward Small
When: 05/12/2019

Function returns the needed spatial and temporal data from each of the historical
data points inside the WMO boxes it is given, excluding the current float being processed.
The WMO boxes are passed in the form of [WMO name, use ctd data, use bottle data, use argo data],
and the configuration is also passed in, which is used to locate the different data locally. The
configuration can be changed in load_configuration.py

N.B. Change to code on xx/xx/2015 - Only load longitude, latitude, and dates

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
from ow_calibration.change_dates.change_dates import change_dates
from ow_calibration.get_region.data_functions.get_data import get_data
from ow_calibration.get_region.data_functions.wrap_longitude import wrap_longitude


def get_region_hist_locations(pa_wmo_numbers, pa_float_name, config):
    """
    Uses the WMO boxes and to return all of the historical data in the given area,
    excluding the float that is currently being analysed.
    :param pa_wmo_numbers: 2D array containing the name of the WMO boxes that cover the area
    of interest, and flags for whether we want to use argo, bottle, and/or CTD data
    :param pa_float_name: string of the name of the float currently being processed
    :param config: Dictionary containing configuration settings. Used to find locations of folders
    and file containing data
    :return: the latitude, longitude and age of each data point we want
    """

    # set up matrices to hold data
    grid_lat = []
    grid_long = []
    grid_dates = []

    # go through each of the WMO boxes
    for wmo_box in pa_wmo_numbers:

        # go through each of the columns denoting whether we should use CTD, bottle, and/or argo
        for data_type in range(1, 4):

            # get the data
            data = get_data(wmo_box, data_type, config, pa_float_name)

            # if we have data, combine it with the other data then reset it
            if data:
                grid_lat = np.concatenate([grid_lat, data['lat'][0]])
                grid_long = np.concatenate([grid_long, data['long'][0]])
                grid_dates = np.concatenate([grid_dates, data['dates'][0]])
                data = []

    if grid_lat.__len__() == 0:

        grid_lat = 999
        grid_long = 999
        grid_dates = 'NaN'

    else:

        grid_long = wrap_longitude(grid_long)
        # decimalise dates
        grid_dates = change_dates(grid_dates)

    return grid_lat, grid_long, grid_dates
