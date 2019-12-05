"""
-----Get Historical Region Data-----

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
import scipy.io as scipy
from ow_calibration.change_dates.change_dates import change_dates


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
    not_use = []
    data = []

    # go through each of the WMO boxes
    for wmo_box in pa_wmo_numbers:

        # go through each of the columns denoting whether we should use CTD, bottle, and/or argo
        for data_type in range(1, 4):

            # check if we should use this data. If so, get the data
            if wmo_box[data_type] == 1 and data_type == 1:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] +
                                     config['HISTORICAL_CTD_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')

            if wmo_box[data_type] == 1 and data_type == 2:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] +
                                     config['HISTORICAL_BOTTLE_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')

            if wmo_box[data_type] == 1 and data_type == 3:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] +
                                     config['HISTORICAL_ARGO_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')

                # remove the argo float being analysed from the data
                for i in range(0, data['lat'][0].__len__()):
                    if str(data['source'][0][i]).find(pa_float_name) != -1:
                        not_use.append(i)

                data['lat'] = [np.delete(data['lat'], not_use)]
                data['long'] = [np.delete(data['long'], not_use)]
                data['dates'] = [np.delete(data['dates'], not_use)]

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
        # convert longitude to 0 to 360 degrees
        neg_long = np.argwhere(grid_long < 0)
        grid_long[neg_long] = grid_long[neg_long] + 360

        # if we have data close to upper boundary (360), then wrap some of the data round
        # so it appears on the map
        top_long = np.argwhere(grid_long >= 320)
        if top_long.__len__() != 0:
            bottom_long = np.argwhere(grid_long <= 40)
            grid_long[bottom_long] = 360 + grid_long[bottom_long]

        # decimalise dates
        grid_dates = change_dates(grid_dates)

    return grid_lat, grid_long, grid_dates
