"""
-----Get Historical Data-----

Written by: Breck Owens
When: xx/12/2006
Converted to python by: Edward Small
When: 16/01/2020

Function returns the needed spatial, temporal, and ocean characteristic (salinity, temp, pressure)
data from each of the historical data points inside the WMO boxes it is given, excluding the
current float being processed. The WMO boxes are passed in the form of [WMO name, use ctd data,
use bottle data, use argo data], and the configuration is also passed in, which is used to locate
the different data locally. The configuration can be changed in load_configuration.py

Isabelle Gaboury, 26 Sep. 2017: Added check on the dimensions of bottle data.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
from ow_calibration.change_dates.change_dates import change_dates
from ow_calibration.get_region.data_functions.get_data import get_data
from ow_calibration.get_region.data_functions.wrap_longitude import wrap_longitude

# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
def get_region_data(pa_wmo_numbers, pa_float_name, config, index, pa_float_pres):
    """
    Get the historical pressure, salinity, and temperature of selected casts
    :param pa_wmo_numbers: 2D array containing the name of the WMO boxes that cover the area
    of interest, and flags for whether we want to use argo, bottle, and/or CTD data
    :param pa_float_name: string of the name of the float currently being processed
    :param config: Dictionary containing configuration settings. Used to find locations of folders
    and file containing data
    :param index: array of indices of selected historical casts
    :param pa_float_pres: array of pressures for the float being processed
    :return: The salinity, potential temperature, pressure, latitude, longitude, and age of each
    historical cast selected to use
    """

    # maximum depth to retrieve data from (deepest float measurement + MAP_P_DELTA
    max_pres = pa_float_pres.max() + config["MAP_P_DELTA"]

    # set up empty arrays to hold the data to return
    grid_sal = []
    grid_ptmp = []
    grid_pres = []
    grid_lat = []
    grid_long = []
    grid_dates = []
    data = []

    # set up current maximum depth and number of columns
    max_depth = 0
    how_many_cols = 0

    # set up variable to save beginning index for each set of data
    starting_index = 0

    # go through each of the WMO boxes
    for wmo_box in pa_wmo_numbers:

        # go through each of the columns denoting whether we should use CTD, bottle, and/or argo
        for data_type in range(1, 4):

            # get the data
            data = get_data(wmo_box, data_type, config, pa_float_name)

            if data:
                # Sometimes the data comes in wrapped in as a 3d array, so convert to 2d
                if data['pres'].__len__() == 1:
                    data['pres'] = data['pres'][0]
                    data['sal'] = data['sal'][0]
                    data['ptmp'] = data['ptmp'][0]
                    data['lat'] = data['lat'][0].reshape(-1,1)
                    data['long'] = data['long'][0].reshape(-1,1)
                    data['dates'] = data['dates'][0].reshape(-1,1)

                #  check the index of each station to see if it should be loaded
                data_length = data['lat'][0].__len__()
                data_indices = np.arange(0, data_length) + starting_index

                # remember location of last entry
                starting_index = starting_index + data_length

                # load each station
                for i in range(0, data_length):

                    good_indices = np.argwhere(index == data_indices[i])

                    if good_indices.__len__() > 0:
                        # only use non-NaN values
                        not_nan = np.argwhere(np.isnan(data['pres'][:, i]) == 0)

                        # get the non-NaN values
                        pres = data['pres'][not_nan, i]
                        sal = data['sal'][not_nan, i]
                        ptmp = data['ptmp'][not_nan, i]

                        # remove values where pressure exceeds the maximum we want
                        too_deep = np.argwhere(pres > max_pres)
                        pres = np.delete(pres, too_deep[:, 0])
                        sal = np.delete(sal, too_deep[:, 0])
                        ptmp = np.delete(ptmp, too_deep[:, 0])
                        new_depth = pres.__len__()
                        how_many_rows = np.max([new_depth, max_depth])

                        # if the new data we are adding is longer than our columns, we need to fill in
                        # NaNs in the other columns
                        if new_depth > max_depth != 0:
                            grid_pres = np.append(grid_pres, np.ones(
                                (how_many_cols, new_depth - max_depth)) * np.nan, axis=1
                                                  ).reshape((how_many_cols, how_many_rows))
                            grid_ptmp = np.append(grid_ptmp, np.ones(
                                (how_many_cols, new_depth - max_depth)) * np.nan, axis=1
                                                  ).reshape((how_many_cols, how_many_rows))
                            grid_sal = np.append(grid_sal, np.ones(
                                (how_many_cols, new_depth - max_depth)) * np.nan, axis=1
                                                 ).reshape((how_many_cols, how_many_rows))

                        # if the new data we are adding is shorter than our columns, then we need to
                        # fill in the rest with NaNs so it's the same length
                        elif new_depth < max_depth:
                            pres = np.append(pres, np.ones((max_depth - new_depth, 1)) * np.nan)
                            ptmp = np.append(ptmp, np.ones((max_depth - new_depth, 1)) * np.nan)
                            sal = np.append(sal, np.ones((max_depth - new_depth, 1)) * np.nan)

                        # if we don't have any data saved yet, create the grid matrix with the
                        # first data set
                        if grid_pres.__len__() == 0:

                            grid_pres = pres.reshape((1, pres.__len__()))
                            grid_ptmp = ptmp.reshape((1, pres.__len__()))
                            grid_sal = sal.reshape((1, pres.__len__()))

                        # if we already have data saved, add the new data to the saved data
                        else:

                            grid_pres = np.append(grid_pres, pres).reshape(
                                how_many_cols + 1, how_many_rows)
                            grid_ptmp = np.append(grid_ptmp, ptmp).reshape(
                                how_many_cols + 1, how_many_rows)
                            grid_sal = np.append(grid_sal, sal).reshape(
                                how_many_cols + 1, how_many_rows)

                        # save the latitude, longitude, and date of the new data
                        grid_lat = np.append(grid_lat, data['lat'][0, i])
                        grid_long = np.append(grid_long, data['long'][0, i])
                        grid_dates = np.append(grid_dates, data['dates'][0, i])

                        # readjust our values so we know what column to add the new data to,
                        # and what shape we should expect the data to be
                        max_depth = grid_pres.shape[1]
                        how_many_cols = grid_pres.shape[0]

    # convert longitude to 0 to 360 degrees
    try:
        grid_long = wrap_longitude(grid_long)

        # make sure salinity, pressure, and potential temperature data have all the same NaNs
        sal_nans = np.argwhere(np.isnan(grid_sal))
        for nan in sal_nans:
            grid_pres[nan[0], nan[1]] = np.nan
            grid_ptmp[nan[0], nan[1]] = np.nan

        pres_nans = np.argwhere(np.isnan(grid_pres))
        for nan in pres_nans:
            grid_sal[nan[0], nan[1]] = np.nan
            grid_ptmp[nan[0], nan[1]] = np.nan

        ptmp_nans = np.argwhere(np.isnan(grid_ptmp))
        for nan in ptmp_nans:
            grid_sal[nan[0], nan[1]] = np.nan
            grid_pres[nan[0], nan[1]] = np.nan

        grid_dates = change_dates(grid_dates)

        # transpose data
        grid_sal = grid_sal.T
        grid_pres = grid_pres.T
        grid_ptmp = grid_ptmp.T

    except:
        raise Exception("NO DATA FOUND")

    return grid_sal, grid_ptmp, grid_pres, grid_lat, grid_long, grid_dates
