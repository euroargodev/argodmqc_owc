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
import scipy.io as scipy
from ow_calibration.change_dates.change_dates import change_dates


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

                # if data dimensions don't match, transpose to avoid indexing issues
                if (data['lat'].size == data['pres'].size and
                        data['lat'].shape[1] == data['pres'].shape[1]):
                    data['pres'] = data['pres'].T
                    data['ptmp'] = data['ptmp'].T
                    data['sal'] = data['sal'].T
                    data['temp'] = data['temp'].T

            if wmo_box[data_type] == 1 and data_type == 3:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] +
                                     config['HISTORICAL_ARGO_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')

                # remove the argo float being analysed from the data
                not_use = []
                for i in range(0, data['lat'][0].__len__()):
                    if str(data['source'][0][i]).find(pa_float_name) != -1:
                        not_use.append(i)

                print(not_use.__len__())
                print(data['lat'].shape)
                data['lat'] = np.delete(data['lat'], not_use)
                data['long'] = np.delete(data['long'], not_use)
                data['dates'] = np.delete(data['dates'], not_use)
                data['sal'] = np.delete(data['sal'], not_use, axis=1)
                data['ptmp'] = np.delete(data['ptmp'], not_use, axis=1)
                data['pres'] = np.delete(data['pres'], not_use, axis=1)
                print(data['lat'].shape)

