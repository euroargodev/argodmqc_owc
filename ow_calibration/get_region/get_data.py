"""
-----Get Data File-----

Written by: Edward Small
When: 30/01/2020

Uses the WMO box numbers to fetch the relevant data for the region.

This is a refactor between get_region_data and get_region_hist_locations to avoid
duplicate code

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
import scipy.io as scipy


def get_data(wmo_box, data_type, config, pa_float_name):
    """
    Gets all the data we highlight we want from regions specified by the WMO boxes
    :param pa_wmo_numbers: 2D array containing the name of the WMO boxes that cover the area
    of interest, and flags for whether we want to use argo, bottle, and/or CTD data
    :param data_type: Which type of data we are checking for
    :param config: Dictionary containing configuration settings. Used to find locations of folders
    and file containing data
    :param pa_float_name: String containing the name of the float being profiled
    :return: array containing all (ctd/bottle/argo) data in a wmo box
    """

    data = []
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

        data['lat'] = [np.delete(data['lat'], not_use)]
        data['long'] = [np.delete(data['long'], not_use)]
        data['dates'] = [np.delete(data['dates'], not_use)]
        data['sal'] = [np.delete(data['sal'], not_use, axis=1)]
        data['ptmp'] = [np.delete(data['ptmp'], not_use, axis=1)]
        data['pres'] = [np.delete(data['pres'], not_use, axis=1)]

    return data
