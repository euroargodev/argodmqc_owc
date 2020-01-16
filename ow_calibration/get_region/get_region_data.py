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

    max_pres = pa_float_pres.max() + config["MAP_P_DELTA"]
    print(max_pres)

