"""
-----Interpolate Climatology-----

Written by: Breck Owens
When: xx/07/2005
Converted to python by: Edward Small
When: 17/01/2020

Routine to interpolate climi=atological salinity and pressure data onto float's
potential temperature.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np


def interp_climatology(grid_sal, grid_theta, grid_pres, float_sal, float_theta, float_pres):
    """
    interpolate historical salinity and pressure data on the float theta
    :param grid_sal: historical salinity
    :param grid_theta: historical potential temperature
    :param grid_pres: historical pressure
    :param float_sal: currentfloat salinity
    :param float_theta: current float potential temperature
    :param float_pres: current float pressure
    :return: matrices [number of floats x number of historical profiles] of
    interpolated salinity and interpolated pressure on float theta surface
    """

    # get the shape of the data inputs (float and climatology)
    grid_level = grid_sal.shape[0]
    grid_station = grid_sal.shape[1]
    float_shape = float_sal.shape[0]

    # check that the climatology data has no infinite (bad) values in the middle
    # of the profiles.
    max_level = 0
    for n in range(0, grid_station):
        # find where data is not missing
        good_sal = np.isfinite(grid_sal[:, n])
        good_theta = np.isfinite(grid_theta[:, n])
        good_pres = np.isfinite(grid_pres[:, n])

        # find indices of good data
        good_data_index = (good_sal == good_theta)
        good_data_index = np.argwhere(good_data_index == good_pres)

        # now find the max level
        good_data_len = good_data_index.__len__()
        if good_data_len != 0:
            for m in range(0, good_data_len):
                grid_sal[m, n] = grid_sal[good_data_index[m], n]
                grid_theta[m, n] = grid_theta[good_data_index[m], n]
                grid_pres[m, n] = grid_pres[good_data_index[m], n]
            max_level = np.maximum(max_level, good_data_len)

    # Truncate the number of levels to the maximum level that has
    # available data
    if max_level > 0:
        grid_sal = grid_sal[:max_level, :]
        grid_theta = grid_theta[:max_level, :]
        grid_pres = grid_pres[:max_level, :]
    else:
        raise ValueError("No good climatological data has been found")
