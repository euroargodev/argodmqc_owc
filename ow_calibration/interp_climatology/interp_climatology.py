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
    :param float_sal: current float salinity
    :param float_theta: current float potential temperature
    :param float_pres: current float pressure
    :return: matrices [number of floats x number of historical profiles] of
    interpolated salinity and interpolated pressure on float theta surface
    """

    # get the shape of the data inputs (float and climatology)
    grid_level = grid_sal.shape[0]
    grid_station = grid_sal.shape[1]
    float_len = float_sal.shape[0]

    # initialise variables to hold the interpolated values
    interp_sal = np.full((float_len, grid_station), np.nan, dtype=np.float64)
    interp_pres = interp_sal

    # check that the climatology data has no infinite (bad) values in the middle
    # of the profiles.
    max_level = 0
    for n in range(0, grid_station):
        # find where data is not missing
        grid_good_sal = np.argwhere(np.isfinite(grid_sal[:, n]))
        grid_good_theta = np.argwhere(np.isfinite(grid_theta[:, n]))
        grid_good_pres = np.argwhere(np.isfinite(grid_pres[:, n]))

        # find indices of good data
        grid_good_data_index = []
        for i in grid_good_sal:
            if i in grid_good_theta and i in grid_good_pres:
                grid_good_data_index.append(i[0])

        # now find the max level
        grid_good_data_len = grid_good_data_index.__len__()
        if grid_good_data_len != 0:
            for m in range(0, grid_good_data_len):
                grid_sal[m, n] = grid_sal[grid_good_data_index[m], n]
                grid_theta[m, n] = grid_theta[grid_good_data_index[m], n]
                grid_pres[m, n] = grid_pres[grid_good_data_index[m], n]
            max_level = np.maximum(max_level, grid_good_data_len)

    # Truncate the number of levels to the maximum level that has
    # available data
    if max_level > 0:
        grid_sal = grid_sal[:max_level, :]
        grid_theta = grid_theta[:max_level, :]
        grid_pres = grid_pres[:max_level, :]
    else:
        raise ValueError("No good climatological data has been found")

    # find where data isn't missing in the float
    float_good_sal = np.isfinite(float_sal)
    float_good_theta = np.isfinite(float_theta)
    float_good_pres = np.isfinite(float_pres)

    # get the indices of the good float data
    float_good_data_index = (float_good_sal == float_good_theta)
    float_good_data_index = np.argwhere(float_good_data_index == float_good_pres)

    # get the number of good float data observations
    float_data_len = float_good_data_index.__len__()

    # compare good float data to the closest climatological data
    for n in range(0, float_data_len):

        # get index of good float data
        index = float_good_data_index[n]

        # Find the difference between the float data and the climatological data
        delta_sal = grid_sal - float_sal[index]
        delta_pres = grid_pres - float_pres[index]
        delta_theta = grid_theta - float_theta[index]

        # Find the indices of the closest pressure value each climatological station has to
        # the float pressures)
        delta_pres_min_index = np.argmin(np.abs(delta_pres), axis=0)

        # go through all the climatological stations
        for m in range(0, grid_station):

            # find if delta_theta is different to delta_theta for closest pressure
            # (equals -1 if different)
            tst = np.sign(delta_theta[:, m]) * np.sign(delta_theta[delta_pres_min_index[n], n])

            # look for a theta match below the float pressure
            grid_theta_below_pres = np.argwhere(tst[delta_pres_min_index[m]:grid_level] < 0)
            min_grid_theta = []

            # look for a theta match above the float pressure
            grid_theta_above_pres = np.argwhere(tst[0:delta_pres_min_index[m]] < 0)
            max_grid_theta = []

            # initialise arrays to hold interpolated pressure and salinity
            interp_pres = []
            interp_sal = []

            # there is a theta value at a deeper level
            if grid_theta_below_pres.__len__() > 0:

                min_grid_theta_index = np.min(grid_theta_below_pres)
                i1 = min_grid_theta_index + delta_pres_min_index[m] - 1
                wt = delta_theta[i1, m] / (delta_theta[i1, m] - delta_theta[i1 - 1, m])
                interp_pres.append(wt * delta_pres[i1 - 1, m] + (1 - wt) * delta_pres[i1, m])
                interp_sal.append(wt * delta_sal[i1 - 1, m] + (1 - wt) * delta_sal[i1, m])

            # there is a theta value at a shallower level
            if grid_theta_above_pres.__len__() > 0:

                i2 = np.max(grid_theta_above_pres) - 1
                wt = delta_theta[i2, m] / (delta_theta[i2, m] - delta_theta[i2 + 1, m])
                interp_pres.append(wt * delta_pres[i2 + 1, m] + (1 - wt) * delta_pres[i2, m])
                interp_sal.append(wt * delta_sal[i2 + 1, m] + (1 - wt) * delta_sal[i2, m])


