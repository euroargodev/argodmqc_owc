""" Functions to manipulate and process data

        Parameters
        ----------

        Returns
        -------

"""

import numpy as np
from ..core.stats import covar_xyt_pv


#pylint: disable=too-many-arguments
#pylint: disable=too-many-locals
def map_data_grid(sal, grid_pos, data_pos, lat, long, age,
                  signal_variance, noise_variance, phi, map_pv_use):
    """ Maps historical float data onto a single float

        An optimal mapping routine, taking data measured in arbitrary geographic locations and
        mapping these data onto a more regular grid. As should happen in every mapping problem
        the data are both mapped onto the prescribed grid, and the data locations (so that
        the mapped field can be checked with the original data to ensure that the statistics
        are valid and consistent).

        Before objectively mapping the data, a mean, using the correlation scales to define the
        weights for the sum (see Bretherton, etal, 1975) is removed.  The error estimate includes
        the contributions from both the mapping and the mean estimation.

        Parameters
        ----------
        sal: array of salinities of the historical float data
        grid_pos: array containing single float data [lat, long, age, depth]
        data_pos: n*4 array containing historical float data [lat, long, age, depth]
        lat: scalar latitude
        long: scalar longitude
        age: scalar age
        signal_variance: scalar signal variance
        noise_variance: scalar noise variance
        phi: scalar cross isobaric scale
        map_pv_use: flag for including vorticity (1=include)

        Returns
        -------
        Tuple containing mapped fields, error estimates of mapped fields, mapped fields on original locations, and their error estimates
    """

    # create the data-data covariance matrix
    data_pos_covar = covar_xyt_pv(data_pos, data_pos, lat, long, age, phi, map_pv_use)
    data_data_covar = np.linalg.inv(signal_variance * data_pos_covar +
                                    noise_variance * np.identity(data_pos.__len__()))

    # estimate the mean field and weights
    sum_data_data_covar = sum(sum(data_data_covar))
    mean_field = sum(np.dot(data_data_covar, sal.transpose())) / sum_data_data_covar
    weight = np.dot(data_data_covar, (sal - mean_field))

    # calculate the objectively mapped fields on position data and grid data
    pos_data_covar = signal_variance * data_pos_covar
    data_weight_covar = np.dot(pos_data_covar, weight) + mean_field

    # include the error in the mean (from Brethrerton, 1975)
    dot_covar_diag = np.diag(np.dot(
        np.dot(pos_data_covar, data_data_covar),
        np.transpose(pos_data_covar)))
    covar_sum = np.sum(np.dot(pos_data_covar, data_data_covar), axis=1)
    data_weight_covar_error = np.sqrt(signal_variance - dot_covar_diag +
                                      ((1 - covar_sum) ** 2) / sum_data_data_covar)

    # now map to the data to the regular grid
    grid_data_covar = (signal_variance * covar_xyt_pv(data_pos, grid_pos, lat, long,
                                                      age, phi, map_pv_use)).transpose()
    grid_weight_covar = np.dot(grid_data_covar, weight) + mean_field
    dot_covar_diag = np.diag(np.dot(
        np.dot(grid_data_covar, data_data_covar), np.transpose(grid_data_covar)))
    covar_sum = np.sum(np.dot(grid_data_covar, data_data_covar), axis=1)
    grid_weight_covar_error = np.sqrt(signal_variance - dot_covar_diag +
                                      ((1 - covar_sum) ** 2) / sum_data_data_covar)

    return grid_weight_covar[0], grid_weight_covar_error[0], \
           data_weight_covar, data_weight_covar_error


#pylint: disable=too-many-arguments
#pylint: disable=too-many-locals
#pylint: disable=too-many-branches
#pylint: disable=too-many-statements
def interp_climatology(grid_sal, grid_theta, grid_pres, float_sal, float_theta, float_pres):
    """ Interpolate historical salinity and pressure data on the float theta

        Routine to interpolate climatological salinity and pressure data onto float's potential temperature.

        Parameters
        ----------
        grid_sal: historical salinity
        grid_theta: historical potential temperature
        grid_pres: historical pressure
        float_sal: current float salinity
        float_theta: current float potential temperature
        float_pres: current float pressure

        Returns
        -------
        Matrices [number of floats x number of historical profiles] of interpolated salinity and interpolated pressure on float theta surface
    """
    # get the shape of the data inputs (float and climatology)
    grid_level = grid_sal.shape[0]
    grid_stations = grid_sal.shape[1]
    float_len = float_sal.shape[0]

    # initialise variables to hold the interpolated values
    interp_sal_final = np.full((float_len, grid_stations), np.nan, dtype=np.float64)
    interp_pres_final = np.full((float_len, grid_stations), np.nan, dtype=np.float64)

    # check that the climatology data has no infinite (bad) values in the middle
    # of the profiles.
    max_level = 0
    for i in range(0, grid_stations):

        # find where data is not missing
        grid_good_sal = np.isfinite(grid_sal[:, i])
        grid_good_theta = np.isfinite(grid_theta[:, i])
        grid_good_pres = np.isfinite(grid_pres[:, i])

        # find indices of good data
        grid_good_data_index = []
        for j in range(0, grid_level):
            if grid_good_sal[j] and grid_good_theta[j] and grid_good_pres[j]:
                grid_good_data_index.append(j)

        # now find the max level
        grid_good_data_len = grid_good_data_index.__len__()
        for j in range(0, grid_good_data_len):
            grid_sal[j, i] = grid_sal[grid_good_data_index[j], i]
            grid_theta[j, i] = grid_theta[grid_good_data_index[j], i]
            grid_pres[j, i] = grid_pres[grid_good_data_index[j], i]

        max_level = np.maximum(max_level, grid_good_data_len)

    # Truncate the number of levels to the maximum level that has
    # available data
    if max_level > 0:
        grid_sal = grid_sal[:max_level, :]
        grid_theta = grid_theta[:max_level, :]
        grid_pres = grid_pres[:max_level, :]
    else:
        print("No good climatological data has been found for this profile")
        return interp_sal_final, interp_pres_final

    # find where data isn't missing in the float
    float_good_sal = np.isfinite(float_sal)
    float_good_theta = np.isfinite(float_theta)
    float_good_pres = np.isfinite(float_pres)

    # get the indices of the good float data
    float_good_data_index = []
    for i in range(0, float_len):
        if float_good_sal[i] and float_good_theta[i] and float_good_pres[i]:
            float_good_data_index.append(i)

    # get the number of good float data observations
    float_data_len = float_good_data_index.__len__()

    # compare good float data to the closest climatological data
    for i in range(0, float_data_len):

        # get index of good float data
        index = float_good_data_index[i]

        # Find the difference between the float data and the climatological data
        delta_sal = np.array(grid_sal - float_sal[index])
        delta_pres = np.array(grid_pres - float_pres[index])
        delta_theta = np.array(grid_theta - float_theta[index])

        # Find the indices of the closest pressure value each climatological station has to
        # the float pressures
        abs_delta_pres = np.abs(delta_pres)
        delta_pres_min_index = np.nanargmin(abs_delta_pres, axis=0)

        # go through all the climatological stations
        for j in range(0, grid_stations):

            # find if delta_theta is different to delta_theta for closest pressure
            # (equals -1 if different)
            tst = np.sign(delta_theta[:, j]) * np.sign(delta_theta[delta_pres_min_index[j], j])

            # look for a theta match below the float pressure
            grid_theta_below_pres = np.argwhere(tst[delta_pres_min_index[j]:grid_level] < 0)

            # look for a theta match above the float pressure
            grid_theta_above_pres = np.argwhere(tst[0:delta_pres_min_index[j]] < 0)
            # initialise arrays to hold interpolated pressure and salinity
            interp_pres = []
            interp_sal = []

            # there is a theta value at a deeper level
            if grid_theta_below_pres.__len__() > 0:
                min_grid_theta_index = np.min(grid_theta_below_pres)
                i_1 = min_grid_theta_index + delta_pres_min_index[j]
                w_t = delta_theta[i_1, j] / (delta_theta[i_1, j] - delta_theta[i_1 - 1, j])
                interp_pres.append(w_t * delta_pres[i_1 - 1, j] + (1 - w_t) * delta_pres[i_1, j])
                interp_sal.append(w_t * delta_sal[i_1 - 1, j] + (1 - w_t) * delta_sal[i_1, j])

            # there is a theta value at a shallower level
            if grid_theta_above_pres.__len__() > 0:
                i_2 = np.max(grid_theta_above_pres)
                w_t = delta_theta[i_2, j] / (delta_theta[i_2, j] - delta_theta[i_2 + 1, j])
                interp_pres.append(w_t * delta_pres[i_2 + 1, j] + (1 - w_t) * delta_pres[i_2, j])
                interp_sal.append(w_t * delta_sal[i_2 + 1, j] + (1 - w_t) * delta_sal[i_2, j])

            if interp_pres.__len__() > 0:
                # if there are two nearby theta values, choose the closest one
                abs_interp_pres = np.abs(interp_pres)
                k = np.argmin(abs_interp_pres)
                interp_sal_final[index, j] = interp_sal[k] + float_sal[index]
                interp_pres_final[index, j] = interp_pres[k] + float_pres[index]

    return interp_sal_final, interp_pres_final
