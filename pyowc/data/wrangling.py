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
    grid_stations = grid_sal.shape[1]
    float_len = float_sal.shape[0]

    # initialise variables to hold the interpolated values
    interp_sal_final = np.full((float_len, grid_stations), np.nan, dtype=np.float64)
    interp_pres_final = np.full((float_len, grid_stations), np.nan, dtype=np.float64)

    # guard against mismatched arrays
    if not grid_sal.shape == grid_theta.shape == grid_pres.shape:
        return interp_sal_final, interp_pres_final

    grid_data = _get_cleaned_grid_data(grid_stations, grid_sal, grid_theta, grid_pres)

    if grid_data is None:
        return interp_sal_final, interp_pres_final

    grid_sal, grid_theta, grid_pres = grid_data

    # get the indices of the good float data
    float_good_data_index = np.nonzero(_get_finite_element_mask(float_sal, float_theta, float_pres))[0]

    grid_theta_diff = np.diff(grid_theta, axis=0)
    grid_pres_diff = np.diff(grid_pres, axis=0)
    grid_sal_diff = np.diff(grid_sal, axis=0)

    output_placeholder = np.full_like(grid_theta_diff, np.nan, dtype=float)

    # compare good float data to the closest climatological data
    for index in float_good_data_index:
        # calculate the weighting factor for all interpolants
        interp_factor, sign_changes = calculate_interpolation_weights(float_theta[index], grid_theta, grid_theta_diff)

        # interpolate values of pres
        output_placeholder[sign_changes] = interpolate_values(interp_factor, sign_changes, grid_pres, grid_pres_diff)

        # no sign changes between consecutive values for a station means there are no values to interpolate
        stations_with_possible_interp = np.unique(sign_changes[1])

        # find the interpolated pres value which is closest to the float pres for each station
        min_pres = np.nanargmin(np.abs(output_placeholder[:, stations_with_possible_interp] - float_pres[index]), axis=0)

        # use the closest pres value for the interpolated value
        interp_pres_final[index, stations_with_possible_interp] = output_placeholder[min_pres, stations_with_possible_interp]

        # clear any values which have been written to the placeholder
        output_placeholder[sign_changes] = np.nan

        # interpolate values of sal
        output_placeholder[sign_changes] = interpolate_values(interp_factor, sign_changes, grid_sal, grid_sal_diff)

        # select the values of sal matching the closest pres values
        interp_sal_final[index, stations_with_possible_interp] = output_placeholder[min_pres, stations_with_possible_interp]

        # clear any values which have been written to the placeholder
        output_placeholder[sign_changes] = np.nan

    return interp_sal_final, interp_pres_final


def calculate_interpolation_weights(value, reference_values, reference_values_diff):
    """Find interpolation weights of a value on some given data, also returning the locations in the data."""
    diffs = reference_values - value

    # boolean array which is true if a sign change occurred between rows in a column
    sign_changes = find_sign_changes_in_columns(diffs)

    # now convert to indices as we will be using it often
    sign_changes = np.nonzero(sign_changes)

    # find the bracketing values
    value_before = reference_values[:-1][sign_changes]

    # calculate the interpolation weights
    interpolation_factors = (value - value_before) / reference_values_diff[sign_changes]

    return interpolation_factors, sign_changes


def interpolate_values(interpolation_weights, sign_changes, reference, reference_diff):
    """Use iterpolation weights and sign changes previously calculated to interpolate reference data."""
    return reference[:-1][sign_changes] + interpolation_weights * reference_diff[sign_changes]


def find_sign_changes_in_columns(values):
    """Find sign changes which occur in consecutive rows.

    Note: returns True if a value is 0.0 since we're using to find points for interpolation,
        in this case the interpolation coincides with a value and so the weight will be 0.0 or 1.0.
    """
    return values[1:, :] * values[:-1, :] <= 0.0


def _get_cleaned_grid_data(grid_stations, grid_sal_in, grid_theta_in, grid_pres_in):
    """Return a copy of grid data where finite values have been moved up to the top of each station."""
    # ensure no changes made to reference data
    grid_sal = np.full_like(grid_sal_in, np.nan)
    grid_theta = np.full_like(grid_theta_in, np.nan)
    grid_pres = np.full_like(grid_pres_in, np.nan)

    grid_good_data = _get_finite_element_mask(grid_sal_in, grid_theta_in, grid_pres_in)

    # find the max number of levels from all stations
    column_counts = np.count_nonzero(grid_good_data, axis=0)
    max_level = np.max(column_counts)

    for station_id in range(grid_stations):
        good_values = column_counts[station_id]
        grid_sal[:good_values, station_id] = grid_sal_in[:, station_id][grid_good_data[:, station_id]]
        grid_theta[:good_values, station_id] = grid_theta_in[:, station_id][grid_good_data[:, station_id]]
        grid_pres[:good_values, station_id] = grid_pres_in[:, station_id][grid_good_data[:, station_id]]

    # Truncate the number of levels to the maximum level that has available data
    if max_level > 0:
        grid_sal = grid_sal[:max_level, :]
        grid_theta = grid_theta[:max_level, :]
        grid_pres = grid_pres[:max_level, :]
    else:
        print("No good climatological data has been found for this profile")
        return None

    return grid_sal, grid_theta, grid_pres


def _get_finite_element_mask(sal, theta, pres):
    """Return a boolean array where True entries indicate all of sal/theta/pres are finite."""
    # check that the climatology data has no infinite (bad) values in the middle of the profiles.
    good_sal = np.isfinite(sal)
    good_theta = np.isfinite(theta)
    good_pres = np.isfinite(pres)

    # create an array where True indicates that sal/theta/pres are all finite
    return good_sal & good_theta & good_pres
