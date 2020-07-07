""" Functions dedicated to search

        Parameters
        ----------

        Returns
        -------

"""

import copy
from scipy.interpolate import interpolate
import math
import numpy as np

from pyowc.core import utils
from pyowc.core import stats


#pylint: disable=too-many-arguments
def find_ellipse(data_long, ellipse_long, ellipse_size_long,
                 data_lat, ellipse_lat, ellipse_size_lat,
                 phi, data_pv=0, ellipse_pv=0):
    """ Finds whether a data point exists inside an ellipse

        Calculates whether or not a value exists inside of an ellipse
        of specifed size. If the answer is <1, it exists in the ellipse

        Used to belong in "find_besthist", but was refactored and removed
        to its own file for neatness.

        Parameters
        ----------
        data_long: longitude of the data point
        ellipse_long: longitude of the centre of the ellipse
        ellipse_size_long: size of the ellipse in the longitudinal direction
        data_lat: latitude of the data point
        ellipse_lat: latitude of the centre of the ellipse
        ellipse_size_lat: size of the ellipse in the latitudinal direction
        phi: cross-isobath scale for ellipse
        data_pv: potential vorticity of the data point
        ellipse_pv: potential vorticity of the centre of the ellipse

        Returns
        -------
        float. If <1, it exists inside the ellipse
    """

    total_pv = 0
    if data_pv != 0 and ellipse_pv != 0:
        total_pv = (ellipse_pv - data_pv) / math.sqrt(ellipse_pv ** 2 + data_pv ** 2) / phi

    ellipse = math.sqrt((data_long - ellipse_long) ** 2 / (ellipse_size_long * 3) ** 2 + \
                        (data_lat - ellipse_lat) ** 2 / (ellipse_size_lat * 3) ** 2 + \
                        total_pv ** 2)

    return ellipse


# pylint: disable=fixme
# TODO: ARGODEV-155
# Refactor this code to take objects and dictionaries instead of a ludicrous amount of arguments
# In fact, this function still requires a serious refactor, because it is doing far too much
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
def find_besthist(grid_lat, grid_long, grid_dates, grid_z_value, lat, long, date, z_value,
                  latitude_large, latitude_small, longitude_large, longitude_small,
                  phi_large, phi_small, age_large, age_small, map_pv_use, max_casts):
    """ Finds ln_max_casts number of unique historical data points that are most strongly correlated with the
    float profile being processed


        Find ln_max_casts unique historical points that are most strongly correlated with the float profile

        Rewritten in December 2006 to only use latitude, longitude, data, and water depth as arguments
        and to return index of the station list

        N.B. Change to code on the xx/06/2013: add age_large when computing correlation_large
        - C Cabanes

        N.B Change during conversion to python on the 31/10/2019: Potential Vorticity, Correlation, and
        ellipse are calculated multiple times, so I moved them into their own function.
        These functions can be vectorised using the numpy library (numpy.vectorize(function))
        to use the functions on arrays.
        - Edward Small

        Parameters
        ----------
        grid_lat: array of latitudes of historical data
        grid_long: array of longitudes of historical data
        grid_dates: array of ages of the historical data
        grid_z_value: array of depths of the historical data
        lat: latitude of the float profile
        long: longitude of the float profile
        date: age of the float profile
        z_value: depth of the float profile
        latitude_large: latitude of large ellipse
        latitude_small: latitude of small ellipse
        longitude_large: longitude of large ellipse
        longitude_small: longitude of small ellipse
        phi_large: cross-isobath scale for large ellipse
        phi_small: cross-isobath scale for small ellipse
        age_large: age of data in large ellipse
        age_small: age of data in small ellipse
        map_pv_use: flag for whether to use potential vorticity (see load_configuration.py)
        max_casts: maximum number of data points wanted

        Returns
        -------
        indices of historical data to use
    """

    # make sure arrays are 1 dimensional
    grid_lat = grid_lat.flatten()
    grid_long = grid_long.flatten()
    grid_z_value = grid_z_value.flatten()
    grid_dates = grid_dates.flatten()

    # set up potential vorticity
    potential_vorticity_vec = np.vectorize(utils.potential_vorticity)
    pv_float = 0
    pv_hist = 0

    # if we are using potential vorticity, calculate it
    if map_pv_use == 1:
        pv_float = utils.potential_vorticity(lat, z_value)
        pv_hist = potential_vorticity_vec(grid_lat, grid_z_value)

    # calculate ellipse
    find_ellipse_vec = np.vectorize(find_ellipse)
    ellipse = find_ellipse_vec(grid_long, long, longitude_large, grid_lat, lat, latitude_large,
                               phi_large, pv_hist, pv_float)

    # find points that lie within the ellipse
    hist_long = []
    hist_lat = []
    hist_dates = []
    hist_z_value = []
    index_ellipse = np.argwhere(ellipse < 1)
    for i in index_ellipse:
        hist_long.append(grid_long[i])
        hist_lat.append(grid_lat[i])
        hist_dates.append(grid_dates[i])
        hist_z_value.append(grid_z_value[i])

    index = index_ellipse
    # check to see if too many data points were found
    if index.__len__() > max_casts:

        # uses pseudo-random numbers so the same random numbers are selected for each run
        np.random.seed(index_ellipse.__len__())

        # pick max_casts/3 random points
        rand_matrix = np.random.rand(math.ceil(max_casts / 3))
        index_rand = np.round(rand_matrix * (index_ellipse.__len__() - 1))

        # make sure the points are all unique
        index_rand = np.unique(index_rand)
        # sort them into ascending order
        index_rand.sort()

        # create an array containing the indices of the remaining reference data
        # (index_ellipse array without index_rand array)
        # Then create arrays containing the remaining reference lat, long, age, and z_value
        index_remain = []
        remain_hist_lat = []
        remain_hist_long = []
        remain_hist_z_value = []
        remain_hist_dates = []

        for i in index_ellipse.flatten():
            if i not in index_rand:
                index_remain.append(i)

        for i in index_remain:
            remain_hist_lat = np.append(remain_hist_lat, grid_lat[int(i)])
            remain_hist_long = np.append(remain_hist_long, grid_long[int(i)])
            remain_hist_z_value = np.append(remain_hist_z_value, grid_z_value[int(i)])
            remain_hist_dates = np.append(remain_hist_dates, grid_dates[int(i)])

        # sort remaining points by large spatial correlations
        # if using potential vorticity, calculate it for remaining data
        if map_pv_use == 1:
            pv_hist = potential_vorticity_vec(remain_hist_lat, remain_hist_z_value)

        # calculate the large spatial correlation for each point
        spatial_correlation_vec = np.vectorize(stats.spatial_correlation)
        correlation_large = spatial_correlation_vec(remain_hist_long, long, longitude_large,
                                                    remain_hist_lat, lat, latitude_large,
                                                    remain_hist_dates, date, age_large,
                                                    pv_hist, pv_float, phi_large)

        # combine the large spatial correlation with their indices in a 2d matrix
        correlation_large_combined = np.stack((index_remain, correlation_large)).transpose()
        # sort the matrix in ascending order of correlation
        correlation_large_combined_sorted = correlation_large_combined[
            correlation_large_combined[:, 1].argsort()
        ]

        # work out how many large spatial data points needed
        lsegment2 = 2 * math.ceil(max_casts / 3) - index_rand.__len__()

        # select the best large spatial data points
        index_large_spatial = []
        for i in range(0, lsegment2):
            index_large_spatial = np.append(
                index_large_spatial, correlation_large_combined_sorted[i][0]
            )

        # create array to hold selected indices from random and spatial selection
        rand_large_index = np.concatenate((index_rand, index_large_spatial))
        # sort into ascending order
        rand_large_index.sort()

        # remove the currently selected indices again
        index_remain = []
        remain_hist_lat = []
        remain_hist_long = []
        remain_hist_z_value = []
        remain_hist_dates = []
        for i in index_ellipse.flatten():
            if i not in rand_large_index:
                index_remain.append(i)

        for i in index_remain:
            remain_hist_lat = np.append(remain_hist_lat, grid_lat[int(i)])
            remain_hist_long = np.append(remain_hist_long, grid_long[int(i)])
            remain_hist_z_value = np.append(remain_hist_z_value, grid_z_value[int(i)])
            remain_hist_dates = np.append(remain_hist_dates, grid_dates[int(i)])

        # sort the remaining points by short spatial and temporal correlations
        # if using potential vorticity, calculate it
        if map_pv_use == 1:
            pv_hist = potential_vorticity_vec(
                remain_hist_lat, remain_hist_z_value)

        # calculate the small spatial correlation for each point
        correlation_small = spatial_correlation_vec(remain_hist_long,
                                                    long, longitude_small,
                                                    remain_hist_lat,
                                                    lat, latitude_small,
                                                    remain_hist_dates,
                                                    date, age_small,
                                                    pv_hist, pv_float, phi_small)

        # combine the small spatial correlation with their indices in a 2d matrix
        correlation_small_combined = np.stack((index_remain, correlation_small)).transpose()
        # sort the matrix in ascending order of correlation
        correlation_small_combined_sorted = correlation_small_combined[
            correlation_small_combined[:, 1].argsort()
        ]

        # select the amount we need for short spatial/temporal correlation
        leftover = max_casts - index_rand.__len__() - lsegment2

        index_small_spatial = []
        for i in range(0, leftover):
            index_small_spatial = np.append(
                index_small_spatial, correlation_small_combined_sorted[i][0]
            )

        # put together the best random data, large spatial data, and small spatial data
        index = np.concatenate((index_rand, index_large_spatial, index_small_spatial))
        if index.__len__() != np.unique(index).__len__():
            print("WARNING: not all points are unique")
            rand_large = np.concatenate((index_rand, index_large_spatial))
            rand_small = np.concatenate((index_rand, index_small_spatial))
            large_small = np.concatenate((index_large_spatial, index_small_spatial))
            print("unique points: ", np.unique(index).__len__())
            print("unique rand large: ", rand_large.__len__(),
                  " : ", np.unique(rand_large).__len__())
            print("unique rand small: ", rand_small.__len__(),
                  " : ", np.unique(rand_small).__len__())
            print("unique large small: ", large_small.__len__(),
                  " : ", np.unique(large_small).__len__())

    # ensure that the index is integers
    index = index.flatten()
    index = index.astype(int)

    return index


def nearest_neighbour(x_axis, y_axis, table, x_input, y_input):
    """  Find the nearest neighbour

        Parameters
        ----------
        x_axis: x-axis values
        y_axis: y-axis values
        table: grid data
        x_input: data to interpolate
        y_input: data to interpolate

        Returns
        -------
        nearest neighbour
    """
    x_output = np.abs(x_axis - x_input).argmin()
    y_output = np.abs(y_axis - y_input).argmin()
    return table[y_output, x_output]


def find_25boxes(pn_float_long, pn_float_lat, pa_wmo_boxes):
    """ Find WMO boxes centered on profile

        Finds the 5 x 5 = 25 WMO boxes with the float profile in the centre
        The WMO box numbers, between 90N and 90S are stored in /data/constants/wmo_boxes.mat
        The structure of the matrix is as such:

        Column 1 - box number
        Column 2 - Do we have CTD data (1 = yes, 0 = no)
        Column 3 - Do we have bottle data (1 = yes, 0 = no)
        Column 4 - do we have Argo data (1 = yes, 0 = no)

        N.B. Change to code on the xx/11/2014: extend la_x so interp2 does not think that longitudes in
        the range [5W 5E] are out-of-bound with matlab version >=R2012b - C Cabanes

        N.B. Change during conversion to python on the 01/10/2019. Struggled to find an interpolation
        function that exactly mirrored Matlab's, so I wrote my own - Edward Small

        First, we need to create a look-up table in the form of

            | -5    5   15  ...  355  365
        ----|----------------------------
         85 | 631   1   19  ...  631   1
         75 | 632   2   20  ...  632   2
         65 | 633   3   21  ...  633   3
         ...|...   ...  ... ...  ...  ...
         -85| 648   18  36  ...  648   18

         We do this using 3 matrices:

         - A 1-D matrix for the x axis (la_lookup_x)
         - A 1-D matrix for the y axis (la_lookup_y)
         - A 2-D matrix for the grid data (la_lookup_no)

        Parameters
        ----------
        pn_float_long: float longitude, float
        pn_float_lat: float latitude, float
        pa_wmo_boxes: wmo boxes (explained above), data frame

        Returns
        -------
        (explained above), 25x4 matrix
    """

    la_lookup_x = np.arange(-5, 366, 10, int)

    la_lookup_y = np.arange(85, -86, -10).transpose()

    la_lookup_no = np.full((1, 648), np.arange(1, 649), dtype=int).reshape(36, 18)
    la_lookup_no = np.insert(
        la_lookup_no, 0, la_lookup_no[la_lookup_no.shape[0] - 1]
    ).reshape(37, 18)
    la_lookup_no = np.insert(la_lookup_no, 666, la_lookup_no[1]).reshape(38, 18).transpose()

    # Set up longitudinal and latitudinal values
    ln_x = []
    ln_y = []
    ln_x.append(pn_float_long + .01)
    ln_x.append(pn_float_long + 10.01)
    ln_x.append(pn_float_long - 9.99)
    ln_x.append(pn_float_long + 20.01)
    ln_x.append(pn_float_long - 19.99)

    ln_y.append(pn_float_lat + .01)
    ln_y.append(pn_float_lat + 10.01)
    ln_y.append(pn_float_lat - 9.99)
    ln_y.append(pn_float_lat + 20.01)
    ln_y.append(pn_float_lat - 19.99)

    # wrap longitudinal values
    if ln_x[2] < 0:
        ln_x[2] += 360

    if ln_x[4] < 0:
        ln_x[4] += 360

    if ln_x[0] >= 360:
        ln_x[0] -= 360

    if ln_x[1] >= 360:
        ln_x[1] -= 360

    if ln_x[3] >= 360:
        ln_x[3] -= 360

    if not np.isnan(pn_float_long) and not np.isnan(pn_float_lat):
        ln_i = []
        for i in range(0, 5):
            for j in range(0, 5):
                ln_i.append(
                    nearest_neighbour(la_lookup_x, la_lookup_y, la_lookup_no, ln_x[j], ln_y[i])
                )

    else:
        ln_i = np.full(25, np.nan)

    pa_wmo_numbers = np.full((25, 4), np.nan)
    for i in range(0, 25):
        if not np.isnan(ln_i[i]):
            pa_wmo_numbers[i] = pa_wmo_boxes.get('la_wmo_boxes')[ln_i[i] - 1]

    return pa_wmo_numbers


# pylint:disable=too-many-arguments
# pylint:disable=too-many-locals
# pylint:disable=too-many-branches
# pylint:disable=too-many-statements
def find_10thetas(sal, ptmp, pres, la_ptmp,
                  use_theta_lt=0, use_theta_gt=0,
                  use_pres_lt=0, use_pres_gt=0, use_percent_gt=0.5):
    """ Find on which theta levels salinity variance is lowest

        Chooses 10 theta levels from the float series for use in the linear fit.
        These 10 theta levels are the ones with the minimum S variance on theta.

        These 10 theta levels are distinct (ie. they don't repeat each other).

        Parameters
        ----------
        sal: float salinity
        ptmp: float potential temperature
        pres: float pressure
        la_ptmp: mapped potential temperature
        use_theta_lt: lower bound for potential temperature
        use_theta_gt: upper bound for potential temperature
        use_pres_lt: lower bound for pressure
        use_pres_gt: upper bound for pressure
        use_percent_gt: use percentage greater than

        Returns
        -------
        Theta levels where salinity varies the least
    """

    # We only want 10 theta levels
    no_levels = 10

    # Find how much data we will need to go through
    profile_no = pres.shape[0]
    profile_depth = pres.shape[1]

    # arrays to hold indices with lowest variance and levels
    index = np.empty((no_levels, profile_depth)) * np.nan
    t_levels = np.empty((no_levels, 1)) * np.nan
    p_levels = np.empty((no_levels, 1)) * np.nan
    var_sal_theta = []
    theta_levels = []

    # exclude unmapped, mixed layer
    unmapped = np.argwhere(np.isnan(la_ptmp))

    for i in unmapped:
        pres[i[0], i[1]] = np.nan
        sal[i[0], i[1]] = np.nan
        ptmp[i[0], i[1]] = np.nan

    # only use theta and pressure from specified range

    if use_theta_lt != 0 and use_theta_gt == 0:
        theta_range = np.argwhere(ptmp > use_theta_lt)
        for i in theta_range:
            pres[i[0], i[1]] = np.nan
            sal[i[0], i[1]] = np.nan
            ptmp[i[0], i[1]] = np.nan

    if use_theta_lt == 0 and use_theta_gt != 0:
        theta_range = np.argwhere(ptmp < use_theta_gt)
        for i in theta_range:
            pres[i[0], i[1]] = np.nan
            sal[i[0], i[1]] = np.nan
            ptmp[i[0], i[1]] = np.nan

    if use_theta_lt != 0 and use_theta_gt != 0:

        if use_theta_lt < use_theta_gt:
            # exclude middle band
            theta_range = np.argwhere(np.logical_and(use_theta_lt < ptmp, ptmp < use_theta_gt))

        else:
            theta_range = np.argwhere(np.logical_xor(ptmp < use_theta_gt, ptmp > use_theta_lt))

        for i in theta_range:

            pres[i[0], i[1]] = np.nan
            sal[i[0], i[1]] = np.nan
            ptmp[i[0], i[1]] = np.nan

    if use_pres_lt != 0 and use_pres_gt == 0:
        pres_range = np.argwhere(pres > use_pres_lt)
        for i in pres_range:
            pres[i[0], i[1]] = np.nan
            sal[i[0], i[1]] = np.nan
            ptmp[i[0], i[1]] = np.nan

    if use_pres_lt == 0 and use_pres_gt != 0:
        pres_range = np.argwhere(pres < use_pres_gt)
        for i in pres_range:
            pres[i[0], i[1]] = np.nan
            sal[i[0], i[1]] = np.nan
            ptmp[i[0], i[1]] = np.nan

    if use_pres_lt != 0 and use_pres_gt != 0:
        if use_pres_lt < use_pres_gt:
            pres_range = np.argwhere(np.logical_and(use_pres_lt < pres, pres < use_pres_gt))

        else:
            pres_range = np.argwhere(np.logical_xor(pres < use_pres_gt, pres > use_pres_lt))

        for i in pres_range:
            pres[i[0], i[1]] = np.nan
            sal[i[0], i[1]] = np.nan
            ptmp[i[0], i[1]] = np.nan

    # find minimum and maximum theta
    min_theta = np.ceil(np.nanmin(ptmp) * 10) / 10
    max_theta = np.floor(np.nanmax(ptmp) * 10) / 10

    # only find levels if we have a valid theta range
    if min_theta < max_theta:

        # get pressure levels
        increment = 50
        max_pres = np.nanmax(pres)
        min_pres = np.nanmin(pres)
        pres_levels = np.arange(min_pres, max_pres, increment)

        # check we can get 10 theta levels. If not, alter pressure increment
        if pres_levels.__len__() < no_levels:
            increment = np.floor((max_pres - min_pres) / no_levels)
            pres_levels = np.arange(min_pres, max_pres, increment)

    # interpolate levels onto pressure increments

    interp_t = np.empty((pres_levels.__len__(), profile_depth)) * np.nan

    for depth in range(profile_depth):
        good = np.argwhere(np.logical_and(~np.isnan(pres[:, depth]),
                                          ~np.isnan(ptmp[:, depth])))

        for pres_i in range(pres_levels.__len__()):
            if np.max(pres[good, depth]) > pres_levels[pres_i] > np.min(pres[good, depth]):
                interp = interpolate.interp1d(pres[good, depth].flatten(),
                                              ptmp[good, depth].flatten())
                interp_t[pres_i, depth] = interp(pres_levels[pres_i])

    # find mean of the interpolated pressure at each level

    theta_levels = np.empty((pres_levels.__len__(), 1)) * np.nan
    theta_level_indices = np.empty((pres_levels.__len__(), profile_depth)) * np.nan

    for pres_i in range(pres_levels.__len__()):
        good_interp_t = np.argwhere(~np.isnan(interp_t[pres_i, :]))

        if good_interp_t.__len__() > 0:
            theta_levels[pres_i] = np.nanmean(interp_t[pres_i, good_interp_t])

    # find profile levels closest to theta levels
    # Areas with temperature inversions (eg Southern Ocean, Gulf of Alaska) PTMP is not unique
    # so this will pick out indices from different depths for the same temperature,
    # thus giving an artificially high salinity variance. This is okay because temperature
    # inversions usually have naturally high variance.

    # find indices that minimise theta on each level

    for depth in range(profile_depth):

        for level in range(theta_levels.__len__()):
            theta_diff = np.array([np.nan])

            if np.nanmax(ptmp[:, depth]) > theta_levels[level] > np.nanmin(ptmp[:, depth]):
                theta_diff = np.abs(ptmp[:, depth] - theta_levels[level])

            if np.all(np.isnan(theta_diff)):
                theta_level_indices[level, depth] = np.nan

            else:
                theta_level_indices[level, depth] = np.min(np.argwhere(theta_diff ==
                                                                       np.nanmin(theta_diff)))

    # find salinity variance on these theta levels

    sal_temp = np.empty((theta_levels.__len__(), profile_depth)) * np.nan

    for level in range(theta_levels.__len__()):
        for depth in range(profile_depth):
            theta_index = theta_level_indices[level, depth]

            # only continue if we have a good index
            if ~np.isnan(theta_index):
                theta_index = int(theta_index)
                interval = np.arange(np.max([theta_index - 1, 0]),
                                     np.min([theta_index + 1, profile_no - 1]) + 1,
                                     dtype=int)

                ptmp_diff = ptmp[theta_index, depth] - ptmp[interval, depth]

                if ptmp[theta_index, depth] > theta_levels[level]:
                    pos_diff = np.argwhere(ptmp_diff > 0)

                    if pos_diff.__len__() > 0:
                        min_diff = np.argwhere(ptmp_diff == np.nanmin(ptmp_diff[pos_diff]))
                        k_index = interval[min_diff]

                    else:
                        k_index = theta_index

                if ptmp[theta_index, depth] < theta_levels[level]:
                    neg_diff = np.argwhere(ptmp_diff < 0)

                    if neg_diff.__len__() > 0:
                        min_diff = np.argwhere(-ptmp_diff == np.nanmin(-ptmp_diff[neg_diff]))
                        k_index = interval[min_diff]

                    else:
                        k_index = theta_index

                # else we only have one profile
                if ptmp[theta_index, depth] == theta_levels[level]:
                    k_index = theta_index

                # interpolate theta level, if possible
                if (k_index != theta_index and ~np.isnan(sal[theta_index, depth]) and
                        ~np.isnan(sal[k_index, depth]) and ~np.isnan(ptmp[theta_index, depth]) and
                        ~np.isnan(ptmp[k_index, depth])):
                    interp_ptmp_sal = interpolate.interp1d([ptmp[theta_index, depth],
                                                            ptmp[k_index, depth]],
                                                           [sal[theta_index, depth],
                                                            sal[k_index, depth]])

                    sal_temp[level, depth] = interp_ptmp_sal(theta_levels[level])

                # else we use the closest points
                else:
                    sal_temp[level, depth] = sal[theta_index, depth]

    num_good = np.empty((theta_levels.__len__(), 1)) * np.nan
    percent_s_profs = np.empty((theta_levels.__len__(), 1)) * np.nan
    var_sal_tlevels = np.empty((theta_levels.__len__(), 1)) * np.nan

    # only use salinities on theta levels that have valid values

    for i in range(theta_levels.__len__()):
        good = np.argwhere(~np.isnan(sal_temp[i, :]))
        num_good[i] = good.__len__()

        if num_good.__len__() > 0:
            var_sal_tlevels[i] = np.nanvar(sal_temp[i, good], ddof=1)

    for j in range(theta_levels.__len__()):
        if np.nanmax(num_good) != 0:
            percent_s_profs[j] = num_good[j] / np.nanmax(num_good)

    bad = np.argwhere(percent_s_profs < use_percent_gt)
    var_sal_tlevels[bad] = np.nan
    var_sal_theta = copy.deepcopy(var_sal_tlevels)

    # select the best 10 theta levels

    for i in range(no_levels):
        min_theta_index = np.argwhere(var_sal_tlevels == np.nanmin(var_sal_tlevels))[0, 0]
        index[i, :] = theta_level_indices[min_theta_index, :]
        t_levels[i] = theta_levels[min_theta_index]
        p_levels[i] = pres_levels[min_theta_index]
        var_sal_tlevels[min_theta_index] = np.nan

    return t_levels, p_levels, index, var_sal_theta, theta_levels
