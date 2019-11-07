"""
-----find best historical casts-----

Written by: Breck Owens and Annie Wong
When: xx/12/2006
Converted to python by: Edward Small
When: 30/10/2019

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

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import math
import numpy as np
from ..potential_vorticity.potential_vorticity import potential_vorticity
from ..spatial_correlation.spatial_correlation import spatial_correlation
from ..find_ellipse.find_ellipse import find_ellipse

# pylint: disable=fixme
# TODO: ARGODEV-155
# Refactor this code to take objects and dictionaries instead of a ludicrous amount of arguments
# In fact, this function still requires a serious refactor, because it is doing far too much
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
def find_besthist(
        grid_lat, grid_long, grid_dates, grid_z_value,
        lat, long, date, z_value,
        latitude_large, latitude_small, longitude_large, longitude_small,
        phi_large, phi_small, age_large, age_small, map_pv_use, max_casts
):
    """
    Finds ln_max_casts number of unique historical data points that are most strongly correlated
    with the float profile being processed
    :param grid_lat: array of latitudes of historical data
    :param grid_long: array of longitudes of historical data
    :param grid_dates: array of ages of the historical data
    :param grid_z_value: array of depths of the historical data
    :param lat: latitude of the float profile
    :param long: longitude of the float profile
    :param date: age of the float profile
    :param z_value: depth of the float profile
    :param latitude_large: latitude of large ellipse
    :param latitude_small: latitude of small ellipse
    :param longitude_large: longitude of large ellipse
    :param longitude_small: longitude of small ellipse
    :param phi_large: cross-isobath scale for large ellipse
    :param phi_small: cross-isobath scale for small ellipse
    :param age_large: age of data in large ellipse
    :param age_small: age of data in small ellipse
    :param map_pv_use: flag for whether to use potential vorticity (see load_configuration.py)
    :param max_casts: maximum number of data points wanted
    :return: indices of historical data to use
    """

    # make sure arrays are 1 dimensional
    grid_lat = grid_lat.flatten()
    grid_long = grid_long.flatten()
    grid_z_value = grid_z_value.flatten()
    grid_dates = grid_dates.flatten()

    # set up potential vorticity
    potential_vorticity_vec = np.vectorize(potential_vorticity)
    pv_float = 0
    pv_hist = 0

    # if we are using potential vorticity, calculate it
    if map_pv_use == 1:
        pv_float = potential_vorticity(lat, z_value)
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

        # calculate the large spatial corellation for each point
        spatial_correlation_vec = np.vectorize(spatial_correlation)
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

    index = np.unique(index)
    return index
