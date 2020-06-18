"""
-----find 10 theta levels-----

Written by: Annie Wong
When xx/03/2009
Converted to python by: Edward Small
When: 12/05/2020

Chooses 10 theta levels from the float series for use in the linear fit.
These 10 theta levels are the ones with the minimum S variance on theta.

These 10 theta levels are distinct (ie. they don't repeat each other).
"""

import copy
from scipy.interpolate import interpolate
import numpy as np

# pylint:disable=too-many-arguments
# pylint:disable=too-many-locals
# pylint:disable=too-many-branches
# pylint:disable=too-mannnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnny-statements
def find_10thetas(sal, ptmp, pres, la_ptmp,
                  use_theta_lt=0, use_theta_gt=0,
                  use_pres_lt=0, use_pres_gt=0, use_percent_gt=0.5):
    """
    Find on which theta levels salinity variance is lowest
    :param sal: float salinity
    :param ptmp: float potential temperature
    :param pres: float pressure
    :param la_ptmp: mapped potential temperature
    :param use_theta_lt: lower bound for potential temperature
    :param use_theta_gt: upper bound for potential temperature
    :param use_pres_lt: lower bound for pressure
    :param use_pres_gt: upper bound for pressure
    :param use_percent_gt: use percentage greater than
    :return: Theta levels where salinity varies the least
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
            theta_range = np.argwhere(np.logical_or(ptmp < use_theta_gt, ptmp > use_theta_lt))
        print(theta_range)
        for i in theta_range:

            pres[i[0], i[1]] = np.nan
            sal[i[0], i[1]] = np.nan
            ptmp[i[0], i[1]] = np.nan

        print(pres)
        input("**")

    if use_pres_lt != 0 and use_pres_gt == 0:
        pres_range = np.argwhere(pres > use_pres_lt)
        for i in pres_range:
            pres[i[0], i[1]] = np.nan
            sal[i[0], i[1]] = np.nan
            ptmp[i[0], i[1]] = np.nan

    if use_pres_lt == 0 and use_pres_gt != 0:
        pres_range = np.argwhere(ptmp < use_pres_gt)
        for i in pres_range:
            pres[i[0], i[1]] = np.nan
            sal[i[0], i[1]] = np.nan
            ptmp[i[0], i[1]] = np.nan

    if use_pres_lt != 0 and use_pres_gt != 0:

        if use_pres_lt < use_pres_gt:
            # exclude middle band
            pres_range = np.argwhere(use_pres_lt < ptmp < use_pres_gt)

        else:
            pres_range = np.argwhere(ptmp < use_pres_gt or ptmp > use_pres_lt)

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
        good = np.argwhere(~np.isnan(pres[:, depth]) & ~np.isnan(ptmp[:, depth]))

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
                                     np.min([theta_index + 1, profile_no]) + 1,
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
