from scipy.interpolate import interpolate
import numpy as np


def find_10thetas(SAL, PTMP, PRES, la_ptmp,
                  use_theta_lt=0, use_theta_gt=0,
                  use_pres_lt=0, use_pres_gt=0, use_percent_gt=0.5):
    """
    Find on which theta levels salinity variance is lowest
    :param SAL: float salinity
    :param PTMP: float potential temperature
    :param PRES: float pressure
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
    profile_no = PRES.shape[0]
    profile_depth = PRES.shape[1]

    # arrays to hold indices with lowest variance and levels
    index = np.empty((no_levels, profile_depth)) * np.nan
    t_levels = np.empty((no_levels, 1)) * np.nan
    p_levels = np.empty((no_levels, 1)) * np.nan
    var_sal_theta = []
    theta_levels = []

    # exclude unmapped, mixed layer
    unmapped = np.argwhere(np.isnan(la_ptmp))

    for i in unmapped:
        PRES[i[0], i[1]] = np.nan
        SAL[i[0], i[1]] = np.nan
        PTMP[i[0], i[1]] = np.nan

    # only use theta and pressure from specified range

    if use_theta_lt != 0 and use_theta_gt == 0:
        theta_range = np.argwhere(PTMP > use_theta_lt)
        for i in theta_range:
            PRES[i[0], i[1]] = np.nan
            SAL[i[0], i[1]] = np.nan
            PTMP[i[0], i[1]] = np.nan

    if use_theta_lt == 0 and use_theta_gt != 0:
        theta_range = np.argwhere(PTMP < use_theta_gt)
        for i in theta_range:
            PRES[i[0], i[1]] = np.nan
            SAL[i[0], i[1]] = np.nan
            PTMP[i[0], i[1]] = np.nan

    if use_theta_lt != 0 and use_theta_gt != 0:

        if use_theta_lt < use_theta_gt:
            # exclude middle band
            theta_range = np.argwhere(use_theta_lt < PTMP < use_theta_gt)

        else:
            theta_range = np.argwhere(PTMP < use_theta_gt or PTMP > use_theta_lt)

        for i in theta_range:
            PRES[i[0], i[1]] = np.nan
            SAL[i[0], i[1]] = np.nan
            PTMP[i[0], i[1]] = np.nan

    if use_pres_lt != 0 and use_pres_gt == 0:
        pres_range = np.argwhere(PRES > use_pres_lt)
        for i in pres_range:
            PRES[i[0], i[1]] = np.nan
            SAL[i[0], i[1]] = np.nan
            PTMP[i[0], i[1]] = np.nan

    if use_pres_lt == 0 and use_pres_gt != 0:
        pres_range = np.argwhere(PTMP < use_pres_gt)
        for i in pres_range:
            PRES[i[0], i[1]] = np.nan
            SAL[i[0], i[1]] = np.nan
            PTMP[i[0], i[1]] = np.nan

    if use_pres_lt != 0 and use_pres_gt != 0:

        if use_pres_lt < use_pres_gt:
            # exclude middle band
            pres_range = np.argwhere(use_pres_lt < PTMP < use_pres_gt)

        else:
            pres_range = np.argwhere(PTMP < use_pres_gt or PTMP > use_pres_lt)

        for i in pres_range:
            PRES[i[0], i[1]] = np.nan
            SAL[i[0], i[1]] = np.nan
            PTMP[i[0], i[1]] = np.nan

    # find minimum and maximum theta
    min_theta = np.ceil(np.nanmin(PTMP) * 10) / 10
    max_theta = np.floor(np.nanmax(PTMP) * 10) / 10

    # only find levels if we have a valid theta range
    if min_theta < max_theta:

        # get pressure levels
        increment = 50
        max_pres = np.nanmax(PRES)
        min_pres = np.nanmin(PRES)
        pres_levels = np.arange(min_pres, max_pres, increment)

        # check we can get 10 theta levels. If not, alter pressure increment
        if pres_levels.__len__() < no_levels:
            increment = np.floor((max_pres - min_pres) / no_levels)
            pres_levels = np.arange(min_pres, max_pres, increment)

    # interpolate levels onto pressure increments

    interp_t = np.empty((pres_levels.__len__(), profile_depth)) * np.nan

    for depth in range(profile_depth):
        good = np.argwhere(~np.isnan(PRES[:, depth]) & ~np.isnan(PTMP[:, depth]))

        for pres in range(pres_levels.__len__()):

            if np.max(PRES[good, depth]) > pres_levels[pres] > np.min(PRES[good, depth]):
                interp = interpolate.interp1d(PRES[good, depth].flatten(),
                                              PTMP[good, depth].flatten())
                interp_t[pres, depth] = interp(pres_levels[pres])

    # find mean of the interpolated pressure at each level

    theta_levels = np.empty((pres_levels.__len__(), 1)) * np.nan
    theta_level_indices = np.empty((pres_levels.__len__(), profile_depth)) * np.nan

    for pres in range(pres_levels.__len__()):
        good_interp_t = np.argwhere(~np.isnan(interp_t[pres, :]))

        if good_interp_t.__len__() > 0:
            theta_levels[pres] = np.nanmean(interp_t[pres, good_interp_t])

    # find profile levels closest to theta levels
    # Areas with temperature inversions (eg Southern Ocean, Gulf of Alaska) PTMP is not unique
    # so this will pick out indices from different depths for the same temperature,
    # thus giving an artificially high salinity variance. This is okay because temperature
    # inversions usually have naturally high variance.

    # find indices that minimise theta on each level

    for depth in range(profile_depth):

        for level in range(theta_levels.__len__()):
            theta_diff = np.array([np.nan])

            if np.nanmax(PTMP[:, depth]) > theta_levels[level] > np.nanmin(PTMP[:, depth]):
                theta_diff = np.abs(PTMP[:, depth] - theta_levels[level])

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

                ptmp_diff = PTMP[theta_index, depth] - PTMP[interval, depth]

                if PTMP[theta_index, depth] > theta_levels[level]:
                    pos_diff = np.argwhere(ptmp_diff > 0)

                    if pos_diff.__len__() > 0:
                        min_diff = np.argwhere(ptmp_diff == np.nanmin(ptmp_diff[pos_diff]))
                        k_index = interval[min_diff]

                    else:
                        k_index = theta_index

                if PTMP[theta_index, depth] < theta_levels[level]:
                    neg_diff = np.argwhere(ptmp_diff < 0)

                    if neg_diff.__len__() > 0:
                        min_diff = np.argwhere(-ptmp_diff == np.nanmin(-ptmp_diff[neg_diff]))
                        k_index = interval[min_diff]

                    else:
                        k_index = theta_index

                # else we only have one profile
                if PTMP[theta_index, depth] == theta_levels[level]:
                    k_index = theta_index

                # interpolate theta level, if possible
                if (k_index != theta_index and ~np.isnan(SAL[theta_index, depth]) and
                        ~np.isnan(SAL[k_index, depth]) and ~np.isnan(PTMP[theta_index, depth]) and
                        ~np.isnan(PTMP[k_index, depth])):
                    interp_ptmp_sal = interpolate.interp1d([PTMP[theta_index, depth],
                                                            PTMP[k_index, depth]],
                                                           [SAL[theta_index, depth],
                                                            SAL[k_index, depth]])

                    sal_temp[level, depth] = interp_ptmp_sal(theta_levels[level])

                # else we use the closest points
                else:
                    sal_temp[level, depth] = SAL[theta_index, depth]

    num_good = np.empty((theta_levels.__len__(), 1)) * np.nan
    percent_s_profs = np.empty((theta_levels.__len__(), 1)) * np.nan
    var_sal_tlevels = np.empty((theta_levels.__len__(), 1)) * np.nan

    # only use salinities on theta levels that have valid values

    for i in range(theta_levels.__len__()):
        good = np.argwhere(~np.isnan(sal_temp[i, :]))
        num_good[i] = good.__len__()

        if num_good.__len__() > 0:
            var_sal_tlevels[i] = np.var(sal_temp[i, good], ddof=1)

