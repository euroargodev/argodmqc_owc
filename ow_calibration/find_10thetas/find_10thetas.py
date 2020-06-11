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
    min_theta = np.ceil(np.nanmin(PTMP)*10)/10
    max_theta = np.floor(np.nanmax(PTMP)*10)/10

    # only find levels if we have a valid theta range
    if min_theta < max_theta:

        # get pressure levels
        increment = 50
        max_pres = np.nanmax(PRES)
        min_pres = np.nanmin(PRES)
        pres_levels = np.arange(min_pres, max_pres, increment)

        print(pres_levels.__len__())



