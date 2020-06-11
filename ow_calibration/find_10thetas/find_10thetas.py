import numpy as np

def find_10thetas(SAL, PTMP, PRES, la_ptmp,
                  use_theta_lt, use_theta_gt,
                  use_pres_lt, use_pres_gt, use_percent_gt):
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



