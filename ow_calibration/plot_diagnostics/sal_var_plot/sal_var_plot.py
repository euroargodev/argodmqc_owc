"""
-----Theta Salinity Plot-----

Written by: Edward Small
When: 10/05/2020

Function for plotting salinity curve aainst potential temperature of all the
data used in the analysis
For information on how to use this file, check the README at either:
https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
from scipy.interpolate import interpolate
from ow_calibration.find_10thetas.find_10thetas import find_10thetas


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
# pylint: disable=too-many-statements
# pylint: disable=too-many-branches
def sal_var_plot(levels, sal, pres, ptmp, map_sal, map_sal_errors,
                 map_ptmp, cal_sal, cal_sal_errors, boundaries, profile_no, float_name):
    """
    Create the salinity variance plot for each level
    :param float_name: name of the float
    :param profile_no: profile numbers
    :param map_ptmp: mapped potential temperature
    :param levels: number of levels to plot
    :param sal: float salinity
    :param pres: float pressure
    :param ptmp: float potential temperature
    :param map_sal: mapped salinity
    :param map_sal_errors: mapped salinity errors
    :param cal_sal: calibrated salinity
    :param cal_sal_errors: calibrated salinity errors
    :param boundaries: pressure and temperature boundaries
    :return: Nothing
    """

    # set up matrices
    profile_no = profile_no.flatten()
    no_profiles = pres.shape[1]
    s_int = np.nan * np.ones((levels, no_profiles))
    s_map = np.nan * np.ones((levels, no_profiles))
    s_map_err = np.nan * np.ones((levels, no_profiles))
    s_cal = np.nan * np.ones((levels, no_profiles))
    s_cal_err = np.nan * np.ones((levels, no_profiles))
    thetalevel_index = np.nan * np.ones((levels, no_profiles))

    # Find levels on which we should plot
    use_theta_lt = boundaries[0]
    use_theta_gt = boundaries[1]
    use_pres_lt = boundaries[2]
    use_pres_gt = boundaries[3]
    use_percent_gt = boundaries[4]

    thetas = find_10thetas(copy.deepcopy(sal),
                           copy.deepcopy(ptmp),
                           copy.deepcopy(pres),
                           copy.deepcopy(map_ptmp),
                           use_theta_lt, use_theta_gt,
                           use_pres_lt, use_pres_gt,
                           use_percent_gt)

    bad = np.argwhere(np.isnan(map_ptmp))
    for i in bad:
        pres[i[0], i[1]] = np.nan
        sal[i[0], i[1]] = np.nan
        ptmp[i[0], i[1]] = np.nan
        map_sal[i[0], i[1]] = np.nan
        map_sal_errors[i[0], i[1]] = np.nan
        cal_sal[i[0], i[1]] = np.nan
        cal_sal_errors[i[0], i[1]] = np.nan

    if use_theta_lt.__len__() > 0 and use_theta_gt.__len__() == 0:
        good = np.argwhere(ptmp > use_theta_lt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_theta_lt.__len__() == 0 and use_theta_gt.__len__() > 0:
        good = np.argwhere(ptmp < use_theta_gt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_theta_lt.__len__() > 0 and use_theta_gt.__len__() > 0:
        theta_range_lt = (ptmp < use_theta_gt)
        theta_range_gt = (ptmp > use_theta_lt)

        if use_theta_lt < use_theta_gt:
            # exclude middle band
            good = np.argwhere(np.logical_and(theta_range_lt, theta_range_gt))

        else:

            good = np.argwhere(np.logical_or(theta_range_gt, theta_range_lt))

        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_pres_lt != 0 and use_pres_gt == 0:
        good = np.argwhere(pres > use_pres_lt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_pres_lt == 0 and use_pres_gt != 0:
        good = np.argwhere(pres < use_pres_gt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_pres_lt != 0 and use_pres_gt != 0:
        pres_range_lt = (pres < use_pres_gt)
        pres_range_gt = (pres > use_pres_lt)

        if use_pres_lt < use_pres_gt:
            # exclude middle band
            good = np.argwhere(np.logical_and(pres_range_lt, pres_range_gt))

        else:

            good = np.argwhere(np.logical_or(pres_range_gt, pres_range_lt))

        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    for i in range(no_profiles):
        for j in range(levels):

            if np.nanmax(ptmp[:, i]) > thetas[0][j] > np.nanmin(ptmp[:, i]):
                diff_theta = np.abs(ptmp[:, i] - thetas[0][j])

                if np.argwhere(~np.isnan(diff_theta)).__len__() == 0:
                    thetalevel_index[j, i] = np.nan

                else:
                    thetalevel_index[j, i] = np.nanmin(np.argwhere(diff_theta ==
                                                                   np.nanmin(diff_theta)))

    # Build s matrix for plotting
    for i in range(levels):
        for j in range(no_profiles):
            theta_index = thetalevel_index[i, j]

            if ~np.isnan(theta_index):
                theta_index = int(theta_index)
                inter = np.arange(np.max([theta_index - 1, 0]),
                                  np.min([theta_index + 1, pres.shape[0] - 1]) + 1,
                                  dtype=int)

                ptmp_diff = ptmp[theta_index, j] - ptmp[inter, j]

                if ptmp[theta_index, j] > thetas[0][i]:
                    pos_diff = np.argwhere(ptmp_diff > 0)

                    if pos_diff.__len__() > 0:
                        min_diff = np.argwhere(ptmp_diff == np.nanmin(ptmp_diff[pos_diff]))
                        k_index = inter[min_diff]

                    else:
                        k_index = theta_index

                if ptmp[theta_index, j] < thetas[0][i]:
                    neg_diff = np.argwhere(ptmp_diff < 0)

                    if neg_diff.__len__() > 0:
                        min_diff = np.argwhere(-ptmp_diff == np.nanmin(-ptmp_diff[neg_diff]))
                        k_index = inter[min_diff]

                    else:
                        k_index = theta_index

                else:
                    k_index = theta_index

                if ((k_index != theta_index and ~np.isnan(sal[theta_index, j])) and
                        ~np.isnan(sal[k_index, j]) and ~np.isnan(ptmp[theta_index, j]) and
                        ~np.isnan(ptmp[k_index, j])):

                    interp_ptmp_sal = interpolate.interp1d([ptmp[theta_index, j],
                                                            ptmp[k_index, j]],
                                                           [sal[theta_index, j],
                                                            sal[k_index, j]])

                    s_int[i, j] = interp_ptmp_sal(thetas[0][i][0])

                else:
                    s_int[i, j] = sal[theta_index, j]

                if ((k_index != theta_index and ~np.isnan(map_sal[theta_index, j])) and
                        ~np.isnan(map_sal[k_index, j]) and ~np.isnan(ptmp[theta_index, j]) and
                        ~np.isnan(ptmp[k_index, j])):

                    interp_map_sal = interpolate.interp1d([ptmp[theta_index, j],
                                                           ptmp[k_index, j]],
                                                          [map_sal[theta_index, j],
                                                           map_sal[k_index, j]])
                    s_map[i, j] = interp_map_sal(thetas[0][i][0])

                    interp_map_sal_err = interpolate.interp1d([ptmp[theta_index, j],
                                                               ptmp[k_index, j]],
                                                              [map_sal_errors[theta_index, j],
                                                               map_sal_errors[k_index, j]])

                    s_map_err[i, j] = interp_map_sal_err(thetas[0][i][0])

                else:
                    s_map[i, j] = map_sal[theta_index, j]
                    s_map_err[i, j] = map_sal_errors[theta_index, j]

                if ((k_index != theta_index and ~np.isnan(cal_sal[theta_index, j])) and
                        ~np.isnan(cal_sal[k_index, j]) and ~np.isnan(ptmp[theta_index, j]) and
                        ~np.isnan(ptmp[k_index, j])):

                    interp_cal_sal = interpolate.interp1d([ptmp[theta_index, j],
                                                           ptmp[k_index, j]],
                                                          [cal_sal[theta_index, j],
                                                           cal_sal[k_index, j]])
                    interp_cal_sal_err = interpolate.interp1d([ptmp[theta_index, j],
                                                               ptmp[k_index, j]],
                                                              [cal_sal_errors[theta_index, j],
                                                               cal_sal_errors[k_index, j]])

                    s_cal[i, j] = interp_cal_sal(thetas[0][i][0])
                    s_cal_err[i, j] = interp_cal_sal_err(thetas[0][i][0])

                else:
                    s_cal[i, j] = cal_sal[theta_index, j]
                    s_cal_err[i, j] = cal_sal_errors[theta_index, j]

    # plot data (one plot for each theta level, as selected by user)
    for i in range(levels):
        plt.figure(1)
        plt.plot(profile_no, s_int[i, :], marker='o', color='b',
                 label='uncalibrated float')
        plt.plot(profile_no, s_map[i, :], color='r',
                 linewidth=4, zorder=0)
        plt.plot(profile_no, s_cal[i, :], color=(0, 1, 0),
                 label='calibrated float w/ 1xerr')
        plt.errorbar(profile_no, s_map[i, :], yerr=s_map_err[i, :], color='r', capsize=2)
        plt.fill_between(profile_no, s_cal[i, :] - s_cal_err[i, :],
                         s_cal[i, :] + s_cal_err[i, :], color=(0, 1, 0))
        plt.plot(profile_no, s_map[i, :], color='r',
                 label='mapped salinity')

        plt.xlabel("Profile number")
        plt.ylabel("PSS-78")

        pl.title(float_name +
                 r" salinities with error on $\theta$=" +
                 str(np.round(thetas[0][i][0], 5)) + r"$\circ$C")

        plt.legend()
        plt.show()
