"""
-----Theta Salinity Plot-----

Written by: Edward Small
When: 10/05/2020

Function for plotting the calibrated salinity curve

For information on how to use this file, check the README at either:
https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
def cal_sal_curve_plot(sal, cal_sal, cal_sal_err, sta_sal, sta_sal_err, sta_mean,
                       pcond_factor, pcond_factor_err, profile_no, float_name):
    """
    Create the calibrated salinity curve plot
    :param sta_mean: mean of the differences
    :param profile_no: profile numbers
    :param sta_sal_err: error in mean difference
    :param cal_sal_err: error in calibrated salinity
    :param sal: float salinity
    :param cal_sal: calibrated salinity
    :param sta_sal: mean difference between salinity and calculated salinity per profile
    :param pcond_factor: slope
    :param pcond_factor_err: slope error
    :param float_name: name of the float
    :return: Nothing
    """

    # only plot this if we have calibrated salinity
    if np.argwhere(~np.isnan(cal_sal)).__len__() > 0:
        # set up matrices
        profiles = sal.shape[1]
        sal_offset = cal_sal - sal
        avg_sal_offset = np.nan * np.ones((1, profiles)).flatten()
        avg_sal_offset_err = np.nan * np.ones((1, profiles)).flatten()
        sta_offset = sta_sal - sal
        avg_sta_offset = np.nan * np.ones((1, profiles)).flatten()
        avg_sta_offset_err = np.nan * np.ones((1, profiles)).flatten()

        for i in range(profiles):
            avg_sal_offset[i] = np.nanmean(sal_offset[:, i])
            avg_sal_offset_err[i] = np.nanmean(cal_sal_err[:, i])
            avg_sta_offset[i] = np.nanmean(sta_offset[:, i])
            avg_sta_offset_err[i] = np.nanmean(sta_sal_err[:, i])

        # begin plotting data (two plots)

        # potential conductivity plot
        plt.figure(1)
        plt.subplot(211)

        plt.plot(profile_no[0], pcond_factor[0], color=(0, 1, 0), linestyle='-', linewidth=1)
        good = np.argwhere(np.isfinite(sta_mean))
        plt.plot(profile_no[0, good[:, 1]], sta_mean[0, good[:, 1]],
                 color=(1, 0, 0), linestyle='-', linewidth=1,
                 label='1-1 profile fit')

        plt.errorbar(profile_no[0], pcond_factor[0], yerr=2 * pcond_factor_err[0],
                     color=(0, 0, 1), linestyle='', capsize=2, label='2 x cal error')
        plt.errorbar(profile_no[0], pcond_factor[0], yerr=pcond_factor_err[0],
                     color=(0, 1, 0), linestyle='', capsize=2, label='1 x cal error')

        plt.legend()
        plt.ylabel('r')
        plt.title(float_name +
                  " potential conductivity (mmho/cm) multiplicative correction r with errors")

        # vertically averaged salinity plot
        plt.subplot(212)
        plt.figure(1)

        plt.plot(profile_no[0], avg_sal_offset, color=(0, 1, 0), linestyle='-', linewidth=1)
        good = np.argwhere(np.isfinite(avg_sta_offset))
        plt.plot(profile_no[0, good[:, 0]], avg_sta_offset[good[:, 0]],
                 color=(1, 0, 0), linestyle='-', linewidth=1, label='1-1 profile fit')
        plt.errorbar(profile_no[0], avg_sal_offset, yerr=2 * avg_sal_offset_err,
                     color=(0, 0, 1), linestyle='', capsize=2, label='2 x cal error')
        plt.errorbar(profile_no[0], avg_sal_offset, yerr=avg_sal_offset_err,
                     color=(0, 1, 0), linestyle='', capsize=2, label='1 x cal error')

        plt.legend()
        plt.ylabel(r'$\Delta$ S')
        plt.xlabel("Profile number")
        plt.title(float_name +
                  r" vertically averaged salinity (PSS-78) additive " +
                  r"correction $\Delta$ S with errors")

        plt.show()
