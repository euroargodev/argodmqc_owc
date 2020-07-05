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
def cal_sal_curve_plot(sal, cal_sal, cal_sal_err, sta_sal, sta_sal_err,
                       pcond_factor, pcond_factor_err, float_name):
    """
    Create the calibrated salinity curve plot
    :param sta_sal_err: error in mean difference
    :param cal_sal_err: error in calibrated salinity
    :param sal: float salinity
    :param cal_sal: calibrated salinity
    :param sta_sal: mean difference between salinity and calculated salinity
    :param pcond_factor: slope
    :param pcond_factor_err: slope error
    :param float_name: name of the float
    :return: Nothing
    """

    # only plot this if we have calibrated salinity
    if np.argwhere(~np.isnan(cal_sal)).__len__() > 0:
        # set up matrices
        profiles = sal.shape[1]
        sal_offset = cal_sal-sal
        avg_sal_offset = np.nan * np.ones((1, profiles)).flatten()
        avg_sal_offset_err = np.nan * np.ones((1, profiles)).flatten()
        sta_offset = sta_sal-sal
        avg_sta_offset = np.nan * np.ones((1, profiles)).flatten()
        avg_sta_offset_err = np.nan * np.ones((1, profiles)).flatten()

        for i in range(profiles):
            avg_sal_offset[i] = np.nanmean(sal_offset[:, i])
            avg_sal_offset_err[i] = np.nanmean(cal_sal_err[:, i])
            avg_sta_offset[i] = np.nanmean(sta_offset[: ,i])
            avg_sta_offset_err[i] = np.nanmean(sta_sal_err[:, i])



