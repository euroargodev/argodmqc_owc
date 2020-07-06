"""
-----plot diagnostics----

Written by: Breck Owens
When: 14/06/2011
 Converted to python by: Edward Small
When: 09/05/2020

Master function for loading and plotting all of the analysed data

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import scipy.io as scipy
import matplotlib.pyplot as plt
import numpy as np
from ow_calibration.find_10thetas.find_10thetas import find_10thetas
from ow_calibration.plot_diagnostics.theta_sal_plot.theta_sal_plot import theta_sal_plot
from ow_calibration.plot_diagnostics.trajectory_plot.trajectory_plot import trajectory_plot
from ow_calibration.plot_diagnostics.trajectory_plot.create_dataframe import create_dataframe
from ow_calibration.plot_diagnostics.cal_sal_curve_plot.cal_sal_curve_plot import cal_sal_curve_plot


def plot_diagnostics(float_dir, float_name, config):
    """
    run the plotting procedures
    :param float_dir: location of float source
    :param float_name: name of the float source
    :param config: user configuration dictionary
    :return: nothing, but will save the plots as PDFs
    """

    grid_data_loc = config['FLOAT_MAPPED_DIRECTORY'] + config['FLOAT_MAPPED_PREFIX'] + \
                    float_name + config['FLOAT_MAPPED_POSTFIX']
    float_data_loc = config['FLOAT_SOURCE_DIRECTORY'] + float_dir + \
                     float_name + config['FLOAT_SOURCE_POSTFIX']
    cal_data_loc = config['FLOAT_CALIB_DIRECTORY'] + float_dir + config['FLOAT_CALIB_PREFIX'] + \
                   float_name + config['FLOAT_SOURCE_POSTFIX']
    cal_series_loc = config['FLOAT_CALIB_DIRECTORY'] + float_dir + "calseries_" + \
                     float_name + config['FLOAT_SOURCE_POSTFIX']

    grid_data = scipy.loadmat(grid_data_loc)
    float_data = scipy.loadmat(float_data_loc)
    cal_data = scipy.loadmat(cal_data_loc)
    cal_series = scipy.loadmat(cal_series_loc)

    # create trajectory plot ------------------------------
    grid, floats = create_dataframe(grid_data, float_data)

    trajectory_plot(0, 0, floats, grid, float_name)

    plt.show()

    # create uncalibrated theta_s curve plot ---------------
    sal = np.array(float_data['SAL'])
    ptmp = np.array(float_data['PTMP'])
    pres = float_data['PRES']
    map_sal = grid_data['la_mapped_sal']
    map_ptmp = grid_data['la_ptmp']
    map_errors = grid_data['la_mapsalerrors']
    use_theta_lt = cal_series['use_theta_lt'][0][0]
    use_theta_gt = cal_series['use_theta_gt'][0][0]
    use_pres_lt = cal_series['use_pres_lt'][0][0]
    use_pres_gt = cal_series['use_pres_gt'][0][0]
    use_percent_gt = cal_series['use_percent_gt'][0][0]

    thetas = find_10thetas(copy.deepcopy(sal),
                           copy.deepcopy(ptmp),
                           copy.deepcopy(pres),
                           copy.deepcopy(map_ptmp),
                           use_theta_lt, use_theta_gt,
                           use_pres_lt, use_pres_gt,
                           use_percent_gt)

    index = thetas[2]

    theta_sal_plot(copy.deepcopy(sal).transpose(),
                   copy.deepcopy(ptmp).transpose(),
                   map_sal, map_ptmp, map_errors, index)

    # plot the calibration curve --------------------------

    cal_sal = cal_data['cal_SAL']
    cal_sal_err = cal_data['cal_SAL_err']
    sta_sal = cal_data['sta_SAL']
    sta_sal_err = cal_data['sta_SAL_err']
    sta_mean = cal_data['sta_mean']
    pcond_factor = cal_data['pcond_factor']
    pcond_factor_err = cal_data['pcond_factor_err']
    profile_no = float_data['PROFILE_NO']

    cal_sal_curve_plot(copy.deepcopy(sal), copy.deepcopy(cal_sal),
                       copy.deepcopy(cal_sal_err), sta_sal,
                       sta_sal_err, sta_mean, pcond_factor,
                       pcond_factor_err, profile_no, float_name)
