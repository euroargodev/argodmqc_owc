""" Functions to create plot dashboards
"""
import os
import copy
from scipy.io import loadmat
import matplotlib.pyplot as plt
import numpy as np

from pyowc.plot.plots import theta_sal_plot, cal_sal_curve_plot, sal_var_plot, t_s_profile_plot, trajectory_plot, sal_anom_plot
from pyowc.plot.utils import create_dataframe
from pyowc.core.finders import find_10thetas


#pylint: disable=too-many-locals
def plot_diagnostics(float_dir, float_name, config, levels=2):
    """ Run the plotting procedures

        Parameters
        ----------
        levels: number of theta level plots wanted (max 10)
        float_dir: location of float source
        float_name: name of the float source
        config: user configuration dictionary

        Returns
        -------
        Nothing, but will save the plots as PDFs
    """

    grid_data_loc = os.path.sep.join([config['FLOAT_MAPPED_DIRECTORY'],
                                      config['FLOAT_MAPPED_PREFIX'] + float_name + config['FLOAT_MAPPED_POSTFIX']])
    float_data_loc = os.path.sep.join([config['FLOAT_SOURCE_DIRECTORY'], float_dir,
                                       float_name + config['FLOAT_SOURCE_POSTFIX']])
    cal_data_loc = os.path.sep.join([config['FLOAT_CALIB_DIRECTORY'], float_dir,
                                     config['FLOAT_CALIB_PREFIX'] + float_name + config['FLOAT_SOURCE_POSTFIX']])
    cal_series_loc = os.path.sep.join([config['FLOAT_CALIB_DIRECTORY'], float_dir,
                                       config['FLOAT_CALSERIES_PREFIX'] + float_name + config['FLOAT_SOURCE_POSTFIX']])

    grid_data = loadmat(grid_data_loc)
    float_data = loadmat(float_data_loc)
    cal_data = loadmat(cal_data_loc)
    cal_series = loadmat(cal_series_loc)

    longi = np.transpose(float_data['LONG'])

    float_long = []
    for i in longi:
        if i > 180:
            i = i - 360
            float_long.append(i)
        else:
            float_long.append(i)

    float_long = np.array(float_long)
    float_data['LONG'] = np.transpose(float_long)
    long_i = grid_data['selected_hist']

    grid_long = []
    for i in long_i[:, 0]:
        if i > 180:
            i = i - 360
            grid_long.append(i)
        else:
            grid_long.append(i)

    long_i[:, 0] = np.transpose(np.array(grid_long))

    # create trajectory plot ------------------------------
    grid, floats = create_dataframe(grid_data, float_data)

    trajectory_plot(1, 0, floats, grid, float_name, config)
    plt.show()

    # get data ---------------
    sal = np.array(float_data['SAL'])
    ptmp = np.array(float_data['PTMP'])
    pres = float_data['PRES']
    map_sal = grid_data['la_mapped_sal']
    map_ptmp = grid_data['la_ptmp']
    map_errors = grid_data['la_mapsalerrors']

    use_theta_lt = cal_series['use_theta_lt']
    use_theta_gt = cal_series['use_theta_gt']
    use_pres_lt = cal_series['use_pres_lt']
    use_pres_gt = cal_series['use_pres_gt']
    use_percent_gt = cal_series['use_percent_gt']
    profile_no = float_data['PROFILE_NO']

    # create uncalibrated theta_s curve plot ---------------

    thetas = find_10thetas(copy.deepcopy(sal), copy.deepcopy(ptmp), copy.deepcopy(pres),
                           copy.deepcopy(map_ptmp), use_theta_lt, use_theta_gt,
                           use_pres_lt, use_pres_gt, use_percent_gt)

    index = thetas[2]

    theta_sal_plot(copy.deepcopy(sal).transpose(),
                   copy.deepcopy(ptmp).transpose(),
                   map_sal, map_ptmp, map_errors, index, profile_no[0],
                   config, float_name)

    # create the calibrated salinity anomaly plot for float

    sal_anom_plot(copy.deepcopy(sal), copy.deepcopy(ptmp), profile_no, config, float_name, "uncalibrated")

    # plot the calibration curve --------------------------

    cal_sal = cal_data['cal_SAL']
    sta_sal = cal_data['sta_SAL']
    sta_sal_err = cal_data['sta_SAL_err']
    cal_sal_err = cal_data['cal_SAL_err']
    sta_mean = cal_data['sta_mean']
    pcond_factor = cal_data['pcond_factor']
    pcond_factor_err = cal_data['pcond_factor_err']

    cal_sal_curve_plot(copy.deepcopy(sal), copy.deepcopy(cal_sal),
                       copy.deepcopy(cal_sal_err), sta_sal,
                       sta_sal_err, sta_mean, pcond_factor,
                       pcond_factor_err, profile_no, float_name, config)

    # plot the calibrated theta-S curve from float ----------

    theta_sal_plot(copy.deepcopy(cal_sal).transpose(),
                   copy.deepcopy(ptmp).transpose(),
                   map_sal, map_ptmp, map_errors, index,
                   profile_no[0], config, float_name, "calibrated")

    # plot the salinity time series on theta levels ----------

    boundaries = [use_theta_lt, use_theta_gt,
                  use_pres_lt, use_pres_gt,
                  use_percent_gt]

    sal_var_plot(levels, copy.deepcopy(sal), copy.deepcopy(pres),
                 copy.deepcopy(ptmp), copy.deepcopy(map_sal),
                 copy.deepcopy(map_errors), copy.deepcopy(map_ptmp),
                 copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_err),
                 boundaries, profile_no, float_name, config)

    # create the calibrated salinity anomaly plot for float

    sal_anom_plot(copy.deepcopy(cal_sal), copy.deepcopy(ptmp), profile_no, config, float_name, "calibrated")

    # plot the analysis plots ----------------------------------

    sal_var = thetas[3]
    theta_levels = thetas[4]
    tlevels = thetas[0]
    plevels = thetas[1]

    t_s_profile_plot(sal, ptmp, pres, sal_var,
                     theta_levels, tlevels, plevels,
                     float_name, config)
