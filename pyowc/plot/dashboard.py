""" Functions to create plot dashboards

        Parameters
        ----------

        Returns
        -------
"""
import copy
from scipy.io import loadmat
import matplotlib.pyplot as plt
import numpy as np

from pyowc import core
from pyowc.plot import plots
from pyowc.plot.utils import create_dataframe


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

    grid_data_loc = config['FLOAT_MAPPED_DIRECTORY'] + config['FLOAT_MAPPED_PREFIX'] + \
                    float_name + config['FLOAT_MAPPED_POSTFIX']
    float_data_loc = config['FLOAT_SOURCE_DIRECTORY'] + float_dir + \
                     float_name + config['FLOAT_SOURCE_POSTFIX']
    cal_data_loc = config['FLOAT_CALIB_DIRECTORY'] + float_dir + config['FLOAT_CALIB_PREFIX'] + \
                   float_name + config['FLOAT_SOURCE_POSTFIX']
    cal_series_loc = config['FLOAT_CALIB_DIRECTORY'] + float_dir + "calseries_" + \
                     float_name + config['FLOAT_SOURCE_POSTFIX']

    grid_data = loadmat(grid_data_loc)
    float_data = loadmat(float_data_loc)
    cal_data = loadmat(cal_data_loc)
    cal_series = loadmat(cal_series_loc)

    # create trajectory plot ------------------------------
    grid, floats = create_dataframe(grid_data, float_data)

    plots.trajectory_plot(0, 0, floats, grid, float_name, config)

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

    thetas = core.finders.find_10thetas(copy.deepcopy(sal), copy.deepcopy(ptmp), copy.deepcopy(pres),
                           copy.deepcopy(map_ptmp), use_theta_lt, use_theta_gt,
                           use_pres_lt, use_pres_gt, use_percent_gt)

    index = thetas[2]

    plots.theta_sal_plot(copy.deepcopy(sal).transpose(),
                   copy.deepcopy(ptmp).transpose(),
                   map_sal, map_ptmp, map_errors, index)

    # plot the calibration curve --------------------------

    cal_sal = cal_data['cal_SAL']
    sta_sal = cal_data['sta_SAL']
    sta_sal_err = cal_data['sta_SAL_err']
    cal_sal_err = cal_data['cal_SAL_err']
    sta_mean = cal_data['sta_mean']
    pcond_factor = cal_data['pcond_factor']
    pcond_factor_err = cal_data['pcond_factor_err']
    profile_no = float_data['PROFILE_NO']

    plots.cal_sal_curve_plot(copy.deepcopy(sal), copy.deepcopy(cal_sal),
                       copy.deepcopy(cal_sal_err), sta_sal,
                       sta_sal_err, sta_mean, pcond_factor,
                       pcond_factor_err, profile_no, float_name)

    # plot the calibrated theta-S curve from float ----------

    plots.theta_sal_plot(copy.deepcopy(cal_sal).transpose(),
                   copy.deepcopy(ptmp).transpose(),
                   map_sal, map_ptmp, map_errors, index, "calibrated")

    # plot the salinity time series on theta levels ----------

    boundaries = [use_theta_lt, use_theta_gt,
                  use_pres_lt, use_pres_gt,
                  use_percent_gt]

    plots.sal_var_plot(levels, copy.deepcopy(sal), copy.deepcopy(pres),
                 copy.deepcopy(ptmp), copy.deepcopy(map_sal),
                 copy.deepcopy(map_errors), copy.deepcopy(map_ptmp),
                 copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_err),
                 boundaries, profile_no, float_name)

    # plot the analysis plots ----------------------------------

    sal_var = thetas[3]
    theta_levels = thetas[4]
    tlevels = thetas[0]
    plevels = thetas[1]

    plots.t_s_profile_plot(sal, ptmp, pres, sal_var,
                     theta_levels, tlevels, plevels, float_name)
