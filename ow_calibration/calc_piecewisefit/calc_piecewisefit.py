"""
-----Calculate Piece wise fit-----

Written by: Annie Wong
When: xx/10/2008
Converted to python by: Edward Small
When: 16/06/2020

Calculate the fit of each break and calibrate salinities

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import numpy as np
import scipy.io as scipy
import scipy.interpolate as interpolate
from ow_calibration.tbase_decoder.tbase_decoder import get_topo_grid


def calc_piecewisefit(float_dir, float_name, system_config):
    """
    calibrate salinities
    :param float_dir: float directory name
    :param float_name: name of float
    :param system_config: configuration parameter set up
    :return: Nothing, save output
    """

    # load in the source data

    float_source_data = scipy.loadmat(system_config['FLOAT_SOURCE_DIRECTORY'] +
                                      float_dir + float_name + system_config['FLOAT_SOURCE_POSTFIX'])

    lat = float_source_data['LAT']
    long = float_source_data['LONG']
    dates = float_source_data['DATES']
    sal = float_source_data['SAL']
    ptmp = float_source_data['PTMP']
    pres = float_source_data['PTMP']
    profile_no = float_source_data['PROFILE_NO']
    x_in = np.tile(profile_no, (10, 1))

    # load in the mapped data
    float_mapped_data = scipy.loadmat(system_config['FLOAT_MAPPED_DIRECTORY'] +
                                      float_dir + system_config['FLOAT_MAPPED_PREFIX'] +
                                      float_name + system_config['FLOAT_MAPPED_POSTFIX'])

    mapped_sal = float_mapped_data['la_mapped_sal']
    mapsalerror = float_mapped_data['la_mapsalerrors']
    mapped_ptmp = float_mapped_data['la_ptmp']
    selected_hist = float_mapped_data['selected_hist']

    # retrieve XYZ of float position used to build covariance
    if selected_hist.__len__() > 0:
        if long.shape[0] > 1:
            long = long.flatten()

        if lat.shape[0] > 1:
            lat = lat.flatten()

        if np.any(long > 180):
            long_1 = copy.deepcopy(long) - 360

        else:
            long_1 = copy.deepcopy(long)

        elev, x_grid, y_grid = get_topo_grid(np.nanmin(long_1) - 1, np.nanmax(long_1) + 1,
                                             np.nanmin(lat) - 1, np.nanmax(lat) + 1)

        grid_interp = interpolate.interp2d(x_grid[0, :], y_grid[:, 0],
                                           elev, kind='linear')

        z_grid = []
        for i in range(long_1[0].__len__()):
            z_grid.append(grid_interp(long_1[0][i], lat[0][i]))

        z_grid = -np.array(z_grid)
        coord_float = np.column_stack((long.T, lat.T, z_grid))

        # load the calibration settings

        float_calseries = scipy.loadmat(system_config['FLOAT_CALIB_DIRECTORY'] + float_dir +
                                        system_config['FLOAT_CALSERIES_PREFIX'] + float_name +
                                        system_config['FLOAT_MAPPED_POSTFIX'])

        calseries = float_calseries['calseries']
        max_breaks = float_calseries['max_breaks']
        breaks = float_calseries['breaks']
        use_theta_gt = float_calseries['use_theta_gt']
        use_theta_lt = float_calseries['use_theta_lt']
        use_pres_gt = float_calseries['use_pres_gt']
        use_pres_lt = float_calseries['use_pres_lt']
        use_percent_gt = float_calseries['use_percent_gt']

        m, n = pres.shape

        cal_sal = np.ones((m, n)) * np.nan
        cal_sal_err = np.ones((m, n)) * np.nan
        cal_cond = np.ones((m, n)) * np.nan
        cal_cond_err = np.ones((m, n)) * np.nan
        pcond_factor = np.ones((1, n)) * np.nan
        pcond_factor_err = np.ones((1, n)) * np.nan
        time_deriv = np.ones((1, n)) * np.nan
        time_deriv_err = np.ones((1, n)) * np.nan
        sta_mean = np.ones((1, n)) * np.nan
        sta_rms = np.ones((1, n)) * np.nan
        sta_sal = np.ones((m, n)) * np.nan
        sta_sal_err = np.ones((m, n)) * np.nan
        fceof = []
        fbreaks = []
