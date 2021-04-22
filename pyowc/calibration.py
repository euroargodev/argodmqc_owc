""" Functions to calibrate data

"""
import os
import time
import copy

import gsw
import numpy as np

from scipy.io import loadmat, savemat
import scipy.interpolate as interpolate

from .core.stats import signal_variance, noise_variance, build_cov, fit_cond
from .core.finders import find_besthist, find_25boxes, find_10thetas
from .data.wrangling import interp_climatology, map_data_grid
from .data.fetchers import get_topo_grid, get_region_data, get_region_hist_locations, frontal_constraint_saf
from .helper import load_varibales_from_file, process_profiles_la_variables, process_profiles_grid_variables, \
    get_float_data, process_profile_hist_variables, remove_statical_outliers, sort_numpy_array, \
    check_and_make_numpy_arry, selected_historical_points


# pylint: disable=too-many-lines

# pylint: disable=too-many-locals
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
# pylint: disable=too-many-nested-blocks
# pylint: disable=invalid-name
# pylint: disable=fixme
def update_salinity_mapping(float_dir, config, float_name):
    """ Calculates values needed for analysis and save on file

        Cecile Cabanes, June. 2013: use of "map_large_scale" (time scale) used
        to map the large scale field

        Function used to store all the resulting variables needed to run the analysis
        for each float profile.

        Also saves all the settings used for the analysis, in case an operator wishes
        to either run the analysis again, or share the analysis with another operator.

        Parameters
        ----------
        float_dir: directory where the float being analysed is stored
        float_name: name of the float being analysed
        config: configuration settings set by the user

        Returns
        -------
        Nothing, but does save values
    """

    # Load float source data ---------------------------------------------

    # Get float file name
    filename = os.path.sep.join([config['FLOAT_SOURCE_DIRECTORY'], float_dir,
                                 float_name + config['FLOAT_SOURCE_POSTFIX']])

    # Load the float data
    float_source_data = loadmat(filename)

    # Get the profile number and size of the data in the profile
    profile_no = float_source_data['PROFILE_NO'][0]
    float_level_count = float_source_data['SAL'].shape[0]

    # Load all the mapping parameters, including the WMO boxes -----------

    wmo_boxes = loadmat(os.path.sep.join([config['CONFIG_DIRECTORY'], config['CONFIG_WMO_BOXES']]))
    max_casts = config['CONFIG_MAX_CASTS']
    map_use_pv = config['MAP_USE_PV']
    map_use_saf = config['MAP_USE_SAF']
    long_large = config['MAPSCALE_LONGITUDE_LARGE']
    long_small = config['MAPSCALE_LONGITUDE_SMALL']
    lat_large = config['MAPSCALE_LATITUDE_LARGE']
    lat_small = config['MAPSCALE_LATITUDE_SMALL']
    phi_large = config['MAPSCALE_PHI_LARGE']
    phi_small = config['MAPSCALE_PHI_SMALL']
    map_age_large = config['MAPSCALE_AGE_LARGE']
    map_age_small = config['MAPSCALE_AGE_SMALL']
    map_p_delta = config['MAP_P_DELTA']
    map_p_exclude = config['MAP_P_EXCLUDE']

    # Display the configuration to the user ------------------------------

    print("\nCONFIGURATION PARAMETERS")
    print("__________________________________________________________")
    print("WMO box file: ", config['CONFIG_WMO_BOXES'])
    print("max_casts: ", max_casts)
    print("map_use_pv: ", map_use_pv)
    print("map_use_saf: ", map_use_saf)
    print("long_large: ", long_large)
    print("long_small: ", long_small)
    print("lat_large: ", lat_large)
    print("lat_small: ", lat_small)
    print("phi_large: ", phi_large)
    print("phi_small: ", phi_small)
    print("map_age_large: ", map_age_large)
    print("map_age_small: ", map_age_small)
    print("map_p_delta: ", map_p_delta)
    print("map_p_exclude: ", map_p_exclude)
    print("__________________________________________________________")

    # Load precalculated mapped data -------------------------------------

    # Check to see if we have any precalculated mapped data
    mapped_data_path = os.path.sep.join([config['FLOAT_MAPPED_DIRECTORY'], float_dir,
                                         config['FLOAT_MAPPED_PREFIX'] +
                                         float_name +
                                         config['FLOAT_MAPPED_POSTFIX']])
    mapped_data_path = os.path.abspath(mapped_data_path)

    # load data from file
    data = load_varibales_from_file(mapped_data_path, float_level_count)

    # Compare profile numbers in the float source against the mapped data matrix -
    profile_index = data['profile_index']
    missing_profile_index = []

    for i in range(0, profile_no.__len__()):
        profiles = np.argwhere(data['la_profile_no'] == profile_no[i])
        if profiles.size == 0:
            missing_profile_index.append(i)

    # update mapped data with missing profiles ---------------------------

    for i in range(0, missing_profile_index.__len__()):
        # start the timer
        start_time = time.time()

        print("UPDATE_SALINITY_MAPPING: Working on profile ", i)

        # Get the current profile being worked on
        missing_profile = missing_profile_index[i]

        # append profile numbers
        data['la_profile_no'] = np.insert(data['la_profile_no'], profile_index, profile_no[missing_profile])
        # Construct elements for this profile

        # initialise matrices to hold and save parameter settings
        data['scale_long_large'].append(np.nan)
        data['scale_lat_large'].append(np.nan)
        data['scale_long_small'].append(np.nan)
        data['scale_lat_small'].append(np.nan)
        data['scale_phi_large'].append(np.nan)
        data['scale_phi_small'].append(np.nan)
        data['scale_age_large'].append(np.nan)
        data['scale_age_small'].append(np.nan)
        data['use_pv'].append(np.nan)
        data['use_saf'].append(np.nan)
        data['p_delta'].append(np.nan)
        data['p_exclude'].append(np.nan)

        # helper method 1
        data = process_profiles_la_variables(data, float_level_count, profile_index)

        # get data from float
        float_data = get_float_data(float_source_data, missing_profile)

        # if we have any good location data and pressure data from the float
        if not np.isnan(float_data['float_long']) and not np.isnan(float_data['float_lat']) \
                and np.argwhere(np.isnan(float_data['float_pres']) == 0).any():

            # tbase.int file requires longitudes from 0 to +/-180
            float_long_tbase = copy.deepcopy(float_data['float_long'])

            if float_long_tbase > 180:
                float_long_tbase -= 360

            # find the depth of the ocean at the float location
            float_elev, float_x, float_y = get_topo_grid(float_long_tbase - 1,
                                                         float_long_tbase + 1,
                                                         float_data['float_lat'] - 1,
                                                         float_data['float_lat'] + 1,
                                                         config)

            float_interp = interpolate.interp2d(float_x[0, :],
                                                float_y[:, 0],
                                                float_elev,
                                                kind='linear')

            float_data['float_z'] = -float_interp(float_long_tbase, float_data['float_lat'])[0]

            # gather data from area surrounding the float location
            wmo_numbers = find_25boxes(float_data['float_long'], float_data['float_lat'], wmo_boxes)

            grid_data = {}
            grid_lat, grid_long, grid_dates = get_region_hist_locations(wmo_numbers, float_name, config)

            grid_data['grid_lat'] = grid_lat
            grid_data['grid_long'] = grid_long
            grid_data['grid_dates'] = grid_dates

            # if we have data in the surrounding area, find depths at these points
            if grid_lat.__len__() > 0:

                grid_data = process_profiles_grid_variables(grid_data, config)

                # make sure that the grid and float longitudes match at the 0-360 mark
                float_long_0 = float_data['float_long']
                if np.argwhere(grid_long > 360).__len__() > 0:
                    if 0 <= float_data['float_long'] <= 20:
                        float_long_0 += 360

                index = find_besthist(grid_data['grid_lat'], grid_data['grid_long'], grid_data['grid_dates'],
                                      grid_data['grid_z'], float_data['float_lat'], float_long_0,
                                      float_data['float_date'], float_data['float_z'],
                                      lat_large, lat_small, long_large, long_small,
                                      phi_large, phi_small, map_age_large, map_age_small,
                                      map_use_pv, max_casts)

                # Now that we have the indices of the best spatial and temporal data
                # we can get the rest of the data from these casts
                best_hist_data = get_region_data(wmo_numbers, float_name, config, index, float_data['float_pres'])

                best_hist_data['grid_z'] = grid_data['grid_z'][index]

                # If we are near the Subantarctic Front we need to figure out if
                # the profile is north or south of it. Then we should remove data not on
                # the same side the profile is on
                if map_use_saf == 1:
                    best_hist_data_2 = frontal_constraint_saf(config, best_hist_data, float_data)

                    # Use frontal separation only if there are at least 5 profiles
                    if len(best_hist_data_2['grid_sal'][0]) > 5:
                        best_hist_data.update(best_hist_data_2)

                # make the float longitude wrap around the 0-360 mark if the historical data has
                if np.argwhere(best_hist_data['grid_long'] > 360).__len__() > 0:
                    if 0 <= float_data['float_long'] <= 20:
                        float_data['float_long'] += 360

                # interpolate historical data onto float theta levels
                hist_interp_sal, hist_interp_pres = interp_climatology(best_hist_data['grid_sal'],
                                                                       best_hist_data['grid_ptmp'],
                                                                       best_hist_data['grid_pres'],
                                                                       float_data['float_sal'],
                                                                       float_data['float_ptmp'],
                                                                       float_data['float_pres'])

                # map float theta levels below the exclusion point

                for n_level in range(float_level_count):

                    # for each float theta level, only use interpolated salinity
                    # that aren't NaN's
                    if not np.isnan(float_data['float_sal'][n_level]) \
                            and float_data['float_pres'][n_level] >= map_p_exclude:

                        hist_data = process_profile_hist_variables(grid_data=best_hist_data,
                                                                   float_pres=float_data['float_pres'],
                                                                   hist_interp_sal=hist_interp_sal,
                                                                   hist_interp_pres=hist_interp_pres,
                                                                   n_level=n_level,
                                                                   map_p_delta=map_p_delta
                                                                   )

                        # only proceed with analysis if we have more than 5 points
                        if hist_data['hist_sal'].__len__() > 5:

                            # Need to check for statistical outliers
                            mean_sal = np.mean(hist_data['hist_sal'])
                            signal_sal = signal_variance(hist_data['hist_sal'])
                            outlier1 = np.argwhere(np.abs(hist_data['hist_sal'] - mean_sal) /
                                                   np.sqrt(signal_sal) > 3)
                            outlier = outlier1[:, 0]

                            # remove the statical outliers
                            hist_data = remove_statical_outliers(outlier, hist_data)

                            # calculate signal and noise for complete data
                            hist_sal_flatten = hist_data['hist_sal'].flatten()
                            noise_sal = noise_variance(hist_sal_flatten,
                                                       hist_data['hist_lat'].flatten(),
                                                       hist_data['hist_long'].flatten())

                            signal_sal = signal_variance(hist_data['hist_sal'])

                            # map residuals

                            float_data_array = np.array([float_data['float_lat'], float_data['float_long'],
                                                         float_data['float_date'], float_data['float_z']]
                                                        ).reshape((-1, 4))

                            hist_data_array = np.column_stack((hist_data['hist_lat'], hist_data['hist_long'],
                                                               hist_data['hist_dates'], hist_data['hist_z']))

                            mapped_values = map_data_grid(hist_sal_flatten,
                                                          float_data_array, hist_data_array, long_large, lat_large,
                                                          map_age_large, signal_sal, noise_sal, phi_large, map_use_pv)

                            # use short scales to map the residuals

                            sal_residual = hist_sal_flatten - mapped_values[2]
                            sal_signal_residual = signal_variance(sal_residual)

                            mapped_residuals = map_data_grid(sal_residual,
                                                             float_data_array,
                                                             hist_data_array,
                                                             long_small, lat_small,
                                                             map_age_small,
                                                             sal_signal_residual, noise_sal,
                                                             phi_small, map_use_pv)

                            data['la_ptmp'][n_level, profile_index] = float_data['float_ptmp'][n_level]
                            data['la_mapped_sal'][n_level, profile_index] = mapped_values[0] + \
                                                                            mapped_residuals[0]

                            dot_map_values = np.dot(mapped_values[1], mapped_values[1])
                            data['la_map_sal_errors'][n_level, profile_index] = np.sqrt(dot_map_values +
                                                                                        mapped_residuals[1] *

                                                                                        mapped_residuals[1])
                            data['la_noise_sal'][n_level, profile_index] = noise_sal
                            data['la_signal_sal'][n_level, profile_index] = signal_sal
                            data['scale_long_large'][profile_index] = long_large
                            data['scale_lat_large'][profile_index] = lat_large
                            data['scale_long_small'][profile_index] = long_small
                            data['scale_lat_small'][profile_index] = lat_small
                            data['scale_phi_large'][profile_index] = phi_large
                            data['scale_phi_small'][profile_index] = phi_small
                            data['scale_age_large'][profile_index] = map_age_large
                            data['scale_age_small'][profile_index] = map_age_small
                            data['use_pv'][profile_index] = map_use_pv
                            data['use_saf'][profile_index] = map_use_saf
                            data['p_delta'][profile_index] = map_p_delta
                            data['p_exclude'][profile_index] = map_p_exclude

                            # only save selected historical points
                            data = selected_historical_points(data, hist_data, profile_index)

        print("time elapsed: ", round(time.time() - start_time, 2), " seconds")
        profile_index += 1

    # as a quality control check, just make sure salinities are between 30 and 40
    bad_sal_30 = np.argwhere(data['la_mapped_sal'] < 30)
    bad_sal_40 = np.argwhere(data['la_mapped_sal'] > 40)
    bad_sal = np.concatenate((bad_sal_30, bad_sal_40))

    for sal in bad_sal:
        data['la_mapped_sal'][sal[0], sal[1]] = np.nan

    # sort the data by profile number
    sorted_profile_index = data['la_profile_no'].argsort()

    # finalise the data sorted by profile number, make data type numpy arrays
    data = sort_numpy_array(data, sorted_profile_index, ['la_ptmp', 'la_mapped_sal',
                                                         'la_map_sal_errors', 'la_noise_sal',
                                                         'la_signal_sal'])
    # check is numpy array
    data = check_and_make_numpy_arry(data)

    # sort value in dic
    data = sort_numpy_array(data, sorted_profile_index)

    if data['selected_hist'].__len__() > 0:
        selected_hist_1 = np.array(sorted(data['selected_hist'], key=lambda x: x[2]))

    data['selected_hist'] = selected_hist_1

    # define the saving location
    save_location = os.path.sep.join([config['FLOAT_MAPPED_DIRECTORY'],
                                      config['FLOAT_MAPPED_PREFIX'] + float_name + config['FLOAT_MAPPED_POSTFIX']])

    # save the data
    savemat(save_location, data)


# pylint: disable=invalid-name
# pylint: disable=too-many-locals
# pylint: disable=unused-variable
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
def calc_piecewisefit(float_dir, float_name, system_config):
    """ Calibrate salinities

        Calculate the fit of each break and calibrate salinities

        Will with data from:
        FLOAT_SOURCE_DIRECTORY
        FLOAT_MAPPED_DIRECTORY

        Parameters
        ----------
        float_dir: float directory name
        float_name: name of float
        system_config: configuration parameter set up

        Returns
        -------
        Nothing, save output

    """

    # load in the source data
    float_source_data_path = os.path.sep.join([system_config['FLOAT_SOURCE_DIRECTORY'], float_dir,
                                               float_name + system_config['FLOAT_SOURCE_POSTFIX']])
    float_source_data = loadmat(float_source_data_path)

    lat = float_source_data['LAT']
    long = float_source_data['LONG']
    sal = float_source_data['SAL']
    ptmp = float_source_data['PTMP']
    pres = float_source_data['PRES']
    profile_no = float_source_data['PROFILE_NO']
    x_in = np.tile(profile_no, (10, 1))

    # load in the mapped data
    float_mapped_data_path = os.path.sep.join([system_config['FLOAT_MAPPED_DIRECTORY'], float_dir,
                                               system_config['FLOAT_MAPPED_PREFIX'] +
                                               float_name + system_config['FLOAT_MAPPED_POSTFIX']])
    float_mapped_data = loadmat(float_mapped_data_path)

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
                                             np.nanmin(lat) - 1, np.nanmax(lat) + 1, system_config)

        grid_interp = interpolate.interp2d(x_grid[0, :], y_grid[:, 0],
                                           elev, kind='linear')

        z_grid = []
        for i in range(long_1[0].__len__()):
            z_grid.append(grid_interp(long_1[0][i], lat[0][i]))

        z_grid = -np.array(z_grid)
        coord_float = np.column_stack((long.T, lat.T, z_grid))

    # load the calibration settings
    float_calseries_path = os.path.sep.join([system_config['FLOAT_CALIB_DIRECTORY'], float_dir,
                                             system_config['FLOAT_CALSERIES_PREFIX'] + float_name +
                                             system_config['FLOAT_MAPPED_POSTFIX']])
    float_calseries = loadmat(float_calseries_path)

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
    cal_sal1 = np.ones((m, n)) * np.nan
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

    sstatus = 1
    unique_cal = np.unique(calseries.flatten())
    # bad profiles are flagged as zero
    bad = np.argwhere(unique_cal == 0)

    if bad.__len__() > 0:
        unique_cal = np.delete(unique_cal, bad)

    n_seq = unique_cal.__len__()
    if n_seq == 1 and max_breaks.__len__() > 1:
        print("Error in specificying number of possible break points")
        print(str(max_breaks), " specified, should be ",
              str([max_breaks.__len__(), n_seq]))

    # we have multiple cal series, make sure that break information is provideed for all segments
    elif n_seq > 1:
        # only one max break specified, specify this break for all segements
        if max_breaks.__len__() == 1:
            max_breaks = np.ones((n_seq, 1)) * max_breaks

        # error in specification of max breaks
        elif max_breaks.__len__() != n_seq:
            print("Error in specifying the number of possible break points")
            print(str(max_breaks), " specified, should be 1 or ",
                  str([max_breaks.__len__(), n_seq]))
            sstatus = 0

    if breaks.__len__() > 0:
        ns, nb = breaks.shape

        # error in specifying breaks
        if ns != n_seq:
            print("Error in specifying break points")
            print("For multiple cal series, need to specify breaks for each series")
            print("Have ", str(n_seq), " or ", str(ns), " sets of breaks")

            sstatus = 0

        for n in range(n_seq):
            nb = np.argwhere(np.isfinite(breaks[n, :])).__len__()

            if nb > max_breaks[n]:
                print("Error, for cal series ", str(unique_cal[n]), "max number of breaks ",
                      str(max_breaks[n]), " less than ", str(nb), "prescribed breaks")
                sstatus = 0

            elif nb < max_breaks[n]:
                print("Specified ", str(nb), " breaks. Will search up to ",
                      str(max_breaks[n]), " breaks")

            else:
                print(str(nb), "fixed breaks prescribed")

    # set_calseries returned a bad status variable, write out file with NaNs
    if sstatus == 0:
        float_calib_filename = (os.path.sep.join([system_config['FLOAT_CALIB_DIRECTORY'], float_dir,
                                                  system_config['FLOAT_CALIB_PREFIX'] + float_name +
                                                  system_config['FLOAT_CALIB_POSTFIX']]))

        savemat(float_calib_filename,
                {'cal_SAL': cal_sal,
                 'cal_SAL_err': cal_sal_err,
                 'pcond_factor': pcond_factor,
                 'pcond_factor_err': pcond_factor_err,
                 'cal_COND': cal_cond,
                 'cal_COND_err': cal_cond_err,
                 'time_deriv': time_deriv,
                 'time_deriv_err': time_deriv_err,
                 'sta_mean': sta_mean,
                 'sta_rms': sta_rms,
                 'sta_SAL': sta_sal,
                 'sta_SAL_err': sta_sal_err,
                 'PROFILE_NO': profile_no,
                 'fcoef': fceof,
                 'fbreaks': fbreaks})

        return

    # loop through sequences of calseries

    for i in range(n_seq):
        calindex = np.argwhere(calseries == unique_cal[i])[:, 1]
        k = calindex.__len__()

        # chose 10 float theta levels to use for the piecewise linear fit
        unique_coord_float = coord_float[calindex, :]
        unique_sal = sal[:, calindex]
        unique_ptmp = ptmp[:, calindex]
        unique_pres = pres[:, calindex]
        unique_mapped_ptmp = mapped_ptmp[:, calindex]
        unique_mapped_sal = mapped_sal[:, calindex]
        unique_mapsalerrors = mapsalerror[:, calindex]

        ten_sal = np.ones((10, k)) * np.nan
        ten_ptmp = np.ones((10, k)) * np.nan
        ten_pres = np.ones((10, k)) * np.nan
        ten_mapped_sal = np.ones((10, k)) * np.nan
        ten_mapsalerrors = np.ones((10, k)) * np.nan

        # make deep copies for calibration layter
        unique_sal_1 = copy.deepcopy(unique_sal)
        unique_ptmp_1 = copy.deepcopy(unique_ptmp)
        unique_pres_1 = copy.deepcopy(unique_pres)
        unique_mapped_ptmp_1 = copy.deepcopy(unique_ptmp)
        unique_mapped_sal_1 = copy.deepcopy(unique_mapped_sal)
        unique_mapsalerrors_1 = copy.deepcopy(unique_mapsalerrors)

        theta, p, index, var_s_th, th = find_10thetas(copy.deepcopy(unique_sal),
                                                      copy.deepcopy(unique_ptmp),
                                                      copy.deepcopy(unique_pres),
                                                      copy.deepcopy(unique_mapped_ptmp),
                                                      use_theta_lt, use_theta_gt,
                                                      use_pres_lt, use_pres_gt,
                                                      use_percent_gt)

        index = np.array(index, dtype=int)
        pp = np.argwhere(np.isnan(index) == 0)
        # only proceed if we have valied theta levels
        if pp.__len__() > 0:
            for ipr in range(k):
                jj = np.argwhere(index[:, ipr] >= 0)
                if jj.__len__() > 0:
                    ten_sal[0:jj.__len__(), ipr] = unique_sal[index[jj, ipr], ipr].flatten()
                    ten_ptmp[0:jj.__len__(), ipr] = unique_ptmp[index[jj, ipr], ipr].flatten()
                    ten_pres[0:jj.__len__(), ipr] = unique_pres[index[jj, ipr], ipr].flatten()
                    ten_mapped_sal[0:jj.__len__(), ipr] = unique_mapped_sal[index[jj, ipr],
                                                                            ipr].flatten()
                    ten_mapsalerrors[0:jj.__len__(), ipr] = unique_mapsalerrors[index[jj, ipr],
                                                                                ipr].flatten()
            # calculate potential conductivites and errors for mapped values and float values
            # calculate pcond error by perturbing salinity
            # (avoids problems caused by non-linearity)

            # constant for conductivity at sal=35, temp=15 and pres=0
            sw_c3515 = 42.914

            icond = gsw.conversions.C_from_SP(ten_sal,
                                              ten_ptmp,
                                              0)
            mapped_cond = gsw.conversions.C_from_SP(ten_mapped_sal,
                                                    ten_ptmp,
                                                    0)

            mapped_cond1 = gsw.conversions.C_from_SP(ten_mapped_sal + ten_mapsalerrors / 100,
                                                     ten_ptmp, 0)

            mapconderrors = 100 * np.abs(mapped_cond - mapped_cond1)

            # independent variable for pieve wise fit (profile number)
            x = x_in[:, calindex]
            y = mapped_cond / icond
            err = mapconderrors / icond

            # calculate off-diagonal terms for error estimate

            covariance = build_cov(ten_ptmp, unique_coord_float, system_config)

            # if no break points are set
            if breaks.__len__() == 0:
                (xfit, condslope, condslope_err,
                 time_deriv_temp, time_deriv_err_temp,
                 sta_mean_temp, sta_rms_temp, ndf,
                 fit_coef, fit_breaks) = fit_cond(x, y, err,
                                                  covariance,
                                                  'max_no_breaks',
                                                  max_breaks[i][0])
                pcond_factor[0][calindex] = condslope
                pcond_factor_err[0][calindex] = condslope_err
                time_deriv[0][calindex] = time_deriv_temp.flatten()
                time_deriv_err[0][calindex] = time_deriv_err_temp.flatten()
                sta_mean[0][calindex], = sta_mean_temp
                sta_rms[0][calindex] = sta_rms_temp

            else:
                breaks_in = breaks[i, :]
                breaks_in = breaks_in[np.argwhere(np.isfinite(breaks_in))]

                if max_breaks[i]:
                    (xfit, condslope, condslope_err,
                     time_deriv_temp, time_deriv_err_temp,
                     sta_mean_temp, sta_rms_temp, ndf,
                     fit_coef, fit_breaks) = fit_cond(x, y, err,
                                                      covariance,
                                                      'breaks',
                                                      breaks_in,
                                                      'max_no_breaks',
                                                      max_breaks[i][0])
                    pcond_factor[0][calindex] = condslope
                    pcond_factor_err[0][calindex] = condslope_err
                    time_deriv[calindex] = time_deriv_temp
                    time_deriv_err[calindex] = time_deriv_err_temp
                    sta_mean[0][calindex], = sta_mean_temp
                    sta_rms[0][calindex] = sta_rms_temp

                else:
                    (xfit, condslope, condslope_err,
                     time_deriv_temp, time_deriv_err_temp,
                     sta_mean_temp, sta_rms_temp, ndf,
                     fit_coef, fit_breaks) = fit_cond(x, y, err,
                                                      covariance,
                                                      'breaks',
                                                      breaks_in)
                    pcond_factor[0][calindex] = condslope
                    pcond_factor_err[0][calindex] = condslope_err
                    time_deriv[calindex] = time_deriv_temp
                    time_deriv_err[calindex] = time_deriv_err_temp
                    sta_mean[0][calindex], = sta_mean_temp
                    sta_rms[0][calindex] = sta_rms_temp

            # apply calibrations to float data

            if pcond_factor[0][calindex].__len__() > 0:
                unique_cond = gsw.conversions.C_from_SP(unique_sal_1, unique_ptmp_1, 0)
                cal_cond[:, calindex] = np.dot(np.ones((m, 1)),
                                               pcond_factor[:, calindex]) * unique_cond
                cal_sal[:, calindex] = gsw.conversions.SP_from_C(cal_cond[:, calindex],
                                                                 unique_ptmp_1,
                                                                 0)
                cal_cond_err[:, calindex] = np.dot(np.ones((m, 1)),
                                                   pcond_factor_err[:, calindex]) * unique_cond
                cal_sal1[:, calindex] = gsw.conversions.SP_from_C((cal_cond[:, calindex] +
                                                                   cal_cond_err[:, calindex]),
                                                                  unique_ptmp, 0)

                cal_sal_err[:, calindex] = np.abs(cal_sal[:, calindex] - cal_sal1[:, calindex])

                # estimate the error in salinity for station by fit

                sta_cond = np.dot(np.ones((m, 1)), sta_mean[:, calindex]) * unique_cond
                sta_sal[:, calindex] = gsw.conversions.SP_from_C(sta_cond, unique_ptmp, 0)
                sta_cond_err = np.dot(np.ones((m, 1)), sta_rms[:, calindex]) * unique_cond
                sta_sal1 = gsw.conversions.SP_from_C(sta_cond + sta_cond_err, unique_ptmp, 0)
                sta_sal_err[:, calindex] = np.abs(sta_sal[:, calindex] - sta_sal1)

                for n in range(fit_coef.__len__()):
                    fceof.append(fit_coef[n])

                if fit_breaks.__len__() > 0:
                    fbreaks.append(fit_breaks)

    # save calibration data

    float_calib_name = (os.path.sep.join([system_config['FLOAT_CALIB_DIRECTORY'], float_dir,
                                          system_config['FLOAT_CALIB_PREFIX'] + float_name
                                          + system_config['FLOAT_CALIB_POSTFIX']]))

    savemat(float_calib_name, {'cal_SAL': cal_sal,
                               'cal_SAL_err': cal_sal_err,
                               'pcond_factor': pcond_factor,
                               'pcond_factor_err': pcond_factor_err,
                               'cal_COND': cal_cond,
                               'cal_COND_err': cal_cond_err,
                               'time_deriv': time_deriv,
                               'time_deriv_err': time_deriv_err,
                               'sta_mean': sta_mean,
                               'sta_rms': sta_rms,
                               'sta_SAL': sta_sal,
                               'sta_SAL_err': sta_sal_err,
                               'PROFILE_NO': profile_no,
                               'fcoef': fceof,
                               'fbreaks': fbreaks})
