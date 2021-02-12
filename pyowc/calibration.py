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

    try:

        # open up mapped data
        float_mapped_data = loadmat(mapped_data_path)

        # Save the data to variables
        la_mapped_sal = float_mapped_data['la_mapped_sal']
        la_mapsalerrors = float_mapped_data['la_mapsalerrors']
        la_noise_sal = float_mapped_data['la_noise_sal']
        la_signal_sal = float_mapped_data['la_signal_sal']
        la_ptmp = float_mapped_data['la_ptmp']
        la_profile_no = float_mapped_data['la_profile_no'].flatten()
        scale_long_large = float_mapped_data['scale_long_large'].flatten()
        scale_lat_large = float_mapped_data['scale_lat_large'].flatten()
        scale_long_small = float_mapped_data['scale_long_small'].flatten()
        scale_lat_small = float_mapped_data['scale_lat_small'].flatten()
        scale_phi_large = float_mapped_data['scale_phi_large'].flatten()
        scale_phi_small = float_mapped_data['scale_phi_small'].flatten()
        scale_age_large = float_mapped_data['scale_age_large'].flatten()
        scale_age_small = float_mapped_data['scale_age_small'].flatten()
        use_pv = float_mapped_data['use_pv'].flatten()
        use_saf = float_mapped_data['use_saf'].flatten()
        p_delta = float_mapped_data['p_delta'].flatten()
        p_exclude = float_mapped_data['p_exclude'].flatten()
        selected_hist = float_mapped_data['selected_hist']

        # Check to see if this is an older version run without the saf constraint
        if not "use_saf" in float_mapped_data:
            use_saf = np.zeros(float_mapped_data['use_pv'].shape)

        # Get mapped data shape
        profile_index = la_mapped_sal.shape[1]
        max_depth = la_mapped_sal.shape[0]
        how_many_cols = la_mapped_sal.shape[1]
        new_depth = float_level_count

        # if we have more data available than in the current mapped data, we need to extend
        # the matrices so we can add this data

        if new_depth > max_depth != 0:
            la_mapped_sal = np.insert(la_mapped_sal, la_mapped_sal.shape[0],
                                      np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                      axis=0)
            la_mapsalerrors = np.insert(la_mapsalerrors, la_mapsalerrors.shape[0],
                                        np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                        axis=0)
            la_noise_sal = np.insert(la_noise_sal, la_noise_sal.shape[0],
                                     np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                     axis=0)
            la_signal_sal = np.insert(la_signal_sal, la_signal_sal.shape[0],
                                      np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                      axis=0)
            la_ptmp = np.insert(la_ptmp, la_ptmp.shape[0],
                                np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                axis=0)

        print("Using precalculated data: ", mapped_data_path)
        print("__________________________________________________________")

    # If we don't have any precalculated mapped data
    except FileNotFoundError:

        # initialise variables
        profile_index = 0
        la_profile_no = np.empty(0)
        selected_hist = []
        la_ptmp = np.empty((float_level_count, 0))
        la_mapped_sal = np.empty((float_level_count, 0))
        la_mapsalerrors = np.empty((float_level_count, 0))
        la_noise_sal = np.empty((float_level_count, 0))
        la_signal_sal = np.empty((float_level_count, 0))
        scale_long_large = []
        scale_lat_large = []
        scale_long_small = []
        scale_lat_small = []
        scale_phi_large = []
        scale_phi_small = []
        scale_age_large = []
        scale_age_small = []
        use_pv = []
        use_saf = []
        p_delta = []
        p_exclude = []

        print("No precalculated data at: %s" % mapped_data_path)
        print("__________________________________________________________\n")

    # Compare profile numbers in the float source against the mapped data matrix -
    missing_profile_index = []

    for i in range(0, profile_no.__len__()):
        profiles = np.argwhere(la_profile_no == profile_no[i])
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
        la_profile_no = np.insert(la_profile_no, profile_index, profile_no[missing_profile])
        # Construct elements for this profile

        # if we are inserting changing a column in existing data
        if profile_index < la_ptmp.shape[1]:
            la_ptmp[:, profile_index] = np.nan * np.ones(float_level_count)
            la_mapped_sal[:, profile_index] = np.nan * np.ones(float_level_count)
            la_mapsalerrors[:, profile_index] = np.nan * np.ones(float_level_count)
            la_noise_sal[:, profile_index] = np.nan * np.ones(float_level_count)
            la_signal_sal[:, profile_index] = np.nan * np.ones(float_level_count)

        # if we are adding a new column
        else:
            la_ptmp = np.hstack((la_ptmp,
                                 np.nan * np.ones((float_level_count, 1))))
            la_mapped_sal = np.hstack((la_mapped_sal,
                                       np.nan * np.ones((float_level_count, 1))))
            la_mapsalerrors = np.hstack((la_mapsalerrors,
                                         np.nan * np.ones((float_level_count, 1))))
            la_noise_sal = np.hstack((la_noise_sal,
                                      np.nan * np.ones((float_level_count, 1))))
            la_signal_sal = np.hstack((la_signal_sal,
                                       np.nan * np.ones((float_level_count, 1))))

        # initialise matrices to hold and save parameter settings
        scale_long_large.append(np.nan)
        scale_lat_large.append(np.nan)
        scale_long_small.append(np.nan)
        scale_lat_small.append(np.nan)
        scale_phi_large.append(np.nan)
        scale_phi_small.append(np.nan)
        scale_age_large.append(np.nan)
        scale_age_small.append(np.nan)
        use_pv.append(np.nan)
        use_saf.append(np.nan)
        p_delta.append(np.nan)
        p_exclude.append(np.nan)

        # get data from float
        float_lat = float_source_data['LAT'][0, missing_profile]
        float_long = float_source_data['LONG'][0, missing_profile]
        float_date = float_source_data['DATES'][0, missing_profile]
        float_sal = float_source_data['SAL'][:, missing_profile]
        float_tmp = float_source_data['TEMP'][:, missing_profile]
        float_ptmp = float_source_data['PTMP'][:, missing_profile]
        float_pres = float_source_data['PRES'][:, missing_profile]

        # if we have any good location data and pressure data from the float
        if not np.isnan(float_long) \
                and not np.isnan(float_lat) \
                and np.argwhere(np.isnan(float_pres) == 0).any():

            # tbase.int file requires longitudes from 0 to +/-180
            float_long_tbase = copy.deepcopy(float_long)
            if float_long_tbase > 180:
                float_long_tbase -= 360

            # find the depth of the ocean at the float location
            float_elev, float_x, float_y = get_topo_grid(float_long_tbase - 1,
                                                         float_long_tbase + 1,
                                                         float_lat - 1,
                                                         float_lat + 1,
                                                         config)
            float_interp = interpolate.interp2d(float_x[0, :],
                                                float_y[:, 0],
                                                float_elev,
                                                kind='linear')
            float_z = -float_interp(float_long_tbase, float_lat)[0]

            # gather data from area surrounding the float location
            wmo_numbers = find_25boxes(float_long, float_lat, wmo_boxes)
            grid_lat, grid_long, grid_dates = get_region_hist_locations(wmo_numbers, float_name, config)

            # if we have data in the surrounding area, find depths at these points
            if grid_lat.__len__() > 0:
                # tbase.int file requires longitudes from 0 to +/-180
                grid_long_tbase = copy.deepcopy(grid_long)
                g_180 = np.argwhere(grid_long_tbase > 180)
                grid_long_tbase[g_180] -= 360

                # find depth of the ocean at historical locations
                grid_elev, grid_x, grid_y = get_topo_grid(np.amin(grid_long_tbase) - 1,
                                                          np.amax(grid_long_tbase) + 1,
                                                          np.amin(grid_lat) - 1,
                                                          np.amax(grid_lat) + 1,
                                                          config)

                grid_interp = interpolate.interp2d(grid_x[0], grid_y[:, 0],
                                                   grid_elev, kind='linear')

                # As a note, the reason we vectorise the function here is because we do not
                # want to compare every longitudinal value to ever latitude. Rather, we simply
                # want to interpolate each pair of longitudes and latitudes.

                grid_z = -1 * np.vectorize(grid_interp)(grid_long_tbase, grid_lat)

                # make sure that the grid and float longitudes match at the 0-360 mark
                float_long_0 = float_long
                if np.argwhere(grid_long > 360).__len__() > 0:
                    if 0 <= float_long <= 20:
                        float_long_0 += 360

                index = find_besthist(grid_lat, grid_long, grid_dates, grid_z,
                                      float_lat, float_long_0, float_date, float_z,
                                      lat_large, lat_small, long_large, long_small,
                                      phi_large, phi_small, map_age_large, map_age_small,
                                      map_use_pv, max_casts)

                # Now that we have the indices of the best spatial and temporal data
                # we can get the rest of the data from these casts
                [best_hist_sal, best_hist_ptmp, best_hist_pres,
                 best_hist_lat, best_hist_long, best_hist_dates] = get_region_data(wmo_numbers,
                                                                                   float_name,
                                                                                   config,
                                                                                   index,
                                                                                   float_pres)

                best_hist_z = grid_z[index]

                # If we are near the Subantarctic Front we need to figure out if
                # the profile is north or south of it. Then we should remove data not on
                # the same side the profile is on
                if map_use_saf == 1:
                    [best_hist_sal2, best_hist_ptmp2, best_hist_pres2,
                     best_hist_lat2, best_hist_long2, best_hist_dates2, best_hist_z2] = frontal_constraint_saf(config, best_hist_sal, best_hist_ptmp,
                                                                                              best_hist_pres, best_hist_lat, best_hist_long,
                                                                                              best_hist_dates, best_hist_z, float_lat,
                                                                                              float_pres, float_tmp, float_sal)
                    # Use frontal separation only if there are at least 5 profiles
                    if len(best_hist_sal2[0]) > 5:
                        best_hist_sal = best_hist_sal2
                        best_hist_ptmp = best_hist_ptmp2
                        best_hist_pres = best_hist_pres2
                        best_hist_lat = best_hist_lat2
                        best_hist_long = best_hist_long2
                        best_hist_dates = best_hist_dates2
                        best_hist_z = best_hist_z2


                # make the float longitude wrap around the 0-360 mark if the historical data has
                if np.argwhere(best_hist_long > 360).__len__() > 0:
                    if 0 <= float_long <= 20:
                        float_long += 360

                # interpolate historical data onto float theta levels
                hist_interp_sal, hist_interp_pres = interp_climatology(best_hist_sal,
                                                                       best_hist_ptmp,
                                                                       best_hist_pres,
                                                                       float_sal,
                                                                       float_ptmp,
                                                                       float_pres)

                # map float theta levels below the exclusion point

                for n_level in range(float_level_count):

                    # for each float theta level, only use interpolated salinity
                    # that aren't NaN's
                    if not np.isnan(float_sal[n_level]) and float_pres[n_level] >= map_p_exclude:

                        max_hist_casts = np.argwhere(np.isnan(hist_interp_sal[n_level, :]) == 0)
                        hist_sal = hist_interp_sal[n_level, max_hist_casts]
                        hist_pres = hist_interp_pres[n_level, max_hist_casts]
                        hist_long = best_hist_long[max_hist_casts]
                        hist_lat = best_hist_lat[max_hist_casts]
                        hist_dates = best_hist_dates[max_hist_casts]
                        hist_z = best_hist_z[max_hist_casts]

                        # Need points +/- map_p_delta of float pressure
                        delta_index = np.argwhere(np.abs(hist_pres - float_pres[n_level])
                                                  < map_p_delta)[:, 0]
                        hist_sal = hist_sal[delta_index]
                        hist_pres = hist_pres[delta_index]
                        hist_long = hist_long[delta_index]
                        hist_lat = hist_lat[delta_index]
                        hist_dates = hist_dates[delta_index]
                        hist_z = hist_z[delta_index]

                        # only proceed with analysis if we have more than 5 points
                        if hist_sal.__len__() > 5:

                            # Need to check for statistical outliers
                            mean_sal = np.mean(hist_sal)
                            signal_sal = signal_variance(hist_sal)
                            outlier = np.argwhere(np.abs(hist_sal - mean_sal) /
                                                  np.sqrt(signal_sal) > 3)

                            # remove the statstical outliers
                            if outlier.__len__() > 0:
                                hist_sal = np.delete(hist_sal, outlier)
                                hist_pres = np.delete(hist_pres, outlier)
                                hist_long = np.delete(hist_long, outlier).reshape((-1, 1))
                                hist_lat = np.delete(hist_lat, outlier).reshape((-1, 1))
                                hist_dates = np.delete(hist_dates, outlier).reshape((-1, 1))
                                hist_z = np.delete(hist_z, outlier).reshape((-1, 1))

                            # calculate signal and noise for complete data

                            noise_sal = noise_variance(hist_sal.flatten(),
                                                       hist_lat.flatten(),
                                                       hist_long.flatten())
                            signal_sal = signal_variance(hist_sal)

                            # map residuals
                            hist_data = np.array([hist_lat, hist_long,
                                                  hist_dates, hist_z])
                            float_data = np.array([float_lat,
                                                   float_long,
                                                   float_date,
                                                   float_z]).reshape((-1, 4))

                            hist_data = np.column_stack((hist_lat, hist_long,
                                                         hist_dates, hist_z))

                            mapped_values = map_data_grid(hist_sal.flatten(),
                                                          float_data,
                                                          hist_data,
                                                          long_large, lat_large,
                                                          map_age_large,
                                                          signal_sal, noise_sal,
                                                          phi_large, map_use_pv)

                            # use short scales to map the residuals

                            sal_residual = hist_sal.flatten() - mapped_values[2]
                            sal_signal_residual = signal_variance(sal_residual)

                            mapped_residuals = map_data_grid(sal_residual,
                                                             float_data,
                                                             hist_data,
                                                             long_small, lat_small,
                                                             map_age_small,
                                                             sal_signal_residual, noise_sal,
                                                             phi_small, map_use_pv)

                            la_ptmp[n_level, profile_index] = float_ptmp[n_level]
                            la_mapped_sal[n_level, profile_index] = mapped_values[0] + \
                                                                    mapped_residuals[0]

                            dot_map_values = np.dot(mapped_values[1], mapped_values[1])
                            la_mapsalerrors[n_level, profile_index] = np.sqrt(dot_map_values +
                                                                              mapped_residuals[1] *
                                                                              mapped_residuals[1])
                            la_noise_sal[n_level, profile_index] = noise_sal
                            la_signal_sal[n_level, profile_index] = signal_sal
                            scale_long_large[profile_index] = long_large
                            scale_lat_large[profile_index] = lat_large
                            scale_long_small[profile_index] = long_small
                            scale_lat_small[profile_index] = lat_small
                            scale_phi_large[profile_index] = phi_large
                            scale_phi_small[profile_index] = phi_small
                            scale_age_large[profile_index] = map_age_large
                            scale_age_small[profile_index] = map_age_small
                            use_pv[profile_index] = map_use_pv
                            use_saf[profile_index] = map_use_saf
                            p_delta[profile_index] = map_p_delta
                            p_exclude[profile_index] = map_p_exclude

                            # only save selected historical points
                            if selected_hist.__len__() == 0:
                                selected_hist = np.array([hist_long[0][0],
                                                          hist_lat[0][0],
                                                          la_profile_no[profile_index]])
                                selected_hist = np.reshape(selected_hist, (1, 3))

                            for j in range(hist_long.__len__()):
                                m = selected_hist.shape[0]
                                b = np.array([hist_long[j][0], hist_lat[j][0]])
                                c = selected_hist[:, 0:2] - np.ones((m, 1)) * b
                                d = np.argwhere(np.abs(c[:, 0]) < 1 / 60)
                                d_1 = np.argwhere(np.abs(c[d, 1]) < 1 / 60)
                                if d_1.__len__() == 0:
                                    add_hist_data = np.array([hist_long[j][0],
                                                              hist_lat[j][0],
                                                              la_profile_no[profile_index]])
                                    selected_hist = np.vstack((selected_hist,
                                                               add_hist_data))

        print("time elapsed: ", round(time.time() - start_time, 2), " seconds")
        profile_index += 1

    # as a quality control check, just make sure salinities are between 30 and 40
    bad_sal_30 = np.argwhere(la_mapped_sal < 30)
    bad_sal_40 = np.argwhere(la_mapped_sal > 40)
    bad_sal = np.concatenate((bad_sal_30, bad_sal_40))

    for sal in bad_sal:
        la_mapped_sal[sal[0], sal[1]] = np.nan

    # sort the data by profile number
    sorted_profile_index = la_profile_no.argsort()

    # finalise the data sorted by profile number, make data type numpy arrays
    la_ptmp = la_ptmp[:, sorted_profile_index]
    la_mapped_sal = la_mapped_sal[:, sorted_profile_index]
    la_mapsalerrors = la_mapsalerrors[:, sorted_profile_index]
    la_noise_sal = la_noise_sal[:, sorted_profile_index]
    la_signal_sal = la_signal_sal[:, sorted_profile_index]

    if not isinstance(scale_long_large, np.ndarray):
        scale_long_large = np.array(scale_long_large)
        scale_lat_large = np.array(scale_lat_large)
        scale_long_small = np.array(scale_long_small)
        scale_lat_small = np.array(scale_lat_small)
        scale_phi_large = np.array(scale_phi_large)
        scale_phi_small = np.array(scale_phi_small)
        scale_age_large = np.array(scale_age_large)
        scale_age_small = np.array(scale_age_small)
        use_pv = np.array(use_pv)
        use_saf = np.array(use_saf)
        p_delta = np.array(p_delta)
        p_exclude = np.array(p_exclude)
        la_profile_no = np.array(la_profile_no)

    scale_long_large = scale_long_large[sorted_profile_index]
    scale_lat_large = scale_lat_large[sorted_profile_index]
    scale_long_small = scale_long_small[sorted_profile_index]
    scale_lat_small = scale_lat_small[sorted_profile_index]
    scale_phi_large = scale_phi_large[sorted_profile_index]
    scale_phi_small = scale_phi_small[sorted_profile_index]
    scale_age_large = scale_age_large[sorted_profile_index]
    scale_age_small = scale_age_small[sorted_profile_index]
    use_pv = use_pv[sorted_profile_index]
    use_saf = use_saf[sorted_profile_index]
    p_delta = p_delta[sorted_profile_index]
    p_exclude = p_exclude[sorted_profile_index]
    la_profile_no = la_profile_no[sorted_profile_index]

    if selected_hist.__len__() > 0:
        selected_hist_1 = np.array(sorted(selected_hist, key=lambda x: x[2]))
    selected_hist = selected_hist_1

    # define the saving location
    save_location = os.path.sep.join([config['FLOAT_MAPPED_DIRECTORY'],
                                      config['FLOAT_MAPPED_PREFIX'] + float_name + config['FLOAT_MAPPED_POSTFIX']])

    # save the data
    savemat(save_location, {'la_ptmp': la_ptmp,
                            'la_mapped_sal': la_mapped_sal,
                            'la_mapsalerrors': la_mapsalerrors,
                            'scale_long_large': scale_long_large,
                            'scale_lat_large': scale_lat_large,
                            'scale_long_small': scale_long_small,
                            'scale_lat_small': scale_lat_small,
                            'scale_phi_large': scale_phi_large,
                            'scale_phi_small': scale_phi_small,
                            'scale_age_large': scale_age_large,
                            'scale_age_small': scale_age_small,
                            'use_pv': use_pv,
                            'use_saf': use_saf,
                            'p_delta': p_delta,
                            'p_exclude': p_exclude,
                            'la_noise_sal': la_noise_sal,
                            'la_signal_sal': la_signal_sal,
                            'la_profile_no': la_profile_no,
                            'selected_hist': selected_hist})


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
