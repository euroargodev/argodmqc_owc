"""
-----Update Salinity Mapping-----

Written by: Breck Owens
When: xx/11/2007
Converted to python by: Edward Small
When: 02/02/2020

Annie Wong, xx June. 2010: Unexplained changes

Cecile Cabanes, June. 2013: use of "map_large_scale" (time scale) used to map the large scale field

Function used to store all the resulting variables needed to run the analysis for each float profile.

Also saves all the settings used for the analysis, in case an operator wishes to either run the analysis
again, or share the analysis with another operator.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import struct
import time
import numpy as np
import scipy.io as scipy
import scipy.interpolate as interpolate
from ow_calibration.find_besthist.find_besthist import find_besthist
from ow_calibration.find_25boxes.find_25boxes import find_25boxes
from ow_calibration.interp_climatology.interp_climatology import interp_climatology
from ow_calibration.get_region.get_region_hist_locations import get_region_hist_locations
from ow_calibration.get_region.get_region_data import get_region_data
from ow_calibration.map_data_grid.map_data_grid import map_data_grid
from ow_calibration.noise_variance.noise_variance import noise_variance
from ow_calibration.signal_variance.signal_variance import signal_variance
from ow_calibration.tbase_decoder.tbase_decoder import get_topo_grid


def update_salinity_mapping(float_dir, float_name, config):
    """
    Calculates values needed for analysis. Saves them to memory to use later

    :param float_dir: directory where the float being analysed is stored
    :param float_name: name of the float being analysed
    :param config: configuration settings set by the user
    :return: Nothing, but does save values
    """

    # Load float source data ---------------------------------------------

    # Get float file name
    filename = config['FLOAT_SOURCE_DIRECTORY'] + float_dir + float_name + config['FLOAT_SOURCE_POSTFIX']

    # Load the float data
    float_source_data = scipy.loadmat(filename)

    # Get the profile number and size of the data in the profile
    profile_no = float_source_data['PROFILE_NO'][0]
    float_level_count = float_source_data['SAL'].shape[0]
    float_profile_count = float_source_data['SAL'].shape[1]

    # Load all the mapping parameters, including the WMO boxes -----------

    wmo_boxes = scipy.loadmat(config['CONFIG_DIRECTORY'] + config['CONFIG_WMO_BOXES'])
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
    try:

        # open up mapped data
        float_mapped_data = scipy.loadmat(config['FLOAT_MAPPED_DIRECTORY'] + float_dir +
                                          config['FLOAT_MAPPED_PREFIX'] + float_name +
                                          config['FLOAT_MAPPED_POSTFIX'])

        # Save the data to variables
        la_mapped_sal = float_mapped_data['la_mapped_sal']
        la_mapsalerrors = float_mapped_data['la_mapsalerrors']
        la_noise_sal = float_mapped_data['la_noise_sal']
        la_signal_sal = float_mapped_data['la_signal_sal']
        la_ptmp = float_mapped_data['la_ptmp']
        la_profile_no = float_mapped_data['la_profile_no']

        # Check to see if this is an older version run without the saf constraint
        if "use_saf" in float_mapped_data:
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

        print("Using precaulcated data: ", config['FLOAT_MAPPED_DIRECTORY'] + float_dir +
              config['FLOAT_MAPPED_PREFIX'] + float_name + config['FLOAT_MAPPED_POSTFIX'])
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

        print("No precaulcated data")
        print("__________________________________________________________\n")

    # Compare profile numbers in the float source against the mapped data matrix -
    missing_profile_index = []

    for i in range(0, profile_no.__len__()):
        profiles = np.argwhere(la_profile_no == profile_no[i])
        if profiles.size == 0:
            missing_profile_index.append(i)

    # initalise vectors for holding parameters
    scale_long_large = []
    scale_lat_large = []
    scale_long__small = []
    scale_lat_small = []
    scale_phi_large = []
    scale_phi_small = []
    scale_age_large = []
    scale_age_small = []
    use_pav = []
    use_saf = []
    p_delta = []
    p_exclude = []

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
            la_ptmp = np.hstack((la_ptmp, np.nan * np.ones((float_level_count, 1))))
            la_mapped_sal = np.hstack((la_mapped_sal, np.nan * np.ones((float_level_count, 1))))
            la_mapsalerrors = np.hstack((la_mapsalerrors, np.nan * np.ones((float_level_count, 1))))
            la_noise_sal = np.hstack((la_noise_sal, np.nan * np.ones((float_level_count, 1))))
            la_signal_sal = np.hstack((la_signal_sal, np.nan * np.ones((float_level_count, 1))))

        # initialise matrices to hold and save parameter settings
        scale_long_large.append(np.nan)
        scale_lat_large.append(np.nan)
        scale_long__small.append(np.nan)
        scale_lat_small.append(np.nan)
        scale_phi_large.append(np.nan)
        scale_phi_small.append(np.nan)
        scale_age_large.append(np.nan)
        scale_age_small.append(np.nan)
        use_pav.append(np.nan)
        use_saf.append(np.nan)
        p_delta.append(np.nan)
        p_exclude.append(np.nan)

        # get data from float
        float_lat = float_source_data['LAT'][0, missing_profile]
        float_long = float_source_data['LONG'][0, missing_profile]
        float_date = float_source_data['DATES'][0, missing_profile]
        float_sal = float_source_data['SAL'][:, missing_profile]
        float_ptmp = float_source_data['PTMP'][:, missing_profile]
        float_pres = float_source_data['PRES'][:, missing_profile]

        # if we have any good location data and pressure data from the float
        if not np.isnan(float_long) \
                and not np.isnan(float_lat) \
                and np.argwhere(np.isnan(float_pres) == 0).any():

            # tbase.int file requires longitudes from 0 to +/-180
            float_long_tbase = float_long
            if float_long_tbase > 180:
                float_long_tbase -= 360

            # find the depth of the ocean at the float location
            float_elev, float_x, float_y = get_topo_grid(float_long_tbase - 1, float_long_tbase + 1,
                                                         float_lat - 1, float_lat + 1)
            float_interp = interpolate.interp2d(float_x[0, :], float_y[:, 0], float_elev, kind='linear')
            float_z = -float_interp(float_long, float_lat)[0]

            # gather data from area surrounding the float location
            wmo_numbers = find_25boxes(float_long, float_lat, wmo_boxes)
            grid_lat, grid_long, grid_dates = get_region_hist_locations(wmo_numbers,
                                                                        float_name,
                                                                        config)

            # if we have data in the surrouding area, find depths at these points
            if grid_lat.__len__() > 0:
                # tbase.int file requires longitudes from 0 to +/-180
                grid_long_tbase = grid_long
                g_180 = np.argwhere(grid_long_tbase > 180)
                grid_long_tbase[g_180] -= 360

                # find depth of the ocean at historical locations
                grid_elev, grid_x, grid_y = get_topo_grid(np.amin(grid_long_tbase) - 1,
                                                          np.amax(grid_long_tbase) + 1,
                                                          np.amin(grid_lat) - 1,
                                                          np.amax(grid_lat) + 1)

                grid_interp = interpolate.interp2d(grid_x[0], grid_y[:, 0],
                                                   grid_elev, kind='linear')
                """
                As a note, the reason we vectorise the function here is because we do not 
                want to compare every longitudinal value to ever latitude. Rather, we simply 
                want to interpolate each pair of longitudes and latitudes.
                """
                grid_z = -np.vectorize(grid_interp)(grid_long_tbase, grid_lat)

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

                # Now that we have the best data historical data for this profile we can reset
                # the grid matrices

                grid_lat = []
                grid_long = []
                grid_dates = []

                # If we are near the Subantarctic Front we need to figure out if
                # the profile is north or south of it. Then we should remove data not on
                # the same side the profile is on
                if map_use_saf == 1:
                    # TODO: Add the SAF functions here
                    print("SAF functions unavailable in the current version")

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

                # now that we have the interpolated data, we can reset the historical data matrices
                best_hist_sal = []
                best_hist_ptmp = []
                best_hist_pres = []
                wmo_numbers = []
                index = []

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

                            mapped_values = map_data_grid(hist_sal.flatten(),
                                                          np.array([float_lat, float_long,
                                                                    float_date, float_z]).reshape((-1, 4)),
                                                          hist_data.reshape((-1, 4)),
                                                          long_large, lat_large,
                                                          map_age_large,
                                                          signal_sal, noise_sal,
                                                          phi_large, map_use_pv)

                            # use short scales to map the residuals

                            sal_residual = hist_sal.flatten() - mapped_values[2]
                            sal_signal_residual = signal_variance(sal_residual)

                            mapped_residuals = map_data_grid(sal_residual,
                                                             np.array([float_lat, float_long,
                                                                       float_date, float_z]).reshape((-1, 4)),
                                                             hist_data.reshape((-1, 4)),
                                                             long_small, lat_small,
                                                             map_age_small,
                                                             sal_signal_residual, noise_sal,
                                                             phi_small, map_use_pv)

                            la_ptmp[n_level, profile_index] = float_ptmp[n_level]
                            la_mapped_sal[n_level, profile_index] = mapped_values[0] + mapped_residuals[0]
                            la_mapsalerrors[n_level, profile_index] = np.sqrt(np.dot(mapped_values[1],
                                                                                     mapped_values[1]) +
                                                                              mapped_residuals[1] * mapped_residuals[1])
                            la_noise_sal[n_level, profile_index] = noise_sal
                            la_signal_sal[n_level, profile_index] = signal_sal
                            scale_long_large[profile_index] = long_large
                            scale_lat_large[profile_index] = lat_large
                            scale_long__small[profile_index] = long_small
                            scale_lat_small[profile_index] = lat_small
                            scale_phi_large[profile_index] = phi_large
                            scale_phi_small[profile_index] = phi_small
                            scale_age_large[profile_index] = map_age_large
                            scale_age_small[profile_index] = map_age_small
                            use_pav[profile_index] = map_use_pv
                            use_saf[profile_index] = map_use_saf
                            p_delta[profile_index] = map_p_delta
                            p_exclude[profile_index] = map_p_exclude

                            # only save selected historical points
                            if selected_hist.__len__() == 0:
                                selected_hist = np.array(
                                    [hist_long[0][0], hist_lat[0][0], la_profile_no[profile_index]])
                                selected_hist = np.reshape(selected_hist, (1, 3))

                            for j in range(hist_long.__len__()):
                                m, n = selected_hist.shape
                                b = np.array([hist_long[j][0], hist_lat[j][0]])
                                c = selected_hist[:, 0:2] - np.ones((m, 1)) * b
                                d = np.argwhere(np.abs(c[:, 0]) < 1 / 60)
                                d_1 = np.argwhere(np.abs(c[d, 1]) < 1 / 60)
                                if d_1.__len__() == 0:
                                    selected_hist = np.vstack((selected_hist,
                                                               np.array([hist_long[j][0], hist_lat[j][0],
                                                                         la_profile_no[profile_index]])))

        print("time elapsed: ", round(time.time() - start_time, 2), " seconds")
        profile_index += 1

    # as a quality control check, just make sure salinities are between 30 and 40
    bad_sal_30 = np.argwhere(la_mapped_sal < 30)
    bad_sal_40 = np.argwhere(la_mapped_sal > 40)
    bad_sal = np.concatenate((bad_sal_30, bad_sal_40))

    for sal in bad_sal:
        la_mapped_sal[sal[0], sal[1]] = np.nan

    # sort the data by profile number
    sorted_profile_no = np.sort(la_profile_no)
    sorted_profile_index = la_profile_no.argsort()

    la_ptmp = la_ptmp[:, sorted_profile_index]
