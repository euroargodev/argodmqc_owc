""" Functions to calibrate data

"""
import os
import time
import copy
import numpy as np
from scipy.io import loadmat, savemat
import scipy.optimize as scipy
import scipy.interpolate as interpolate

from pyowc import core
from pyowc import data
from pyowc.data.wrangling import interp_climatology


# pylint: disable=too-many-locals
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
# pylint: disable=too-many-nested-blocks
# pylint: disable=invalid-name
# pylint: disable=fixme
def update_salinity_mapping(float_dir, float_name, config):
    """ Calculates values needed for analysis. Saves them to memory to use later

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
    filename = config['FLOAT_SOURCE_DIRECTORY'] + float_dir + \
               float_name + config['FLOAT_SOURCE_POSTFIX']

    # Load the float data
    float_source_data = loadmat(filename)

    # Get the profile number and size of the data in the profile
    profile_no = float_source_data['PROFILE_NO'][0]
    float_level_count = float_source_data['SAL'].shape[0]

    # Load all the mapping parameters, including the WMO boxes -----------

    wmo_boxes = loadmat(config['CONFIG_DIRECTORY'] + config['CONFIG_WMO_BOXES'])
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
    mapped_data_path = config['FLOAT_MAPPED_DIRECTORY'] + float_dir \
                       + config['FLOAT_MAPPED_PREFIX'] + float_name + config['FLOAT_MAPPED_POSTFIX']
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
            float_elev, float_x, float_y = data.fetchers.get_topo_grid(float_long_tbase - 1,
                                                         float_long_tbase + 1,
                                                         float_lat - 1,
                                                         float_lat + 1,
                                                         config)
            float_interp = interpolate.interp2d(float_x[0, :],
                                                float_y[:, 0],
                                                float_elev,
                                                kind='linear')
            float_z = -float_interp(float_long, float_lat)[0]

            # gather data from area surrounding the float location
            wmo_numbers = core.finders.find_25boxes(float_long, float_lat, wmo_boxes)
            grid_lat, grid_long, grid_dates = data.fetchers.get_region_hist_locations(wmo_numbers,
                                                                        float_name,
                                                                        config)

            # if we have data in the surrounding area, find depths at these points
            if grid_lat.__len__() > 0:
                # tbase.int file requires longitudes from 0 to +/-180
                grid_long_tbase = copy.deepcopy(grid_long)
                g_180 = np.argwhere(grid_long_tbase > 180)
                grid_long_tbase[g_180] -= 360

                # find depth of the ocean at historical locations
                grid_elev, grid_x, grid_y = data.fetchers.get_topo_grid(np.amin(grid_long_tbase) - 1,
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

                index = core.finders.find_besthist(grid_lat, grid_long, grid_dates, grid_z,
                                      float_lat, float_long_0, float_date, float_z,
                                      lat_large, lat_small, long_large, long_small,
                                      phi_large, phi_small, map_age_large, map_age_small,
                                      map_use_pv, max_casts)

                # Now that we have the indices of the best spatial and temporal data
                # we can get the rest of the data from these casts
                [best_hist_sal, best_hist_ptmp, best_hist_pres,
                 best_hist_lat, best_hist_long, best_hist_dates] = data.fetchers.get_region_data(wmo_numbers,
                                                                                   float_name,
                                                                                   config,
                                                                                   index,
                                                                                   float_pres)

                best_hist_z = grid_z[index]

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
                            signal_sal = core.stats.signal_variance(hist_sal)
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

                            noise_sal = core.stats.noise_variance(hist_sal.flatten(),
                                                       hist_lat.flatten(),
                                                       hist_long.flatten())
                            signal_sal = core.stats.signal_variance(hist_sal)

                            # map residuals
                            hist_data = np.array([hist_lat, hist_long,
                                                  hist_dates, hist_z])
                            float_data = np.array([float_lat,
                                                   float_long,
                                                   float_date,
                                                   float_z]).reshape((-1, 4))

                            hist_data = np.column_stack((hist_lat, hist_long,
                                                         hist_dates, hist_z))

                            mapped_values = data.wrangling.map_data_grid(hist_sal.flatten(),
                                                          float_data,
                                                          hist_data,
                                                          long_large, lat_large,
                                                          map_age_large,
                                                          signal_sal, noise_sal,
                                                          phi_large, map_use_pv)

                            # use short scales to map the residuals

                            sal_residual = hist_sal.flatten() - mapped_values[2]
                            sal_signal_residual = core.stats.signal_variance(sal_residual)

                            mapped_residuals = data.wrangling.map_data_grid(sal_residual,
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
        index = selected_hist[:, 2].argsort()
        selected_hist = selected_hist[index, :]

    # define the saving location
    save_location = config['FLOAT_MAPPED_DIRECTORY'] + config['FLOAT_MAPPED_PREFIX'] + \
                    float_name + config['FLOAT_MAPPED_POSTFIX']

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


# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
# pylint: disable=too-many-arguments
def set_calseries(float_dir, float_name, system_config):
    """ Set the calseries parameters for analysis and line fitting

        Parameters
        ----------
        use_pres_gt: pressure greater than
        use_theta_gt: ptmp greater than
        use_theta_lt: ptmp less than
        use_pres_lt: pressure less than
        float_dir: location of float
        float_name: float source name
        system_config: configuration settings

        Returns
        -------
        Nothing, but save parameters
    """

    # load float source data
    float_source = loadmat(system_config['FLOAT_SOURCE_DIRECTORY'] +
                                 float_dir + float_name +
                                 system_config['FLOAT_SOURCE_POSTFIX'])

    profile_no = float_source['PROFILE_NO'].flatten()
    no_profiles = profile_no.__len__()

    # Check if we already have a calseries file

    calseries_filename = (system_config['FLOAT_CALIB_DIRECTORY'] +
                          float_dir +
                          system_config['FLOAT_CALSERIES_PREFIX'] +
                          float_name +
                          system_config['FLOAT_CALIB_POSTFIX'])

    # if we already have a calseries file, use those values. Else, use new ones
    try:
        calseries_data = loadmat(calseries_filename)
        breaks = calseries_data['breaks']
        max_breaks = calseries_data['max_breaks']
        calseries = calseries_data['calseries'].flatten()
        calib_profile_no = calseries_data['calib_profile_no'].flatten()
        use_theta_lt = calseries_data['use_theta_lt']
        use_theta_gt = calseries_data['use_theta_gt']
        use_pres_lt = calseries_data['use_pres_lt']
        use_pres_gt = calseries_data['use_pres_gt']

        # use percent may not exist, as it was added later
        try:
            use_percent_gt = calseries_data['use_percent_gt']

        except NameError:
            use_percent_gt = 0.5

        print("Using parameters found in ", calseries_filename,
              "\nTo use new parameters, delete this file")

    except FileNotFoundError:

        # Config calseries parameters

        breaks = []
        max_breaks = 4  # 0 for linear trend, -1 for offset
        calseries = np.ones((1, no_profiles)).flatten()
        # example for splitting time series at profile 33
        # calseries = np.concatenate((np.ones((1, 33)), 2 * np.ones(1,no_profiles - 33)))
        calib_profile_no = profile_no
        use_percent_gt = 0.5
        use_theta_lt = 0
        use_theta_gt = 0
        use_pres_lt = 0
        use_pres_gt = 0

    # ensure values are in a realistic range

    if use_theta_lt > 9999:
        use_theta_lt = 0

    if use_theta_gt > 9999:
        use_theta_gt = 0

    if use_pres_lt > 9999:
        use_pres_lt = 0

    if use_pres_gt > 9999:
        use_pres_gt = 0

    # Check that there are no missing profiles between source and calseries files

    missing_profiles_index = []

    for i in range(no_profiles):
        profiles = np.argwhere(calib_profile_no == profile_no[i])
        if profiles.__len__() == 0:
            missing_profiles_index.append(i)

    # Add the missing profiles to the data set
    for i in range(missing_profiles_index.__len__()):
        missing = missing_profiles_index[i]
        calib_profile_no.append(profile_no[missing])
        # set flag as the same as previous entry
        calseries = np.append(calseries, calseries(max(missing - 1, 1)))

    # sort the calseries file by profile number

    sorted_profile_no = np.argsort(calib_profile_no)
    calib_profile_no = calib_profile_no[sorted_profile_no]
    calseries = calseries[sorted_profile_no]

    # Check that we have good salinity, temperature, and pressure data

    sal = float_source['SAL']
    temp = float_source['TEMP']
    pres = float_source['PRES']

    for i in range(no_profiles):
        sal_nan = np.argwhere(~np.isnan(sal[:, i]))
        temp_nan = np.argwhere(~np.isnan(temp[:, i]))
        pres_nan = np.argwhere(~np.isnan(pres[:, i]))

        # if no good data for this profile, remove it from calseries
        if sal_nan.__len__() == 0 or \
                temp_nan.__len__() == 0 or \
                pres_nan.__len__() == 0:
            calseries[i] = 0

    savemat(calseries_filename, {'breaks': breaks,
                                       'max_breaks': max_breaks,
                                       'calseries': calseries,
                                       'calib_profile_no': calib_profile_no,
                                       'use_theta_lt': use_theta_lt,
                                       'use_theta_gt': use_theta_gt,
                                       'use_pres_lt': use_pres_lt,
                                       'use_pres_gt': use_pres_gt,
                                       'use_percent_gt': use_percent_gt})


# pylint: disable=too-many-locals
# pylint: disable=invalid-name
# pylint: disable=global-variable-undefined
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
def nlbpfun(ubrk_i):
    """ Find residual

        Parameters
        ----------
        ubrk_i: input

        Returns
        -------
        residual
    """

    global A, breaks, nbr1, ubrk_g, xf, yf, w_i, xblim

    if nbr1 > 1:
        ubrk = ubrk_g[0:nbr1 - 1]
        for i in range(nbr1, ubrk_g.__len__()):
            ubrk.append(ubrk_i)

    else:
        ubrk = ubrk_i

    m_b = ubrk.__len__()
    fnumer = np.zeros(ubrk.shape)
    fnumer[0] = np.exp(ubrk[0])

    for i in range(1, m_b):
        fnumer[i] = fnumer[i - 1] + np.exp(ubrk[i])

    fdenom = 1 + fnumer[m_b - 1]

    ftem = (xblim[1] - xblim[0]) / fdenom

    breaks = xblim[0] + ftem * fnumer

    if np.argwhere(np.diff(breaks) == 0).__len__() > 0:
        difference = np.argwhere(np.diff(breaks) == 0)
        breaks[difference + 1] = breaks[difference + 1] + 0.00001

    A, residual = core.stats.brk_pt_fit(xf, yf, w_i, breaks)

    return residual


def fit_cond(x, y, n_err, lvcov, *args):
    """ Get optimal fit

        To decide which fit is optimal, we will use the small sample variation of the
        Akaike Information Criterion.   Having chosen the
        number of parameters, we then use the F-test to see if the reduction in
        variance relative to the original variance is statistically significant.

        We have also used the correlation matrix for the horizontal and vertical
        scales to estimate the number of effective degrees of freedom for the
        fits and for estimating the uncertainties of the fits
        This function implements a non-linear fit of a piecewise linear fit.  The
        methodology is described in:
        Jones, R.H. and I. Dey, 1995, Determining one or more change points.
        Chemistry and Physics of Lipids, 76, 1-6.

        Cecile Cabanes, 2017: force the fit to an offset only if NDF <13. and display
        a warning : to track change see change config 129

        Parameters
        ----------
        x: observations
        y: observations
        n_err: error estimate for each observation
        lvcov: covariance matrix
        param: parameter for our breaks
        br: break points

        Returns
        -------
        xfit              profile number (unique of x)
        condslope         fit estimate for each profile
        condslope_err     estimated rms error for each profile
        time_deriv        estimated change per profile
        time_deriv_err    estimated rms error in change per profile
        sta_mean          mean difference between estimate and actual values
                          averaged over each profile
        sta_rms           rms difference between estimates and actual values
                          averaged over each profile
                          the off-diagonal coariance (lvcov) is taken into
                          account
    """

    # define global variables needed for line fitting
    global A, breaks, nbr1, ubrk_g, xf, yf, w_i, xblim

    # Set up some default values
    tol = 1e-06
    max_brk_dflt = 4
    max_brk_in = []
    nbr1 = -1
    brk_init = []  # guesses for break point
    setbreaks = 0

    nloops = 200  # number of loops to fit profile error

    # parameters for optimisation
    max_fun_evals = 1000

    # form x and y variables for the fit
    xfit = np.unique(x)
    nfit = xfit.__len__()

    # exclude bad points before the fit. Need to reformat matrices so that they match matlab
    x = x.T.flatten()
    y = y.T.flatten()
    n_err = n_err.T.flatten()

    good = np.argwhere(np.isfinite(y) & np.isfinite(x))

    x = x[good].flatten()
    y = y[good].flatten()
    n_err = n_err[good].flatten()

    temp_lvcov = np.empty((good.__len__(), good.__len__()))

    itx = 0
    ity = 0
    for i in good:
        for j in good:
            temp_lvcov[ity, itx] = lvcov[i, j]
            itx += 1
        itx = 0
        ity += 1

    lvcov = temp_lvcov

    npts = x.__len__()

    # check that we actually have some good data

    if npts == 0:
        print("Failed to find any good data")

        condslope = np.nan
        condslope_err = np.nan
        time_deriv = np.nan
        time_deriv_err = np.nan
        sta_mean = np.nan
        sta_rms = np.nan
        NDF = []
        fit_coef = []
        fit_breaks = []

        return (xfit, condslope, condslope_err, time_deriv, time_deriv_err, sta_mean,
                sta_rms, NDF, fit_coef, fit_breaks)

    # condition the series so that the fit is well behaved
    # sort by the independent variable

    x = np.sort(x)
    sorted_index = np.argsort(x, kind='stable')
    y = y[sorted_index]
    n_err = n_err[sorted_index]

    # scale x from -1 to 1

    x_0 = (x[npts - 1] + x[0]) / 2

    if x[0] != x[npts - 1]:
        x_scale = (x[npts - 1] - x[0]) / 2

    else:
        x_scale = 1

    # remove the mean of y and scale by the standard deviation

    y_0 = np.mean(y)
    y_scale = np.std(y)

    if y_scale == 0:
        y_scale = 1

    # calculate x and y used for fitting routine
    xf = (x - x_0) / x_scale
    yf = (y - y_0) / y_scale
    n_err = n_err / y_scale
    xfit = (xfit - x_0) / x_scale

    # get profile times that will be used as independent variables to get
    # error statstics. xp could be different to xfit if there is a profile
    # with no good data

    x_unique, index_unique = np.unique(xf, return_index=True)
    n_prof = x_unique.__len__()

    # convert errors from rms to variance
    err_var = (n_err) ** 2

    # weights for weighted least squares fit

    w_i = np.diag(np.mean(err_var) / err_var)

    # use correlation matrix to compute degrees of freedom

    ndf = np.sum(np.ones((npts, 1)) / (np.dot(lvcov, np.ones((npts, 1)))))

    # calculate the residual sum of squares for the initial series

    if x_unique.__len__() > 3:
        # find 2nd and 2nd to last profile and use them as limits for the break points
        xblim = [x_unique[1], x_unique[n_prof - 2]]

    else:
        # too few profiles
        xblim = [1, 1]

    no_args = args.__len__()

    if np.remainder(no_args, 2) != 0:
        raise ValueError("FIT_COND ERROR - inputs are incorrect")

    if no_args > 0:
        for n in range(int(no_args / 2)):
            parm = args[2 * n]
            value = args[(n * 2) + 1]

            if not isinstance(parm, str):
                raise ValueError("FIT_COND ERROR - inputs are incorrect")

            param = str.lower(parm)

            if param == 'initial_breaks':
                # initial guess for breakpoints
                brk_init = value

                # rescale
                brk_init = (brk_init - x_0) / x_scale
                brk_init = (brk_init - xblim[0]) / np.diff(xblim)

            elif param == 'max_no_breaks':
                max_brk_in = value
                nbr1 = -1

            elif param == 'number_breaks':
                pbrk = value
                nbr1 = pbrk
                max_brk_in = pbrk

            elif param == 'nloops':
                nloops = value

            elif param == 'breaks':
                if value.__len__() > 0:
                    breaks = value
                    breaks = (breaks - x_0) / x_scale
                    nbr = breaks.__len__()
                    setbreaks = 1

            else:
                raise ValueError("Paramater " + param + " not found in parameter list")

    # intialise variable for search over number of break points
    max_brk_in = int(max_brk_in)
    b_pts = np.ones((max_brk_in, max_brk_in + 1)) * np.nan
    b_A = np.ones((max_brk_in + 2, max_brk_in + 1)) * np.nan
    rss = np.ones((1, max_brk_in + 2)) * np.nan
    aic = np.ones((1, max_brk_in + 2)) * np.nan

    # check to see if we have set break points
    if setbreaks:
        if max_brk_in == 0:
            max_brk_in = nbr
            nbr1 = nbr

        # we have fixed break points
        elif max_brk_in > nbr:
            nbr1 = nbr + 1
            A, residual = core.stats.brk_pt_fit(xf, yf, w_i, breaks)
            b_pts[0:nbr, nbr + 1] = breaks.T
            b_A[0:nbr + 2, nbr + 2] = A[0:nbr + 2]
            rss[0, nbr + 2] = np.sum(residual ** 2 / err_var)
            no_param = 2 * (nbr + 1)
            aic[0, nbr + 2] = ndf * np.log(rss[0, nbr + 2] / npts) + \
                              ndf * (ndf + no_param) / (ndf - no_param - 2)

        # we have the same number of specified breaks
        else:
            nbr1 = nbr

        max_brk = max_brk_in
        pbrk = np.arange(nbr1, max_brk + 1)

    else:
        # no break points set
        if isinstance(max_brk_in, list):
            max_brk_in = max_brk_dflt

        max_brk = max_brk_in
        pbrk = np.arange(nbr1, max_brk + 1)

    if ndf < 2 * (max_brk + 2) + 1:
        if ndf > 2 * (nbr1 + 2) + 1:
            pbrk = np.arange(nbr1, np.floor((ndf - 1) / 2 - 2) + 1)
            print("WARNING: only have " + str(ndf) + " degrees of freedom")
            print("Maximum breakpoints to be tried: " + str(np.max(pbrk)))

        else:
            if setbreaks == 1:
                pbrk = np.array([nbr])
                max_brk = nbr
                nbr1 = nbr
                print("WARNING: Only have " + str(ndf) + " degrees of freedom")
                print("Estimate fit with fixed breakpoints")

            else:
                pbrk = np.array([-1])
                print("WARNING: Only have " + str(ndf) + " degrees of freedom")
                print("Estimate offset only")

    for nbr in pbrk:
        if nbr == -1:
            # offset only
            # since this is an error weighted average, yx won't necessarily be 0
            ones_column = np.ones((npts, 1))
            b_A[0, 0] = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(ones_column.T, w_i),
                                                                  ones_column)),
                                             ones_column.T),
                                      w_i),
                               yf)

            residual = yf - (ones_column * b_A[0, 0]).flatten()
            rss[0, 0] = np.sum(residual ** 2 / err_var)
            aic[0, 0] = ndf * np.log(rss[0, 0] / npts) + ndf * (ndf + 1) / (ndf - 3)

        elif nbr == 0:
            # linear fit, no break points
            A, residual = core.stats.brk_pt_fit(xf, yf, w_i)
            b_A[0:2, 1] = A[0:2]
            rss[0, 1] = np.sum(residual ** 2 / err_var)
            aic[0, 1] = ndf * np.log(rss[0, 1] / npts) + ndf * (ndf + 2) / (ndf - 4)

        else:
            nbr2 = brk_init.__len__()

            # Check if there are enough initial guesses
            if nbr2 >= nbr:
                if brk_init.shape[0] > brk_init.shape[1]:
                    brk_init = brk_init.T

                b_guess = brk_init[0:nbr]

            # first guess for breaks as evenly distributed between 2nd and 2nd to last point
            else:
                b_guess = -1 + 2 * np.arange(1, nbr + 1) / (nbr + 1)

            b_g = np.concatenate(([-1], b_guess))
            ubrk_g = []

            nbr = int(nbr)
            for n in range(nbr):
                ubrk_g.append(np.log((b_g[n + 1] - b_g[n]) / (1 - b_g[nbr])))

            if setbreaks:
                # break points are already set
                if nbr1 == max_brk:
                    A, residual = core.stats.brk_pt_fit(xf, yf, w_i, breaks)

                # fit over limited number of breaks
                else:
                    optim = scipy.least_squares(nlbpfun, ubrk_g[nbr1:nbr],
                                                method='lm', ftol=tol, max_nfev=max_fun_evals)
                    ubrk = optim['x'][0]
                    residual = optim['fun']

                    ubrk = np.concatenate((ubrk_g[0:nbr1 - 1], ubrk))
            # get non-linear least squares for break points
            else:
                ubrk_g = np.array(ubrk_g)
                optim = scipy.least_squares(nlbpfun, ubrk_g,
                                            method='lm', ftol=tol, max_nfev=max_fun_evals)
                ubrk = optim['x'][0]
                residual = optim['fun']

            b_pts[0:nbr, nbr] = breaks.T
            b_A[0:nbr + 2, nbr + 1] = A[0:nbr + 2]
            rss[0, nbr + 1] = np.sum(residual ** 2 / err_var)
            p = 2 * (nbr + 1)
            aic[0, nbr + 1] = ndf * np.log(rss[0, nbr + 1] / npts) + ndf * (ndf + p) / (ndf - p - 2)

    if setbreaks and nbr1 == max_brk:
        best = pbrk + 2

    # decide which fit to use (offset, linear, piecewise)
    else:
        if nbr1 > 1:
            pbrk = np.arange((nbr1 - 1), max_brk)

        good = np.array(pbrk + 1, dtype=int)
        best = np.argmin(aic[0, good])

        if isinstance(good, np.ndarray):
            best = good[best] + 1
        else:
            best = good + 1

    if setbreaks & nbr1 == max_brk:
        comment = "Fit evaluated "

    else:
        comment = "Best model found with "

    if best > 2:
        comment = comment + str(best - 2) + " break points"

    elif best == 2:
        comment = comment + "linear fit"

    else:
        comment = comment + "offset value only"

    print(comment)

    if best > 2:
        breaks = b_pts[np.arange(0, best - 2), best - 2].T

    else:
        breaks = []

    best = int(best)
    A = b_A[0:best, best - 1]
    btem = np.concatenate(([xf[0]], breaks))
    E = np.zeros((npts, best))
    E[:, 0] = np.ones(npts).T
    ixb = core.utils.sorter(btem, xf)

    if best > 1:
        for j in range(best - 1):
            ib = np.argwhere(ixb == j)
            E[ib, j + 1] = xf[ib] - btem[j]
            ii = np.argwhere(ixb > j)

            if ii.__len__() > 0:
                E[ii, j + 1] = btem[j + 1] - btem[j]

    # get uncertainnties in fit parameters
    B_i = (np.dot(np.dot(E.T, w_i), E)) ** -1
    P = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(B_i, E.T),
                                           w_i), np.diag(err_var)),
                             w_i),
                      E),
               B_i)
    P_1 = np.diag(P)
    P = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(B_i, E.T),
                                                  w_i),
                                           np.diag(err_var)),
                                    lvcov),
                             w_i),
                      E),
               B_i)
    P_2 = np.diag(P)

    # reduce matrix to have only one value per profile
    btem = np.concatenate(([xfit[0]], breaks))
    E = np.zeros((xfit.__len__(), best))
    E[:, 0] = np.ones((xfit.__len__(), 1)).T
    ixb = core.utils.sorter(btem, xfit)

    if best >= 2:
        for j in range(best - 1):
            # point to x values greater than break j
            ib = np.argwhere(ixb == j)
            E[ib, j + 1] = xfit[ib] - btem[j]
            # point to x values less than the one just below x
            ii = np.argwhere(ixb > j)

            if ii.__len__() > 0:
                E[ii, j + 1] = btem[j + 1] - btem[j]

    # fit the values
    yfit = np.dot(E, A)

    # factor to increase monte carlo error estimate, taking into account
    # that we have correlated noise
    P_3 = np.dot(E, P_2) / np.dot(E, P_1)

    real_A = copy.deepcopy(A)
    real_E = copy.deepcopy(E)
    real_breaks = copy.deepcopy(breaks)
    real_xf = copy.deepcopy(xf)
    real_yf = copy.deepcopy(yf)
    ubrk_g = []

    if best == 1:
        err = np.ones((nfit, 1)) * P_2[0]

    else:
        err = 0

        for i in range(nloops):

            yf = real_yf + n_err * np.random.randn(yf.size)

            if best == 2:
                # E for linear case is already calculated
                A, residual = core.stats.brk_pt_fit(xf, yf, w_i)

            elif setbreaks:
                # E stays fixed if breaks are specified
                A, residual = core.stats.brk_pt_fit(xf, yf, w_i, breaks)

            else:
                # give an initial guess as the fitted break points to speed up calculation
                nbr = real_breaks.__len__()
                b_g = np.concatenate(([-1], real_breaks))

                for n in range(nbr):
                    ubrk_g.append(np.log((b_g[n + 1] - b_g[n]) / (1 - b_g[nbr])))

                optim = scipy.least_squares(nlbpfun, ubrk_g,
                                            method='lm', ftol=tol, max_nfev=max_fun_evals)

                ubrk = optim['x'][0]
                residual = optim['fun']

                btem = np.concatenate([xfit[0]], breaks)
                E = np.zeros((xfit.__len__(), best))
                E[:, 0] = np.ones((xfit.__len__(), 1)).T
                ixb = core.utils.sorter(btem, xfit)

                for j in range(best - 1):
                    # pointer to x values greater than break point j
                    ib = np.argwhere(ixb == j)
                    E[ib, j + 1] = xfit[ib] - btem[j]
                    # pointer to break points less than the one just below x
                    ii = np.argwhere(ixb > j)

                    if ii.__len__() > 0:
                        E[ii, j + 1] = btem[j + 1] - btem[j]

            err = err + (yfit - np.dot(E, A) ** 2)

        err = err / nloops

        # rescale error to reflect the decrease to the off diagonal covariances

        err = err * P_3

    A = copy.deepcopy(real_A)
    E = copy.deepcopy(real_E)
    breaks = copy.deepcopy(real_breaks)
    xf = copy.deepcopy(real_xf)
    yf = copy.deepcopy(real_yf)

    # get residual statistics for each profile

    w = np.diag(w_i)
    sta_mean = np.ones((1, nfit)) * np.nan
    sta_rms = np.ones((1, nfit)) * np.nan
    ip_1 = np.concatenate(([-1], index_unique.flatten()))

    for n in range(nfit):
        index = np.argwhere(x_unique == xfit[n])

        if index.__len__() > 0:
            w_values = w[np.arange(ip_1[index][0, 0] + 1, ip_1[index + 1][0, 0] + 1)]
            yf_values = yf[np.arange(ip_1[index][0, 0] + 1, ip_1[index + 1][0, 0] + 1)]

            sta_mean[0, n] = np.sum(w_values * yf_values) / np.sum(w_values)

            sta_rms[0, n] = np.sum(w_values * (sta_mean[0, n] - yf_values) ** 2) / np.sum(w_values)
            sta_rms[0, n] = np.sqrt(sta_rms[0, n])

    # convert back to original units

    x_unique = x_unique * x_scale + x_0
    xfit = xfit * x_scale + x_0

    # pass the coeffecients and break points back
    break_pts = []
    if breaks.__len__() > 0:
        break_pts = breaks * x_scale + x_0

    A[0] = A[0] * y_scale * y_0
    P[0] = P[0] * y_scale

    if A.__len__() > 1:
        A[1:best] = A[1:best] * y_scale / x_scale
        P[1:best] = P[1:best] * y_scale / x_scale

    yfit = yfit * y_scale + y_0
    err = np.sqrt(err) * y_scale
    yerr = err
    n_err = n_err * y_scale
    sta_mean = sta_mean * y_scale + y_0
    sta_rms = sta_rms * y_scale

    # get time derivatives and time derivatives errors
    ixb = core.utils.sorter(np.concatenate(([np.min(xfit)], break_pts)), xfit)

    if best == 1:
        time_deriv = np.zeros((nfit, 1))
        time_deriv_err = np.ones((nfit, 1)) * np.nan

    elif best == 2:
        time_deriv = A[1] * np.ones((nfit, 1))
        time_deriv_err = P_2[1] * np.ones((nfit, 1))

    else:
        time_deriv = np.ones((nfit, 1)) * np.nan
        time_deriv_err = np.ones((nfit, 1)) * np.nan
        for j in range(best - 1):
            ib = np.argwhere(ixb == j)
            time_deriv[ib] = A[j + 1]
            time_deriv_err[ib] = P_2[j + 1]

    condslope = yfit.T
    condslope_err = yerr.T

    # return fit parameters
    fit_coef = A
    fit_breaks = break_pts

    return (xfit, condslope, condslope_err, time_deriv, time_deriv_err,
            sta_mean, sta_rms, ndf, fit_coef, fit_breaks)
