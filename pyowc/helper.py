from copy import deepcopy
from pathlib import Path

import numpy as np
from scipy.interpolate import interp2d
from scipy.io import loadmat
from .data.fetchers import get_topo_grid

def load_varibales_from_file(mapped_data_path, float_level_count) -> dict:
    float_mapped_data = {}
    if Path(mapped_data_path).is_file():

        # open up mapped data
        float_mapped_data = loadmat(mapped_data_path)

        # flatten the data
        float_mapped_data["la_profile_no"] = float_mapped_data['la_profile_no'].flatten()
        float_mapped_data["scale_long_large"] = float_mapped_data['scale_long_large'].flatten()
        float_mapped_data["scale_lat_large"] = float_mapped_data['scale_lat_large'].flatten()
        float_mapped_data["scale_long_small"] = float_mapped_data['scale_long_small'].flatten()
        float_mapped_data["scale_lat_small"] = float_mapped_data['scale_lat_small'].flatten()
        float_mapped_data["scale_phi_large"] = float_mapped_data['scale_phi_large'].flatten()
        float_mapped_data["scale_phi_small"] = float_mapped_data['scale_phi_small'].flatten()
        float_mapped_data["scale_age_large"] = float_mapped_data['scale_age_large'].flatten()
        float_mapped_data["scale_age_small"] = float_mapped_data['scale_age_small'].flatten()
        float_mapped_data["use_pv"] = float_mapped_data['use_pv'].flatten()
        float_mapped_data["use_saf"] = float_mapped_data['use_saf'].flatten()
        float_mapped_data["p_delta"] = float_mapped_data['p_delta'].flatten()
        float_mapped_data["p_exclude"] = float_mapped_data['p_exclude'].flatten()

        # Check to see if this is an older version run without the saf constraint
        if not "use_saf" in float_mapped_data:
            float_mapped_data["use_saf"] = np.zeros(float_mapped_data['use_pv'].shape)

        # Get mapped data shape
        profile_index = float_mapped_data["la_mapped_sal"].shape[1]
        max_depth = float_mapped_data["la_mapped_sal"].shape[0]
        how_many_cols = float_mapped_data["la_mapped_sal"].shape[1]
        new_depth = float_level_count

        # if we have more data available than in the current mapped data, we need to extend
        # the matrices so we can add this data

        if new_depth > max_depth != 0:
            float_mapped_data["la_mapped_sal"] = np.insert(float_mapped_data["la_mapped_sal"],
                                                           float_mapped_data["la_mapped_sal"].shape[0],
                                                           np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                                           axis=0)
            float_mapped_data["la_map_sal_errors"] = np.insert(float_mapped_data["la_map_sal_errors"],
                                                               float_mapped_data["la_map_sal_errors"].shape[0],
                                                               np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                                               axis=0)
            float_mapped_data["la_noise_sal"] = np.insert(float_mapped_data["la_noise_sal"],
                                                          float_mapped_data["la_noise_sal"].shape[0],
                                                          np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                                          axis=0)
            float_mapped_data["la_signal_sal"] = np.insert(float_mapped_data["la_signal_sal"],
                                                           float_mapped_data["la_signal_sal"].shape[0],
                                                           np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                                           axis=0)
            float_mapped_data["la_ptmp"] = np.insert(float_mapped_data["la_ptmp"],
                                                     float_mapped_data["la_ptmp"].shape[0],
                                                     np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                                                     axis=0)

        print("Using precalculated data: ", mapped_data_path)
        print("__________________________________________________________")

    # If we don't have any precalculated mapped data

    else:
        # initialise variables
        profile_index = 0
        float_mapped_data["la_mapped_sal"] = np.empty((float_level_count, 0))
        float_mapped_data["la_map_sal_errors"] = np.empty((float_level_count, 0))
        float_mapped_data["la_noise_sal"] = np.empty((float_level_count, 0))
        float_mapped_data["la_signal_sal"] = np.empty((float_level_count, 0))
        float_mapped_data["la_ptmp"] = np.empty((float_level_count, 0))
        float_mapped_data["la_profile_no"] = np.empty(0)
        float_mapped_data["scale_long_large"] = []
        float_mapped_data["scale_lat_large"] = []
        float_mapped_data["scale_long_small"] = []
        float_mapped_data["scale_lat_small"] = []
        float_mapped_data["scale_phi_large"] = []
        float_mapped_data["scale_phi_small"] = []
        float_mapped_data["scale_age_large"] = []
        float_mapped_data["scale_age_small"] = []
        float_mapped_data["use_pv"] = []
        float_mapped_data["use_saf"] = []
        float_mapped_data["p_delta"] = []
        float_mapped_data["p_exclude"] = []
        float_mapped_data["selected_hist"] = []

        print("No precalculated data at: %s" % mapped_data_path)
        print("__________________________________________________________\n")

        return float_mapped_data


def get_float_data(float_source_data, missing_profile) -> dict:
    data = {}
    # get data from float
    data['float_lat'] = float_source_data['LAT'][0, missing_profile]
    data['float_long'] = float_source_data['LONG'][0, missing_profile]
    data['float_date'] = float_source_data['DATES'][0, missing_profile]
    data['float_sal'] = float_source_data['SAL'][:, missing_profile]
    data['float_tmp'] = float_source_data['TEMP'][:, missing_profile]
    data['float_ptmp'] = float_source_data['PTMP'][:, missing_profile]
    data['float_pres'] = float_source_data['PRES'][:, missing_profile]

    return data


def process_profiles_la_variables(data, float_level_count, profile_index):
    # if we are inserting changing a column in existing data
    if profile_index < data['la_ptmp'].shape[1]:
        data['la_ptmp'][:, profile_index] = np.nan * np.ones(float_level_count)
        data['la_mapped_sal'][:, profile_index] = np.nan * np.ones(float_level_count)
        data['la_map_sal_errors'][:, profile_index] = np.nan * np.ones(float_level_count)
        data['la_noise_sal'][:, profile_index] = np.nan * np.ones(float_level_count)
        data['la_signal_sal'][:, profile_index] = np.nan * np.ones(float_level_count)

    # if we are adding a new column
    else:
        data['la_ptmp'] = np.hstack((data['la_ptmp'],
                                     np.nan * np.ones((float_level_count, 1))))
        data['la_mapped_sal'] = np.hstack((data['la_mapped_sal'],
                                           np.nan * np.ones((float_level_count, 1))))
        data['la_map_sal_errors'] = np.hstack((data['la_map_sal_errors'],
                                               np.nan * np.ones((float_level_count, 1))))
        data['la_noise_sal'] = np.hstack((data['la_noise_sal'],
                                          np.nan * np.ones((float_level_count, 1))))
        data['la_signal_sal'] = np.hstack((data['la_signal_sal'],
                                           np.nan * np.ones((float_level_count, 1))))
    return data


def process_profiles_grid_variables(grid_data, config):
    # tbase.int file requires longitudes from 0 to +/-180
    grid_long_tbase = deepcopy(grid_data['grid_long'])

    g_180 = np.argwhere(grid_long_tbase > 180)

    grid_long_tbase[g_180] -= 360

    # find depth of the ocean at historical locations
    grid_elev, grid_x, grid_y = get_topo_grid(np.amin(grid_long_tbase) - 1,
                                              np.amax(grid_long_tbase) + 1,
                                              np.amin(grid_data['grid_lat']) - 1,
                                              np.amax(grid_data['grid_lat']) + 1,
                                              config)

    grid_interp = interp2d(grid_x[0], grid_y[:, 0],
                           grid_elev, kind='linear')

    # As a note, the reason we vectorise the function here is because we do not
    # want to compare every longitudinal value to ever latitude. Rather, we simply
    # want to interpolate each pair of longitudes and latitudes.

    grid_z = -1 * np.vectorize(grid_interp)(grid_long_tbase, grid_data['grid_lat'])

    grid_data['grid_z'] = grid_z
    grid_data['grid_x'] = grid_x
    grid_data['grid_y'] = grid_y
    grid_data['grid_elev'] = grid_elev

    return grid_data


def process_profile_best_hist():
    # If we are near the Subantarctic Front we need to figure out if
    # the profile is north or south of it. Then we should remove data not on
    # the same side the profile is on
    if map_use_saf == 1:
        [best_hist_sal2, best_hist_ptmp2, best_hist_pres2,
         best_hist_lat2, best_hist_long2, best_hist_dates2,
         best_hist_z2] = frontal_constraint_saf(config, best_hist_sal, best_hist_ptmp, best_hist_pres,
                                                best_hist_lat, best_hist_long, best_hist_dates, best_hist_z,
                                                float_lat, float_pres, float_tmp, float_sal)
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
                outlier1 = np.argwhere(np.abs(hist_sal - mean_sal) /
                                       np.sqrt(signal_sal) > 3)
                outlier = outlier1[:, 0]

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
