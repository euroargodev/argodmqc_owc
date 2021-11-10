"""
A set of helper method to size of update_salinity_mapping
"""
from copy import deepcopy
from pathlib import Path

import numpy as np
from scipy.interpolate import interp2d
from scipy.io import loadmat

from .data.fetchers import get_topo_grid


# pylint: disable=too-many-statements
def load_varibales_from_file(mapped_data_path, float_level_count) -> dict:
    """

    Parameters
    ----------
    mapped_data_path :
    float_level_count :

    Returns
    -------

    """
    float_mapped_data = {}
    if Path(mapped_data_path).is_file():

        # open up mapped data
        float_mapped_data = loadmat(mapped_data_path)

        if 'la_mapsalerrors' in float_mapped_data.keys():
            float_mapped_data['la_map_sal_errors'] = float_mapped_data['la_mapsalerrors']
            float_mapped_data.pop('la_mapsalerrors')

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
        float_mapped_data["profile_index"] = float_mapped_data["la_mapped_sal"].shape[1]
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
        float_mapped_data["profile_index"] = 0
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

        print(f"No precalculated data at: {mapped_data_path}")
        print("__________________________________________________________\n")

    return float_mapped_data


def get_float_data(float_source_data, missing_profile) -> dict:
    """

    Parameters
    ----------
    float_source_data :
    missing_profile :

    Returns
    -------

    """
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
    """

    Parameters
    ----------
    data :
    float_level_count :
    profile_index :

    Returns
    -------

    """
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
    """

    Parameters
    ----------
    grid_data :
    config :

    Returns
    -------

    """
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


# pylint: disable=too-many-arguments
def process_profile_hist_variables(grid_data, float_pres, hist_interp_sal, hist_interp_pres, n_level, map_p_delta):
    """

    Parameters
    ----------
    grid_data :
    float_pres :
    hist_interp_sal :
    hist_interp_pres :
    n_level :
    map_p_delta :

    Returns
    -------

    """
    max_hist_casts = np.argwhere(np.isnan(hist_interp_sal[n_level, :]) == 0)
    hist_sal = hist_interp_sal[n_level, max_hist_casts]
    hist_pres = hist_interp_pres[n_level, max_hist_casts]
    hist_long = grid_data['grid_long'][max_hist_casts]
    hist_lat = grid_data['grid_lat'][max_hist_casts]
    hist_dates = grid_data['grid_dates'][max_hist_casts]
    hist_z = grid_data['grid_z'][max_hist_casts]

    # Need points +/- map_p_delta of float pressure
    delta_index = np.argwhere(np.abs(hist_pres - float_pres[n_level]) < map_p_delta)[:, 0]
    hist_sal = hist_sal[delta_index]
    hist_pres = hist_pres[delta_index]
    hist_long = hist_long[delta_index]
    hist_lat = hist_lat[delta_index]
    hist_dates = hist_dates[delta_index]
    hist_z = hist_z[delta_index]
    return {'hist_sal': hist_sal, 'hist_pres': hist_pres, 'hist_long': hist_long,
            'hist_lat': hist_lat, 'hist_dates': hist_dates, 'hist_z': hist_z}


def remove_statical_outliers(outlier, hist_data):
    """

    Parameters
    ----------
    outlier :
    hist_data :

    Returns
    -------

    """
    if outlier.__len__() > 0:
        hist_data['hist_sal'] = np.delete(hist_data['hist_sal'], outlier)
        hist_data['hist_pres'] = np.delete(hist_data['hist_pres'], outlier)
        hist_data['hist_long'] = np.delete(hist_data['hist_long'], outlier).reshape((-1, 1))
        hist_data['hist_lat'] = np.delete(hist_data['hist_lat'], outlier).reshape((-1, 1))
        hist_data['hist_dates'] = np.delete(hist_data['hist_dates'], outlier).reshape((-1, 1))
        hist_data['hist_z'] = np.delete(hist_data['hist_z'], outlier).reshape((-1, 1))

    return hist_data


def check_and_make_numpy_arry(data):
    """

    Parameters
    ----------
    data : dictory that should contain

    Returns
    -------
    a dictionary where all elements are numpy arrays

    """
    for key, value in data.items():
        if not isinstance(value, np.ndarray):
            data[key] = np.array(value)

    return data


def sort_numpy_array(data, index, keys=None):
    """

    Parameters
    ----------
    keys : subset of elements to sort
    data : dictorain of value
    index : index to sort over

    Returns
    -------
    a dictionary where all elements are sorted within themselfes

    """
    if keys:
        data.update({key: data[key][:, index] for key in keys})
        # for key in keys:
        #    data[key] = data[key][:, index]
    else:
        # data.update({key: value[index] for (key, value) in data.items() if value.size > 1})
        for key, value in data.items():
            if value.size > 1:
                if len(value.shape) > 1:
                    data[key] = value[:, index]
                else:
                    data[key] = value[index]

    return data


def selected_historical_points(data, hist_data, profile_index):
    """

    Parameters
    ----------
    data :
    hist_data :
    profile_index :

    Returns
    -------

    """
    # only save selected historical points
    if data['selected_hist'].__len__() == 0:
        selected_hist = np.array([hist_data['hist_long'][0][0], hist_data['hist_lat'][0][0],
                                  data['la_profile_no'][profile_index]])

        selected_hist = np.reshape(selected_hist, (1, 3))
        data['selected_hist'] = selected_hist

    count = len(hist_data['hist_long'])

    for j in range(count):
        length = data['selected_hist'].shape[0]
        lon_lat = np.array([hist_data['hist_long'][j][0], hist_data['hist_lat'][j][0]])
        new_object = data['selected_hist'][:, 0:2] - np.ones((length, 1)) * lon_lat
        d_0 = np.argwhere(np.abs(new_object[:, 0]) < 1 / 60)
        d_1 = np.argwhere(np.abs(new_object[d_0, 1]) < 1 / 60)
        if len(d_1) == 0:
            add_hist_data = np.array([hist_data['hist_long'][j][0], hist_data['hist_lat'][j][0],
                                      data['la_profile_no'][profile_index]])
            if len(data['selected_hist']) == 0:
                data['selected_hist'] = add_hist_data
            else:
                data['selected_hist'] = np.vstack((data['selected_hist'], add_hist_data))

    return data
