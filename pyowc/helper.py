"""A set of helper method to size of update_salinity_mapping"""

import os
from copy import deepcopy
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat, savemat

from .data.fetchers import get_topo_grid


# pylint: disable=too-many-statements
def load_varibales_from_file(mapped_data_path, float_level_count) -> dict:
    """Parameters
    ----------
    mapped_data_path :
    float_level_count :

    Returns:
    -------

    """
    float_mapped_data = {}
    if Path(mapped_data_path).is_file():
        # open up mapped data
        float_mapped_data = loadmat(mapped_data_path)

        if "la_mapsalerrors" in float_mapped_data:
            float_mapped_data["la_map_sal_errors"] = float_mapped_data["la_mapsalerrors"]
            float_mapped_data.pop("la_mapsalerrors")

        # flatten the data
        float_mapped_data["la_profile_no"] = float_mapped_data["la_profile_no"].flatten()
        float_mapped_data["scale_long_large"] = float_mapped_data["scale_long_large"].flatten()
        float_mapped_data["scale_lat_large"] = float_mapped_data["scale_lat_large"].flatten()
        float_mapped_data["scale_long_small"] = float_mapped_data["scale_long_small"].flatten()
        float_mapped_data["scale_lat_small"] = float_mapped_data["scale_lat_small"].flatten()
        float_mapped_data["scale_phi_large"] = float_mapped_data["scale_phi_large"].flatten()
        float_mapped_data["scale_phi_small"] = float_mapped_data["scale_phi_small"].flatten()
        float_mapped_data["scale_age_large"] = float_mapped_data["scale_age_large"].flatten()
        float_mapped_data["scale_age_small"] = float_mapped_data["scale_age_small"].flatten()
        float_mapped_data["use_pv"] = float_mapped_data["use_pv"].flatten()
        float_mapped_data["use_saf"] = float_mapped_data["use_saf"].flatten()
        float_mapped_data["p_delta"] = float_mapped_data["p_delta"].flatten()
        float_mapped_data["p_exclude"] = float_mapped_data["p_exclude"].flatten()

        # Check to see if this is an older version run without the saf constraint
        if "use_saf" not in float_mapped_data:
            float_mapped_data["use_saf"] = np.zeros(float_mapped_data["use_pv"].shape)

        # Get mapped data shape
        float_mapped_data["profile_index"] = float_mapped_data["la_mapped_sal"].shape[1]
        max_depth = float_mapped_data["la_mapped_sal"].shape[0]
        how_many_cols = float_mapped_data["la_mapped_sal"].shape[1]
        new_depth = float_level_count

        # if we have more data available than in the current mapped data, we need to extend
        # the matrices so we can add this data

        if new_depth > max_depth != 0:
            float_mapped_data["la_mapped_sal"] = np.insert(
                float_mapped_data["la_mapped_sal"],
                float_mapped_data["la_mapped_sal"].shape[0],
                np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                axis=0,
            )
            float_mapped_data["la_map_sal_errors"] = np.insert(
                float_mapped_data["la_map_sal_errors"],
                float_mapped_data["la_map_sal_errors"].shape[0],
                np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                axis=0,
            )
            float_mapped_data["la_noise_sal"] = np.insert(
                float_mapped_data["la_noise_sal"],
                float_mapped_data["la_noise_sal"].shape[0],
                np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                axis=0,
            )
            float_mapped_data["la_signal_sal"] = np.insert(
                float_mapped_data["la_signal_sal"],
                float_mapped_data["la_signal_sal"].shape[0],
                np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                axis=0,
            )
            float_mapped_data["la_ptmp"] = np.insert(
                float_mapped_data["la_ptmp"],
                float_mapped_data["la_ptmp"].shape[0],
                np.ones((new_depth - max_depth, how_many_cols)) * np.nan,
                axis=0,
            )

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
    """Parameters
    ----------
    float_source_data :
    missing_profile :

    Returns:
    -------

    """
    data = {}
    # get data from float
    data["float_lat"] = float_source_data["LAT"][0, missing_profile]
    data["float_long"] = float_source_data["LONG"][0, missing_profile]
    data["float_date"] = float_source_data["DATES"][0, missing_profile]
    data["float_sal"] = float_source_data["SAL"][:, missing_profile]
    data["float_tmp"] = float_source_data["TEMP"][:, missing_profile]
    data["float_ptmp"] = float_source_data["PTMP"][:, missing_profile]
    data["float_pres"] = float_source_data["PRES"][:, missing_profile]

    return data


def process_profiles_la_variables(data, float_level_count, profile_index):
    """Parameters
    ----------
    data :
    float_level_count :
    profile_index :

    Returns:
    -------

    """
    # if we are inserting changing a column in existing data
    if profile_index < data["la_ptmp"].shape[1]:
        data["la_ptmp"][:, profile_index] = np.nan * np.ones(float_level_count)
        data["la_mapped_sal"][:, profile_index] = np.nan * np.ones(float_level_count)
        data["la_map_sal_errors"][:, profile_index] = np.nan * np.ones(float_level_count)
        data["la_noise_sal"][:, profile_index] = np.nan * np.ones(float_level_count)
        data["la_signal_sal"][:, profile_index] = np.nan * np.ones(float_level_count)

    # if we are adding a new column
    else:
        data["la_ptmp"] = np.hstack((data["la_ptmp"], np.nan * np.ones((float_level_count, 1))))
        data["la_mapped_sal"] = np.hstack((data["la_mapped_sal"], np.nan * np.ones((float_level_count, 1))))
        data["la_map_sal_errors"] = np.hstack((data["la_map_sal_errors"], np.nan * np.ones((float_level_count, 1))))
        data["la_noise_sal"] = np.hstack((data["la_noise_sal"], np.nan * np.ones((float_level_count, 1))))
        data["la_signal_sal"] = np.hstack((data["la_signal_sal"], np.nan * np.ones((float_level_count, 1))))
    return data


def process_profiles_grid_variables(grid_data, config):
    """Parameters
    ----------
    grid_data :
    config :

    Returns:
    -------

    """
    # tbase.int file requires longitudes from 0 to +/-180
    grid_long_tbase = deepcopy(grid_data["grid_long"])

    g_180 = np.argwhere(grid_long_tbase > 180)

    grid_long_tbase[g_180] -= 360

    # find depth of the ocean at historical locations
    grid_elev, grid_x, grid_y = get_topo_grid(
        np.amin(grid_long_tbase) - 1,
        np.amax(grid_long_tbase) + 1,
        np.amin(grid_data["grid_lat"]) - 1,
        np.amax(grid_data["grid_lat"]) + 1,
        config,
    )

    grid_interp = RegularGridInterpolator((grid_y[:, 0], grid_x[0]), grid_elev, method="linear")

    points = np.column_stack((grid_data["grid_lat"].ravel(), grid_long_tbase.ravel()))

    grid_z = -1 * grid_interp(points).reshape(grid_data["grid_lat"].shape)

    grid_data["grid_z"] = grid_z
    grid_data["grid_x"] = grid_x
    grid_data["grid_y"] = grid_y
    grid_data["grid_elev"] = grid_elev

    return grid_data


# pylint: disable=too-many-arguments
def process_profile_hist_variables(grid_data, float_pres, hist_interp_sal, hist_interp_pres, n_level, map_p_delta):
    """Parameters
    ----------
    grid_data :
    float_pres :
    hist_interp_sal :
    hist_interp_pres :
    n_level :
    map_p_delta :

    Returns:
    -------

    """
    max_hist_casts = np.argwhere(np.isnan(hist_interp_sal[n_level, :]) == 0)
    hist_sal = hist_interp_sal[n_level, max_hist_casts]
    hist_pres = hist_interp_pres[n_level, max_hist_casts]
    hist_long = grid_data["grid_long"][max_hist_casts]
    hist_lat = grid_data["grid_lat"][max_hist_casts]
    hist_dates = grid_data["grid_dates"][max_hist_casts]
    hist_z = grid_data["grid_z"][max_hist_casts]

    # Need points +/- map_p_delta of float pressure
    delta_index = np.argwhere(np.abs(hist_pres - float_pres[n_level]) < map_p_delta)[:, 0]
    hist_sal = hist_sal[delta_index]
    hist_pres = hist_pres[delta_index]
    hist_long = hist_long[delta_index]
    hist_lat = hist_lat[delta_index]
    hist_dates = hist_dates[delta_index]
    hist_z = hist_z[delta_index]
    return {
        "hist_sal": hist_sal,
        "hist_pres": hist_pres,
        "hist_long": hist_long,
        "hist_lat": hist_lat,
        "hist_dates": hist_dates,
        "hist_z": hist_z,
    }


def remove_statical_outliers(outlier, hist_data):
    """Parameters
    ----------
    outlier :
    hist_data :

    Returns:
    -------

    """
    if outlier.__len__() > 0:
        hist_data["hist_sal"] = np.delete(hist_data["hist_sal"], outlier)
        hist_data["hist_pres"] = np.delete(hist_data["hist_pres"], outlier)
        hist_data["hist_long"] = np.delete(hist_data["hist_long"], outlier).reshape((-1, 1))
        hist_data["hist_lat"] = np.delete(hist_data["hist_lat"], outlier).reshape((-1, 1))
        hist_data["hist_dates"] = np.delete(hist_data["hist_dates"], outlier).reshape((-1, 1))
        hist_data["hist_z"] = np.delete(hist_data["hist_z"], outlier).reshape((-1, 1))

    return hist_data


def check_and_make_numpy_arry(data):
    """Parameters
    ----------
    data : dictory that should contain

    Returns:
    -------
    a dictionary where all elements are numpy arrays

    """
    for key, value in data.items():
        if not isinstance(value, np.ndarray):
            data[key] = np.array(value)

    return data


def sort_numpy_array(data, index, keys=None):
    """Parameters
    ----------
    keys : subset of elements to sort
    data : dictorain of value
    index : index to sort over

    Returns:
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
    """Parameters
    ----------
    data :
    hist_data :
    profile_index :

    Returns:
    -------

    """
    # only save selected historical points
    if data["selected_hist"].__len__() == 0:
        selected_hist = np.array(
            [hist_data["hist_long"][0][0], hist_data["hist_lat"][0][0], data["la_profile_no"][profile_index]]
        )

        selected_hist = np.reshape(selected_hist, (1, 3))
        data["selected_hist"] = selected_hist

    count = len(hist_data["hist_long"])

    for j in range(count):
        length = data["selected_hist"].shape[0]
        lon_lat = np.array([hist_data["hist_long"][j][0], hist_data["hist_lat"][j][0]])
        new_object = data["selected_hist"][:, 0:2] - np.ones((length, 1)) * lon_lat
        d_0 = np.argwhere(np.abs(new_object[:, 0]) < 1 / 60)
        d_1 = np.argwhere(np.abs(new_object[d_0, 1]) < 1 / 60)
        if len(d_1) == 0:
            add_hist_data = np.array(
                [hist_data["hist_long"][j][0], hist_data["hist_lat"][j][0], data["la_profile_no"][profile_index]]
            )
            if len(data["selected_hist"]) == 0:
                data["selected_hist"] = add_hist_data
            else:
                data["selected_hist"] = np.vstack((data["selected_hist"], add_hist_data))

    return data


def create_la_wmo_boxes_file(config: dict):
    """Converted from MATLAB version by Delphine Dobler, August 2024."""
    clim_ctd_= os.path.sep.join([config["HISTORICAL_DIRECTORY"], config["HISTORICAL_CTD_PREFIX"]])
    clim_argo_= os.path.sep.join([config["HISTORICAL_DIRECTORY"], config["HISTORICAL_ARGO_PREFIX"]])
    clim_bottle_= os.path.sep.join([config["HISTORICAL_DIRECTORY"], config["HISTORICAL_BOTTLE_PREFIX"]])
    la_wmo_boxes_file= os.path.sep.join([config["CONFIG_DIRECTORY"], config["CONFIG_WMO_BOXES"]])

    la_wmo_boxes = np.zeros((648, 4))
    la_wmo_boxes[:,1] = [
        1800,1700,1600,1500,1400,1300,1200,1100,1000,3000,3100,
        3200,3300,3400,3500,3600,3700,3800,1801,1701,1601,1501,
        1401,1301,1201,1101,1001,3001,3101,3201,3301,3401,3501,
        3601,3701,3801,1802,1702,1602,1502,1402,1302,1202,1102,
        1002,3002,3102,3202,3302,3402,3502,3602,3702,3802,1803,
        1703,1603,1503,1403,1303,1203,1103,1003,3003,3103,3203,
        3303,3403,3503,3603,3703,3803,1804,1704,1604,1504,1404,
        1304,1204,1104,1004,3004,3104,3204,3304,3404,3504,3604,
        3704,3804,1805,1705,1605,1505,1405,1305,1205,1105,1005,
        3005,3105,3205,3305,3405,3505,3605,3705,3805,1806,1706,
        1606,1506,1406,1306,1206,1106,1006,3006,3106,3206,3306,
        3406,3506,3606,3706,3806,1807,1707,1607,1507,1407,1307,
        1207,1107,1007,3007,3107,3207,3307,3407,3507,3607,3707,
        3807,1808,1708,1608,1508,1408,1308,1208,1108,1008,3008,
        3108,3208,3308,3408,3508,3608,3708,3808,1809,1709,1609,
        1509,1409,1309,1209,1109,1009,3009,3109,3209,3309,3409,
        3509,3609,3709,3809,1810,1710,1610,1510,1410,1310,1210,
        1110,1010,3010,3110,3210,3310,3410,3510,3610,3710,3810,
        1811,1711,1611,1511,1411,1311,1211,1111,1011,3011,3111,
        3211,3311,3411,3511,3611,3711,3811,1812,1712,1612,1512,
        1412,1312,1212,1112,1012,3012,3112,3212,3312,3412,3512,
        3612,3712,3812,1813,1713,1613,1513,1413,1313,1213,1113,
        1013,3013,3113,3213,3313,3413,3513,3613,3713,3813,1814,
        1714,1614,1514,1414,1314,1214,1114,1014,3014,3114,3214,
        3314,3414,3514,3614,3714,3814,1815,1715,1615,1515,1415,
        1315,1215,1115,1015,3015,3115,3215,3315,3415,3515,3615,
        3715,3815,1816,1716,1616,1516,1416,1316,1216,1116,1016,
        3016,3116,3216,3316,3416,3516,3616,3716,3816,1817,1717,
        1617,1517,1417,1317,1217,1117,1017,3017,3117,3217,3317,
        3417,3517,3617,3717,3817,7817,7717,7617,7517,7417,7317,
        7217,7117,7017,5017,5117,5217,5317,5417,5517,5617,5717,
        5817,7816,7716,7616,7516,7416,7316,7216,7116,7016,5016,
        5116,5216,5316,5416,5516,5616,5716,5816,7815,7715,7615,
        7515,7415,7315,7215,7115,7015,5015,5115,5215,5315,5415,
        5515,5615,5715,5815,7814,7714,7614,7514,7414,7314,7214,
        7114,7014,5014,5114,5214,5314,5414,5514,5614,5714,5814,
        7813,7713,7613,7513,7413,7313,7213,7113,7013,5013,5113,
        5213,5313,5413,5513,5613,5713,5813,7812,7712,7612,7512,
        7412,7312,7212,7112,7012,5012,5112,5212,5312,5412,5512,
        5612,5712,5812,7811,7711,7611,7511,7411,7311,7211,7111,
        7011,5011,5111,5211,5311,5411,5511,5611,5711,5811,7810,
        7710,7610,7510,7410,7310,7210,7110,7010,5010,5110,5210,
        5310,5410,5510,5610,5710,5810,7809,7709,7609,7509,7409,
        7309,7209,7109,7009,5009,5109,5209,5309,5409,5509,5609,
        5709,5809,7808,7708,7608,7508,7408,7308,7208,7108,7008,
        5008,5108,5208,5308,5408,5508,5608,5708,5808,7807,7707,
        7607,7507,7407,7307,7207,7107,7007,5007,5107,5207,5307,
        5407,5507,5607,5707,5807,7806,7706,7606,7506,7406,7306,
        7206,7106,7006,5006,5106,5206,5306,5406,5506,5606,5706,
        5806,7805,7705,7605,7505,7405,7305,7205,7105,7005,5005,
        5105,5205,5305,5405,5505,5605,5705,5805,7804,7704,7604,
        7504,7404,7304,7204,7104,7004,5004,5104,5204,5304,5404,
        5504,5604,5704,5804,7803,7703,7603,7503,7403,7303,7203,
        7103,7003,5003,5103,5203,5303,5403,5503,5603,5703,5803,
        7802,7702,7602,7502,7402,7302,7202,7102,7002,5002,5102,
        5202,5302,5402,5502,5602,5702,5802,7801,7701,7601,7501,
        7401,7301,7201,7101,7001,5001,5101,5201,5301,5401,5501,
        5601,5701,5801,7800,7700,7600,7500,7400,7300,7200,7100,
        7000,5000,5100,5200,5300,5400,5500,5600,5700,5800
    ]

    for i_hist_type in [2,3,4]:
        if i_hist_type==2:
            clim_dir=clim_ctd_
            clim_type="ctd_"
            prefix=config["HISTORICAL_CTD_PREFIX"]
        elif i_hist_type==4:
            clim_dir=clim_argo_
            clim_type="argo_"
            prefix=config["HISTORICAL_ARGO_PREFIX"]
        elif i_hist_type==3:
            clim_dir=clim_bottle_
            clim_type="bot_"
            prefix=config["HISTORICAL_BOTTLE_PREFIX"]

        prefix_str = prefix.split("/")
        clim_dir = os.path.sep.join([config["HISTORICAL_DIRECTORY"], os.path.sep.join(prefix_str[:-1])])

        for file in Path(clim_dir).iterdir():
            file_name = file.name
            if len(file_name) < len(clim_type):
                continue
            if file_name.startswith(clim_type):
                box_no = file_name.replace(clim_type, "").replace(".mat", "")
                box_no = float(box_no)
                i_box=(la_wmo_boxes[:,1]==box_no)
                la_wmo_boxes[i_box,i_hist_type - 1] = 1

    savemat(la_wmo_boxes_file,'la_wmo_boxes')
