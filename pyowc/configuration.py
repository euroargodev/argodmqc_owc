""" Configuration loader
"""
import collections
import os

import numpy as np
from scipy.io import loadmat, savemat


def print_cfg(config):
    """ Return a string of the configuration dictionnary """
    cfg_str = []
    for k in config:
        cfg_str.append("%30s: %s\n" % (k, str(config[k])))
    return "\n" + "".join(cfg_str)


def load():
    """
    Takes in no arguments
    :return: A dictionary containing mapped pairs of
    {PARAMETER_NAME: parameter_value}
    """

    cfg = {

        # ===============================
        #
        #    Climatology Data Input Paths

        'HISTORICAL_DIRECTORY': "data/climatology",
        'HISTORICAL_CTD_PREFIX': "/historical_ctd/ctd_",
        'HISTORICAL_BOTTLE_PREFIX': "/historical_bot/bot_",
        'HISTORICAL_ARGO_PREFIX': "/historical_argo/argo_",

        # ===============================
        #
        #    Float Input Path
        #

        'FLOAT_SOURCE_DIRECTORY': "data/float_source/",
        'FLOAT_SOURCE_POSTFIX': ".mat",

        # ===============================
        #
        #    Mapping Output Path
        #

        'FLOAT_MAPPED_DIRECTORY': "data/float_mapped/",
        'FLOAT_MAPPED_PREFIX': "map_",
        'FLOAT_MAPPED_POSTFIX': ".mat",

        # ===============================
        #
        #    Calibration Output Path
        #

        'FLOAT_CALIB_DIRECTORY': "data/float_calib/",
        'FLOAT_CALIB_PREFIX': "cal_",
        'FLOAT_CALSERIES_PREFIX': "calseries_",
        'FLOAT_CALIB_POSTFIX': ".mat",

        # ===============================
        #
        #    Diagnostic Plots Output Path
        #

        'FLOAT_PLOTS_DIRECTORY': "data/float_plots/",
        'FLOAT_PLOTS_FORMAT': "eps",
        # ===============================
        #
        #    Constants File Path
        #

        'CONFIG_DIRECTORY': "data/constants/",
        'CONFIG_COASTLINES': "coastdat.mat",
        'CONFIG_WMO_BOXES': "wmo_boxes.mat",
        'CONFIG_SAF': "TypicalProfileAroundSAF.mat",

        # ===============================
        #
        #    Objective Mapping Parameters
        #

        # max number of historical casts used in objective mapping
        'CONFIG_MAX_CASTS': 300,

        # 1=use PV constraint, 0=don't use PV constraint, in objective mapping
        'MAP_USE_PV': 0,

        # 1=use SAF separation criteria, 0=don't use SAF separation criteria, in objective mapping
        'MAP_USE_SAF': 0,

        # spatial decorrelation scales, in degrees
        'MAPSCALE_LONGITUDE_LARGE': 8,
        'MAPSCALE_LONGITUDE_SMALL': 4,
        'MAPSCALE_LATITUDE_LARGE': 4,
        'MAPSCALE_LATITUDE_SMALL': 2,

        # cross-isobath scales, dimensionless, see BS(2005)
        'MAPSCALE_PHI_LARGE': 0.1,
        'MAPSCALE_PHI_SMALL': 0.02,

        # temporal decorrelation scale, in years
        'MAPSCALE_AGE_LARGE': 20,
        'MAPSCALE_AGE_SMALL': 5,

        # exclude the top xxx dbar of the water column
        'MAP_P_EXCLUDE': 100,

        # only use historical data that are within +/- yyy dbar from float data
        'MAP_P_DELTA': 200,

        # ===============================
        #
        #    Plotting Parameters
        #

        # Theta bounds for salinity anomaly plot
        'THETA_BOUNDS': [[0, 5], [5, 20]]
    }
    return collections.OrderedDict(sorted(cfg.items()))


# pylint: disable=too-many-locals
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
# pylint: disable=too-many-nested-blocks
# pylint: disable=invalid-name
# pylint: disable=fixme
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
    float_source = loadmat(os.path.sep.join([system_config['FLOAT_SOURCE_DIRECTORY'], float_dir,
                                             float_name + system_config['FLOAT_SOURCE_POSTFIX']]))

    profile_no = float_source['PROFILE_NO'].flatten()
    no_profiles = profile_no.__len__()

    # Check if we already have a calseries file

    calseries_filename = (os.path.sep.join([system_config['FLOAT_CALIB_DIRECTORY'], float_dir,
                                            system_config['FLOAT_CALSERIES_PREFIX'] +
                                            float_name +
                                            system_config['FLOAT_CALIB_POSTFIX']]))

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
        # calseries = np.concatenate((
        #   np.ones((1, 18)).flatten(),
        #    2 * np.ones((1, no_profiles - 18)).flatten()))
        calib_profile_no = profile_no
        use_percent_gt = 0.5
        use_theta_lt = []
        use_theta_gt = []
        use_pres_lt = []
        use_pres_gt = []

    # ensure values are in a realistic range

    if use_theta_lt.__len__() > 1:
        print("More than one potential temperature boundary used, removing boundary...")
        use_theta_lt = []

    if use_theta_gt.__len__() > 1:
        print("More than one potential temperature boundary used, removing boundary...")
        use_theta_gt = []

    if use_pres_lt.__len__() > 1:
        print("More than one pressure boundary used, removing boundary...")
        use_pres_lt = []

    if use_pres_gt.__len__() > 1:
        print("More than one pressure boundary used, removing boundary...")
        use_pres_gt = []

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
