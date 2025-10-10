"""Configuration loader"""

import collections
import os

import numpy as np
from scipy.io import loadmat, savemat

# pylint: disable=too-many-locals
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
# pylint: disable=too-many-nested-blocks
# pylint: disable=invalid-name
# pylint: disable=fixme
def set_calseries(float_dir, float_name, system_config):
    """Set the calseries parameters for analysis and line fitting

    Parameters
    ----------
    use_pres_gt: pressure greater than
    use_theta_gt: ptmp greater than
    use_theta_lt: ptmp less than
    use_pres_lt: pressure less than
    float_dir: location of float
    float_name: float source name
    system_config: configuration settings

    Returns:
    -------
    Nothing, but save parameters
    """
    # load float source data
    float_source = loadmat(
        os.path.sep.join(
            [system_config["FLOAT_SOURCE_DIRECTORY"], float_dir, float_name + system_config["FLOAT_SOURCE_POSTFIX"]]
        )
    )

    profile_no = float_source["PROFILE_NO"].flatten()
    no_profiles = profile_no.__len__()

    # Check if we already have a calseries file

    calseries_filename = os.path.sep.join(
        [
            system_config["FLOAT_CALIB_DIRECTORY"],
            float_dir,
            system_config["FLOAT_CALSERIES_PREFIX"] + float_name + system_config["FLOAT_CALIB_POSTFIX"],
        ]
    )

    # if we already have a calseries file, use those values. Else, use new ones
    try:
        calseries_data = loadmat(calseries_filename)
        breaks = calseries_data["breaks"]
        max_breaks = calseries_data["max_breaks"]
        calseries = calseries_data["calseries"].flatten()
        calib_profile_no = calseries_data["calib_profile_no"].flatten()
        use_theta_lt = calseries_data["use_theta_lt"]
        use_theta_gt = calseries_data["use_theta_gt"]
        use_pres_lt = calseries_data["use_pres_lt"]
        use_pres_gt = calseries_data["use_pres_gt"]

        # use percent may not exist, as it was added later
        try:
            use_percent_gt = calseries_data["use_percent_gt"]

        except NameError:
            use_percent_gt = 0.5

        print("Using parameters found in ", calseries_filename, "\nTo use new parameters, delete this file")

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

    sal = float_source["SAL"]
    temp = float_source["TEMP"]
    pres = float_source["PRES"]

    for i in range(no_profiles):
        sal_nan = np.argwhere(~np.isnan(sal[:, i]))
        temp_nan = np.argwhere(~np.isnan(temp[:, i]))
        pres_nan = np.argwhere(~np.isnan(pres[:, i]))

        # if no good data for this profile, remove it from calseries
        if sal_nan.__len__() == 0 or temp_nan.__len__() == 0 or pres_nan.__len__() == 0:
            calseries[i] = 0

    savemat(
        calseries_filename,
        {
            "breaks": system_config["BREAKS"],
            "max_breaks": system_config["MAX_BREAKS"],
            "calseries": system_config["CALSERIES"],
            "calib_profile_no": system_config["CALIB_PROFILE_NO"],
            "use_theta_lt": system_config["USE_THETA_LT"],
            "use_theta_gt": system_config["USE_THETA_GT"],
            "use_pres_lt": system_config["USE_PRES_LT"],
            "use_pres_gt": system_config["USE_PRES_GT"],
            "use_percent_gt": system_config["USE_PERCENT_GT"],
        },
    )
