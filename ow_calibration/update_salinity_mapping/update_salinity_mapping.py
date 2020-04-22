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

import scipy.io as scipy
import numpy as np
import time


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

    wmo_boxes = scipy.loadmat(config['CONFIG_DIRECTORY'] + config['CONFIG_WMO_BOXES'])['la_wmo_boxes']
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

        print("No precaulcated data")
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
        start_time = time.perf_counter()

        print("UPDATE_SALINITY_MAPPING: Working on profile ", i)

        # Get the current profile being worked on
        missing_profile = missing_profile_index[i]

        # append profile numbers
        la_profile_no = np.insert(la_profile_no, profile_index, profile_no[missing_profile])
        # Construct elements for this profile

        la_ptmp = np.hstack((la_ptmp, np.nan * np.ones((float_level_count, 1))))
        print(la_ptmp)
        input("----------------------")

        profile_index += 1
