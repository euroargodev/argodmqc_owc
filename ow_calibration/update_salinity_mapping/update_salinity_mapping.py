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

    # Load precalculated mapped data -------------------------------------

    try:
        float_mapped_filename = scipy.loadmat(config['FLOAT_MAPPED_DIRECTORY'] + float_dir +
                                              config['FLOAT_MAPPED_PREFIX'] + float_name +
                                              config['FLOAT_MAPPED_POSTFIX'])

    except FileNotFoundError:

        print("__________________________________________________________")
        print("No precaulcated data\n")

        float_mapped_filename = np.nan
