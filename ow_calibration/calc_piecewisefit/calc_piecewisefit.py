"""
-----Calculate Piece wise fit-----

Written by: Annie Wong
When: xx/10/2008
Converted to python by: Edward Small
When: 16/06/2020

Calculate the fit of each break and calibrate salinities

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
import scipy.io as scipy

def calc_piecewisefit(float_dir, float_name, system_config):
    """
    calibrate salinities
    :param float_dir: float directory name
    :param float_name: name of float
    :param system_config: configuration parameter set up
    :return: Nothing, save output
    """

    # load in the source data

    float_source_data = scipy.loadmat(system_config['FLOAT_SOURCE_DIRECTORY'] +
                                      float_dir + float_name + system_config['FLOAT_SOURCE_POSTFIX'])

    lat = float_source_data['LAT']
    long = float_source_data['LONG']
    dates = float_source_data['DATES']
    sal = float_source_data['SAL']
    ptmp = float_source_data['PTMP']
    pres = float_source_data['PTMP']
    profile_no = float_source_data['PROFILE_NO']
    x_in = np.tile(profile_no, (10,1))

    # load in the mapped data
    float_mapped_data = scipy.loadmat(system_config['FLOAT_MAPPED_DIRECTORY'] +
                                      float_dir + system_config['FLOAT_MAPPED_PREFIX'] +
                                      float_name + system_config['FLOAT_MAPPED_POSTFIX'])

    mapped_sal = float_mapped_data['la_mapped_sal']
    mapsalerror = float_mapped_data['la_mapsalerrors']
    mapped_ptmp = float_mapped_data['la_ptmp']