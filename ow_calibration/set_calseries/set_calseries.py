"""
-----noise_variance-----

Written by: Breck Owens
When: xx/09/2006
Converted to python by: Edward Small
When: 12/05/2020

Set the calseries parameters for analysis and line fitting

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import scipy.io as scipy


def set_calseries(float_dir, float_name, system_config):
    """
    Set up the parameters for line fitting
    :param float_dir:
    :param float_name:
    :param system_config:
    :return: Nothing, but save parameters
    """

    # load float source data
    float_source = scipy.loadmat(system_config['FLOAT_SOURCE_DIRECTORY'] +
                                 float_dir + float_name +
                                 system_config['FLOAT_SOURCE_POSTFIX'])

    profile_no = float_source['PROFILE_NO']
    no_profiles = profile_no.shape[1]

    # Check if we already have a calseries file

    calseries_filename = (system_config['FLOAT_CALIB_DIRECTORY'] +
                          float_dir +
                          system_config['FLOAT_CALSERIES_PREFIX'] +
                          float_name +
                          system_config['FLOAT_CALIB_POSTFIX'])