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

from ow_calibration.load_configuration.load_configuration import load_configuration

def update_salinity_mapping(float_dir, float_name, config):
    """
    Calculates values needed for analysis. Saves them to memory to use later

    :param float_dir: directory where the float being analysed is stored
    :param float_name: name of the float being analysed
    :param config: configuration settings set by the user
    :return: Nothing, but does save values
    """

    # Get float file name
    filename = config['FLOAT_SOURCE_DIRECTORY'] + float_dir + float_name + config['FLOAT_SOURCE_POSTFIX']
    print(filename)

