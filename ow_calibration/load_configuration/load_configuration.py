"""
-----Load Configuration File-----

Written by: Edward Small
When: 18/09/2019
Converted to python by: N/A
When: N/A

Contains the configuration settings of the user for an analysis.
The variables are saved in a dictionary.
This means that values are accessible both via index and via key.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""


def load_configuration():
    """
    Takes in no arguments
    :return: A dictionary containing mapped pairs of
    {PARAMETER_NAME: parameter_value}
    """

    return {

        # ===============================
        #
        #    Climatology Data Input Paths

        'HISTORICAL_DIRECTORY': "../../data/climatology",
        'HISTORICAL_CTD_PREFIX': "/historical_ctd/ctd_",
        'HISTORICAL_BOTTLE_PREFIX': "/historical_bot/bot_",
        'HISTORICAL_ARGO_PREFIX': "/historical_argo/argo_",

        # ===============================
        #
        #    Float Input Path
        #

        'FLOAT_SOURCE_DIRECTORY': "data/float_source",
        'FLOAT_SOURCE_POSTFIX': ".mat",

        # ===============================
        #
        #    Mapping Output Path
        #

        'FLOAT_MAPPED_DIRECTORY': "/home/awong/argo/OW/data/float_mapped/",
        'FLOAT_MAPPED_PREFIX': "map_",
        'FLOAT_MAPPED_POSTFIX': ".mat",

        # ===============================
        #
        #    Calibration Output Path
        #

        'FLOAT_CALIB_DIRECTORY': "/home/awong/argo/OW/data/float_calib/",
        'FLOAT_CALIB_PREFIX': "cal_",
        'FLOAT_CALSERIES_PREFIX': "calseries_",
        'FLOAT_CALIB_POSTFIX': ".mat",

        # ===============================
        #
        #    Diagnostic Plots Output Path
        #

        'FLOAT_PLOTS_DIRECTORY': "/home/awong/argo/OW/data/float_plots/",

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
        'MAPSCALE_PHI_LARGE': 0.5,
        'MAPSCALE_PHI_SMALL': 0.1,

        # temporal decorrelation scale, in years
        'MAPSCALE_AGE_LARGE': 20,
        'MAPSCALE_AGE_SMALL': 10,

        # exclude the top xxx dbar of the water column
        'MAP_P_EXCLUDE': 200,

        # only use historical data that are within +/- yyy dbar from float data
        'MAP_P_DELTA': 250
    }
