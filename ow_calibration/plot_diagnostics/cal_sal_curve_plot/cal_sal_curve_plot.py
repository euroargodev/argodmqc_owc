"""
-----Theta Salinity Plot-----

Written by: Edward Small
When: 10/05/2020

Function for plotting the calibrated salinity curve

For information on how to use this file, check the README at either:
https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
def cal_sal_curve_plot(sal, cal_sal, sta_sal, pcond_factor, pcond_factor_err, float_name):
    """
    Create the calibrated salinity curve plot
    :param sal: float salinity
    :param cal_sal: calibrated salinity
    :param sta_sal: mean difference between salinity and calculated salinity
    :param pcond_factor: slope
    :param pcond_factor_err: slope error
    :param float_name: name of the float
    :return: Nothing
    """

