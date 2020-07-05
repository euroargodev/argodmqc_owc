"""
-----Theta Salinity Plot-----
Written by: Edward Small
When: 10/05/2020
Function for plotting salinity curve aainst potential temperature of all the
data used in the analysis
For information on how to use this file, check the README at either:
https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
from ow_calibration.find_10thetas.find_10thetas import find_10thetas


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
def sal_var_plot(levels, sal, pres, ptmp, map_sal, map_sal_errors,
                 map_ptmp, cal_sal, cal_sal_errors, boundaries):
    """
    Create the salinity variance plot for each level
    :param map_ptmp: mapped potential temperature
    :param levels: number of levels to plot
    :param sal: float salinity
    :param pres: float pressure
    :param ptmp: float potential temperature
    :param map_sal: mapped salinity
    :param map_sal_errors: mapped salinity errors
    :param cal_sal: calibrated salinity
    :param cal_sal_errors: calibrated salinity errors
    :param boundaries: pressure and temperature boundaries
    :return: Nothing
    """

    use_theta_lt = boundaries[0]
    use_theta_gt = boundaries[1]
    use_pres_lt = boundaries[2]
    use_pres_gt = boundaries[3]
    use_percent_gt = boundaries[4]

    thetas = find_10thetas(copy.deepcopy(sal),
                           copy.deepcopy(ptmp),
                           copy.deepcopy(pres),
                           copy.deepcopy(map_ptmp),
                           use_theta_lt, use_theta_gt,
                           use_pres_lt, use_pres_gt,
                           use_percent_gt)
    print("called")