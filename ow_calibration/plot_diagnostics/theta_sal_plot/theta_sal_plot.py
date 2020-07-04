"""
-----Theta Salinity Plot-----

Written by: Edward Small
When: 10/05/2020

Function for plotting Theta salinity of all the data used in the analysis

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import matplotlib.pyplot as plt

def theta_sal_plot(sal, theta, map_sal, map_theta, map_errors):
    """
    Create the salinity theta curve
    :param sal: float salinity
    :param theta: float potential temperature
    :param map_sal: mapped salinity
    :param map_theta: mapped potential temperature
    :param map_errors: mapped salinity errors
    :return: Nothing
    """


