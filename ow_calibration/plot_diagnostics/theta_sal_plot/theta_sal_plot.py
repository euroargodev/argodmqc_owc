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
import matplotlib.pylab as pl
import numpy as np


def theta_sal_plot(sal, theta, map_sal, map_theta, map_errors, title='uncalibrated'):
    """
    Create the salinity theta curve
    :param title: Addition to the title
    :param sal: float salinity
    :param theta: float potential temperature
    :param map_sal: mapped salinity
    :param map_theta: mapped potential temperature
    :param map_errors: mapped salinity errors
    :return: Nothing
    """

    fig, ax = plt.subplots()
    color_n = sal.__len__()
    colors = pl.cm.jet(np.linspace(0, 1, color_n))

    for i in range(sal.__len__()):
        # plot salinities
        ax = plt.plot(sal[i], theta[i], color=colors[i])

    plt.title(title + " float data with mapped salinity and objective errors")
    plt.xlabel("Salinity (PSS-78)")
    plt.ylabel(r"$\theta$ $^\circ$C")
    plt.show()
