"""
-----profile Plot-----

Written by: Edward Small
When: 10/05/2020

Plots useful data, such as salinity variance, pressure curves, etc

For information on how to use this file, check the README at either:
https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import matplotlib.pyplot as plt


def t_s_profile_plot(sal, ptmp, pres, sal_var, theta_levels, tlevels, plevels, float_name):
    """
    plots profile plots
    :param sal: salinity
    :param ptmp: potential temperature
    :param pres: pressure
    :param sal_var: salininty variance on theta levels
    :param theta_levels: theta levels
    :param tlevels: temp at theta levels
    :param plevels: pressure at theta levels
    :param float_name: name of the float
    :return: Nothing
    """

    plt.figure(1)

    # plot t-s profile
    plt.subplot(222)

    plt.plot(sal, ptmp, color='b', linestyle='-', linewidth=0.5)
    plt.xlabel("PSS-78")
    plt.title(" OW chosen levels - " + float_name)

    for i in tlevels:
        plt.axhline(y=i, color=(0, 1, 0), linestyle='-')

    #plot s_variance on t
    plt.subplot(221)

    plt.plot(sal_var, theta_levels, color='b', linewidth=0.5)
    plt.xlabel("Salinity variance")
    plt.ylabel(r"Potential temp ($^{\circ}$C")
    plt.title("Salinity Variance on Theta")

    for i in tlevels:
        plt.axhline(y=i, color=(0, 1, 0), linestyle='-')


    plt.show()
