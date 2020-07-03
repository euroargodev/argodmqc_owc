"""
-----trajectory plot-----

Written by: Edward Small
When: 10/05/2020

function for plotting locations of all the data used in the analysis, including:

historical climatology
profile locations and order

Can also plot reef data and bathymetry by passing in a 1 (True) into the function

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""


import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import harmonica as hm


# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
def trajectory_plot(mapped_data, float_lat, float_long):
    """
    Reveal trajectory diagonstic plot
    :param mapped_data: historical climatology data (mapped)
    :param float_lat:
    :param float_long:
    :return: Nothing
    """

    fig, ax = plt.subplots()
    projection = ccrs.PlateCarree()
    ax = plt.axes(projection=projection)

    plt.show()

