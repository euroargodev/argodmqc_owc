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
import cartopy.feature as feature
import harmonica as hm


# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
def trajectory_plot(mapped_data, float_long, float_lat):
    """
    Reveal trajectory diagnostic plot
    :param mapped_data: historical climatology data (mapped)
    :param float_lat: float latitude
    :param float_long: float longitude
    :return: Nothing
    """

    # set up plot and axis
    fig, ax = plt.subplots()
    projection = ccrs.PlateCarree()
    ax = plt.axes(projection=projection)
    topo = hm.datasets.fetch_topography_earth()

    # plot land
    plt.contourf(topo.coords['longitude'], topo.coords['latitude'], topo.variables['topography'],
                 levels=[50, 10000], colors='bisque', antialiased=False,
                 linestyles='solid', linewidths=0.3)

    # plot climatology
    ax.plot(mapped_data[:, 0], mapped_data[:, 1],
            color='#FF6347', marker='s',
            linestyle='None', markersize=2,
            transform=projection, label="Climatology"
            )

    con = plt.contour(topo.coords['longitude'], topo.coords['latitude'], topo.variables['topography'],
                      levels=[-8000, -6000, -4000, -2000, -1000, -800, -600, -400, -200], colors='black',
                      linestyles='solid', linewidths=0.3)

    plt.show()
