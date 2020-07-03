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
def trajectory_plot(mapped_data, float_long, float_lat,
                    levels=[-8000, -6000, -4000, -2000, -1000, -800, -600, -400, -200, 0],
                    bathy=False, cmap=None):
    """
    Reveal trajectory diagnostic plot
    :param levels: levels to plot bathymetry data
    :param cmap: shading for bathymetry ('block' or 'shade')
    :param bathy: Boolean to plot bathymetry
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
                 levels=[50, 10000], colors='bisque', antialiased=False, linestyles='solid')

    # plot climatology
    ax.plot(mapped_data[:, 0], mapped_data[:, 1],
            color='#FF6347', marker='s',
            linestyle='None', markersize=2,
            transform=projection, label="Climatology"
            )

    # plot float data
    ax.plot(float_long, float_lat,
            color='red', marker='x',
            linestyle='-', markersize=4,
            transform=projection, label="Float Profiles"
            )

    # plot bathymetry
    if bathy:
        levels = levels
        con = plt.contour(topo.coords['longitude'], topo.coords['latitude'], topo.variables['topography'],
                          levels=levels, colors='black',
                          linestyles='solid', linewidths=0.3)
        plt.clabel(con, inline=1, fontsize=7)

        # shade, if wanted
        if cmap == "shade":
            pc = topo.topography.plot.pcolormesh(
                ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False, cmap='winter',
                levels=levels, antialiased=True
            )
            plt.colorbar(
                pc, label="metres", orientation="horizontal", aspect=50, pad=0.1, shrink=0.6
            )

        elif cmap == "block":
            pc = plt.contourf(topo.coords['longitude'], topo.coords['latitude'], topo.variables['topography'],
                              levels=[-8000, -6000, -4000, -2000, -1000, -800, -600, -400, -200, 0],
                              cmap='winter', linestyles='solid', linewidths=0.3, antialiased=True)
            plt.colorbar(
                pc, label="metres", orientation="horizontal", aspect=50, pad=0.1, shrink=0.6
            )

    plt.show()
