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
# pyline: disable=too-many-arguments
def trajectory_plot(float_name, mapped_data, float_long, float_lat,
                    depths=[-8000, -6000, -4000, -2000, -1000, -500, -200, 0],
                    bathy=False, cmap='gray', style='block'):
    """
    Reveal trajectory diagnostic plot
    :param float_name: name of the float
    :param style: shading style for color map
    :param levels: levels to plot bathymetry data
    :param cmap: shading for bathymetry ('block' or 'shade')
    :param bathy: Boolean to plot bathymetry
    :param mapped_data: historical climatology data (mapped)
    :param float_lat: float latitude
    :param float_long: float longitude
    :return: Nothing
    """

    # set up plot and axis
    projection = ccrs.PlateCarree()
    traj_plot = plt.axes(projection=projection)
    topo = hm.datasets.fetch_topography_earth()
    map_long = mapped_data[:, 0]
    map_lat = mapped_data[:, 1]

    # plot land
    plt.contour(topo.coords['longitude'], topo.coords['latitude'], topo.variables['topography'],
                 levels=[50, 10000], colors='black', antialiased=True, linestyles='solid', linewidths=0.5)
    plt.contourf(topo.coords['longitude'], topo.coords['latitude'], topo.variables['topography'],
                 levels=[50, 10000], colors='#228B22', antialiased=True, linestyles='solid')

    # plot climatology
    traj_plot.plot(map_long, map_lat,
                   color='#3633FF', marker='s',
                   linestyle='None', markersize=2,
                   transform=projection, label="Climatology"
                   )

    # plot float data
    traj_plot.plot(float_long, float_lat,
                   color='red', marker='x',
                   linestyle='-', markersize=4,
                   transform=projection, label="Float Profiles"
                   )

    # plot bathymetry
    if bathy:
        con = plt.contour(topo.coords['longitude'], topo.coords['latitude'], topo.variables['topography'],
                          levels=depths, colors='black',
                          linestyles='solid', linewidths=0.3)
        plt.clabel(con, inline=1, fontsize=7)

        # shade, if wanted
        if style == "shade":
            pc = topo.topography.plot.pcolormesh(
                ax=traj_plot, transform=ccrs.PlateCarree(), add_colorbar=False, cmap=cmap,
                levels=depths, antialiased=True
            )
            plt.colorbar(
                pc, label="metres", orientation="horizontal", aspect=50, pad=0.1, shrink=0.6
            )

        elif style == "block":
            pc = plt.contourf(topo.coords['longitude'], topo.coords['latitude'],
                              topo.variables['topography'],
                              levels=depths,
                              cmap=cmap, linestyles='solid', antialiased=True)
            plt.colorbar(
                pc, label="metres", orientation="horizontal", aspect=50, pad=0.1, shrink=0.6
            )

    # annotate float data (every 5, plus first and last float))
    color = plt.get_cmap('jet')

    for i in range(float_lat.__len__()):
        if i == 0 or i % 5 == 0 or i == float_lat.__len__() - 1:
            plt.annotate(i + 1, (float_long[i], float_lat[i]),
                         color=color((i + 1) / float_lat.__len__()))

    # Set location focus and gridlines
    traj_plot.gridlines(crs=projection, draw_labels=True)
    traj_plot.set_extent([np.min(map_long) - 10, np.max(map_long) + 10, np.min(map_lat) - 10, np.max(map_lat) + 10],
                         crs=projection)

    # add legend
    plt.legend()

    # titles and axis

    plt.ylabel("Latitude")
    plt.xlabel("Longitude")
    plt.title(float_name + " profile locations with historical data", pad=25)

    # display the plot
    plt.show()
