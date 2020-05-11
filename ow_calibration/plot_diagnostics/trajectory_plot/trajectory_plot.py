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


import geopandas as gdp
import matplotlib.pyplot as plt
import numpy as np



def trajectory_plot(bath, reef, floats, climatology, float_name):

    # load in the coastline data
    coastline = "data/constants/coastline/ne_10m_coastline.shp"
    map_coast = gdp.read_file(coastline)
    traj_map = map_coast.plot(color='black', label='coastline')

    # if wanted, load in bathymetric data and plot it
    if bath == 1:
        bathymetry0 = "data/constants/bathymetry/ne_10m_bathymetry_L_0.shp"
        bathymetry200 = "data/constants/bathymetry/ne_10m_bathymetry_K_200.shp"
        bathymetry1000 = "data/constants/bathymetry/ne_10m_bathymetry_J_1000.shp"
        bathymetry2000 = "data/constants/bathymetry/ne_10m_bathymetry_I_2000.shp"
        bathymetry3000 = "data/constants/bathymetry/ne_10m_bathymetry_H_3000.shp"
        bathymetry4000 = "data/constants/bathymetry/ne_10m_bathymetry_G_4000.shp"
        bathymetry5000 = "data/constants/bathymetry/ne_10m_bathymetry_F_5000.shp"
        bathymetry6000 = "data/constants/bathymetry/ne_10m_bathymetry_E_6000.shp"
        bathymetry7000 = "data/constants/bathymetry/ne_10m_bathymetry_D_7000.shp"
        bathymetry8000 = "data/constants/bathymetry/ne_10m_bathymetry_C_8000.shp"
        bathymetry9000 = "data/constants/bathymetry/ne_10m_bathymetry_B_9000.shp"
        bathymetry10000 = "data/constants/bathymetry/ne_10m_bathymetry_A_10000.shp"
        map_bath0 = gdp.read_file(bathymetry0)
        map_bath200 = gdp.read_file(bathymetry200)
        map_bath1000 = gdp.read_file(bathymetry1000)
        map_bath2000 = gdp.read_file(bathymetry2000)
        map_bath3000 = gdp.read_file(bathymetry3000)
        map_bath4000 = gdp.read_file(bathymetry4000)
        map_bath5000 = gdp.read_file(bathymetry5000)
        map_bath6000 = gdp.read_file(bathymetry6000)
        map_bath7000 = gdp.read_file(bathymetry7000)
        map_bath8000 = gdp.read_file(bathymetry8000)
        map_bath9000 = gdp.read_file(bathymetry9000)
        map_bath10000 = gdp.read_file(bathymetry10000)
        traj_map = map_bath0.plot(ax=traj_map, color='#BEBEBE', label='>200m', linewidth=2)
        traj_map = map_bath200.plot(ax=traj_map, color='#B8B8B8', linewidth=2)
        traj_map = map_bath1000.plot(ax=traj_map, color='#B0B0B0', label='1000m', linewidth=2)
        traj_map = map_bath2000.plot(ax=traj_map, color='#A9A9A9')
        traj_map = map_bath3000.plot(ax=traj_map, color='#A8A8A8')
        traj_map = map_bath4000.plot(ax=traj_map, color='#A0A0A0')
        traj_map = map_bath5000.plot(ax=traj_map, color='#989898')
        traj_map = map_bath6000.plot(ax=traj_map, color='#909090', label='6000m')
        traj_map = map_bath7000.plot(ax=traj_map, color='#888888')
        traj_map = map_bath8000.plot(ax=traj_map, color='#808080')
        traj_map = map_bath9000.plot(ax=traj_map, color='#787878')
        traj_map = map_bath10000.plot(ax=traj_map, color='#707070')

    # if we want reef data, load it in and plot it
    if reef == 1:
        reef = "data/constants/reefs/ne_10m_reefs.shp"
        map_reef = gdp.read_file(reef)
        traj_map = map_reef.plot(ax=traj_map, color='green', label='reef')

    # set up the latitude and longitude data
    geo_floats = gdp.GeoDataFrame(floats,
                                  geometry=gdp.points_from_xy(floats.Longitude,
                                                              floats.Latitude))
    geo_climatology = gdp.GeoDataFrame(climatology,
                                       geometry=gdp.points_from_xy(climatology.Longitude,
                                                                   climatology.Latitude))

    traj_map = geo_floats.plot(ax=traj_map, color='red', marker="+", label='profile')
    geo_climatology.plot(ax=traj_map, color='#00008B', marker="s", markersize=12, label='climatology')
    plt.plot(floats['Longitude'], floats['Latitude'], color='red', linestyle='-')
    plt.title(("Locations of float " + float_name + " with historical data"))
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.axis([-180, 180, -90, 90])
    plt.xlim(np.min(climatology['Longitude']) - 10, np.max(climatology['Longitude']) + 10)
    plt.ylim(np.min(climatology['Latitude']) - 10, np.max(climatology['Latitude']) + 10)

    for i, txt in enumerate(floats['number']):
        plt.annotate(txt, (floats['Longitude'][i], floats['Latitude'][i]))

    plt.legend(loc='center left', bbox_to_anchor=(0.85, 0.85))
    plt.draw()
