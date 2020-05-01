import pandas as pd
import geopandas as gdp
import numpy as np
import matplotlib.pyplot as plt

hist_lat = np.ones((25)) * + (np.random.rand(25)*2 - 12)
hist_long = np.ones((25)) * + (np.random.rand(25)*2 + 144)
climatology = range(25)
print(hist_long)

float1 = pd.DataFrame(
    {'climatology': climatology,
     'Latitude': hist_lat,
     'Longitude': hist_long})
df = pd.DataFrame(
    {'climatology': ['1', '2', '3', '4', '5'],
     'Latitude': [-10.5, -10.2, -11.1, -11, -11.6],
     'Longitude': [146, 146.2, 144.8, 143.8, 143.9]})
coastline = "../../../data/constants/coastline/ne_10m_coastline.shp"
reef = "../../../data/constants/reefs/ne_10m_reefs.shp"
bathymetry0 = "../../../data/constants/bathymetry/ne_10m_bathymetry_L_0.shp"
bathymetry200 = "../../../data/constants/bathymetry/ne_10m_bathymetry_K_200.shp"
bathymetry1000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_J_1000.shp"
bathymetry2000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_I_2000.shp"
bathymetry3000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_H_3000.shp"
bathymetry4000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_G_4000.shp"
bathymetry5000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_F_5000.shp"
bathymetry6000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_E_6000.shp"
bathymetry7000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_D_7000.shp"
bathymetry8000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_C_8000.shp"
bathymetry9000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_B_9000.shp"
bathymetry10000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_A_10000.shp"
floats = gdp.GeoDataFrame(df, geometry=gdp.points_from_xy(df.Longitude, df.Latitude))
float2 = gdp.GeoDataFrame(float1, geometry=gdp.points_from_xy(float1.Longitude, float1.Latitude))



map_coast = gdp.read_file(coastline)
map_reef = gdp.read_file(reef)
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

ax = map_coast.plot(color='black')
ax = map_reef.plot(ax=ax, color='green')
ax = map_bath0.plot(ax=ax, color='#A8ADFF')
ax = map_bath200.plot(ax=ax, color='#626BFE')
ax = map_bath1000.plot(ax=ax, color='#4E59FE')
ax = map_bath2000.plot(ax=ax, color='#414CFF')
ax = map_bath3000.plot(ax=ax, color='#2E3EFE')
ax = map_bath4000.plot(ax=ax, color='#1826FD')
ax = map_bath5000.plot(ax=ax, color='#000FFF')
ax = map_bath6000.plot(ax=ax, color='#000CD3')
ax = map_bath7000.plot(ax=ax, color='#000BB5')
ax = map_bath8000.plot(ax=ax, color='#00088D')
ax = map_bath9000.plot(ax=ax, color='#000668')
ax = map_bath10000.plot(ax=ax, color='#000442')


ax = floats.plot(ax=ax, color='red', marker="+")
ax = float2.plot(ax=ax, color='#00008B', marker="s", markersize=12)
plt.plot(df['Longitude'], df['Latitude'], color='red', linestyle='-')
plt.title("Locations of float 1904865 with historical data")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.axis([-180, 180, -90, 90])
plt.xlim(np.min(df['Longitude']) - 2, np.max(df['Longitude']) + 2)
plt.ylim(np.min(df['Latitude']) - 2, np.max(df['Latitude']) + 2)


plt.show()