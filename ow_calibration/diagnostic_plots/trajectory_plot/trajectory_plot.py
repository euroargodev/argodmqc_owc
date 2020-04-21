import pandas as pd
import geopandas as gdp
import matplotlib.pyplot as plt


float1 = pd.DataFrame(
    {'climatology': ['0', '1', '2', '3'],
     'Latitude': [-52.5, -52.3, -52.4, -52.5],
     'Longitude': [75, 74.9, 74.9, 75.2]})
df = pd.DataFrame(
    {'climatology': ['1', '2', '3', '4', '5'],
     'Latitude': [-51.5, -52.2, -52.1, -51, -52.6],
     'Longitude': [75, 74.2, 75.8, 75.8, 74.9]})
coastline = "../../../data/constants/coastline/ne_10m_coastline.shp"
reef = "../../../data/constants/reefs/ne_10m_reefs.shp"
bathymetry200 = "../../../data/constants/bathymetry/ne_10m_bathymetry_K_200.shp"
bathymetry1000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_J_1000.shp"
bathymetry2000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_I_2000.shp"
bathymetry3000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_H_3000.shp"
bathymetry4000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_G_4000.shp"
bathymetry5000 = "../../../data/constants/bathymetry/ne_10m_bathymetry_F_5000.shp"
floats = gdp.GeoDataFrame(df, geometry=gdp.points_from_xy(df.Longitude, df.Latitude))
float2 = gdp.GeoDataFrame(float1, geometry=gdp.points_from_xy(float1.Longitude, float1.Latitude))

map_df = gdp.read_file(coastline)
map_ref = gdp.read_file(reef)
map_bath200 = gdp.read_file(bathymetry200)
map_bath1000 = gdp.read_file(bathymetry1000)
map_bath2000 = gdp.read_file(bathymetry2000)
map_bath3000 = gdp.read_file(bathymetry3000)
map_bath4000 = gdp.read_file(bathymetry4000)
map_bath5000 = gdp.read_file(bathymetry5000)

ax = map_df.plot(color='black')
ax = map_ref.plot(ax=ax, color='green')
ax = map_bath200.plot(ax=ax, color='lightblue')
ax = map_bath1000.plot(ax=ax, color='blue')
ax = map_bath2000.plot(ax=ax, color='darkblue')
ax = map_bath3000.plot(ax=ax, color='#00008b')
ax = map_bath4000.plot(ax=ax, color='#000058')
ax = map_bath5000.plot(ax=ax, color='#00003f')

#floats.plot(ax=ax, color='red')
#float2.plot(ax=ax, color='black')

plt.show()