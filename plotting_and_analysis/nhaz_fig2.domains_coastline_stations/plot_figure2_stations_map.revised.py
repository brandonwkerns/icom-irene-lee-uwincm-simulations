import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import json

import helpers
from helpers import plot_map_background0, plot_map_background

plt.close('all')


##
## Function for plotting the map background.
##
def plot_map_background_inset(plot_area = [-180, 180, -70, 70], ax = plt.gca()
                        , projection = 'cyl', resolution = 'l'
                        , draw_meridians = np.arange(0,360,2), meridians_line_width = 0.25
                        , meridians_labels = [0,0,0,1], meridians_fontsize=10
                        , meridians_dashes = [4,2]
                        , draw_parallels = np.arange(-90,90,2), parallels_line_width = 0.25
                        , parallels_labels = [1,0,0,0], parallels_fontsize=10
                        , parallels_dashes = [4,2]
                        , coastline_line_width = 0.5
                        , countries_line_width = 0.25
                        , x_tick_rotation = 45.0
                        , y_tick_rotation = 45.0):

    map1=Basemap(projection = projection,
                 llcrnrlon = plot_area[0],
                 urcrnrlon = plot_area[1],
                 llcrnrlat = plot_area[2],
                 urcrnrlat = plot_area[3],
                 resolution = resolution, ax = ax)

    # draw lat/lon grid lines degrees.
    meridians = map1.drawmeridians(draw_meridians, linewidth = meridians_line_width
                                   , labels = meridians_labels, fontsize = meridians_fontsize
                                   , dashes = meridians_dashes)
    parallels = map1.drawparallels(draw_parallels, linewidth = parallels_line_width
                                   , labels = parallels_labels, fontsize = parallels_fontsize
                                   , dashes = parallels_dashes)
    # map1.drawcoastlines(linewidth = coastline_line_width)
    # map1.drawcountries(linewidth = countries_line_width)

    for m in meridians:
        try:
            meridians[m][1][0].set_rotation(x_tick_rotation)
        except:
            pass

    for p in parallels:
        try:
            parallels[p][1][0].set_rotation(y_tick_rotation)
        except:
            pass

    return map1


plot_area_large_map = [-105.0, -39.0, 4.0, 48.0]

plot_area = [-78.0, -71.0, 36.0, 42.2]

##
## Set up Matplotlib figure.
##
fig = plt.figure(figsize=(8,6))



##
## Map with domains and storm tracks.
##
ax1= fig.add_subplot(2,2,1) 
map1 = plot_map_background0(plot_area_large_map, ax = ax1,
                draw_meridians = np.arange(0,360,10),
                draw_parallels = np.arange(-90,90,10))


## Domains.
H_waves = helpers.plot_umwm_domain(ax1)
H_ocean = helpers.plot_hycom_domain(ax1, color='r')
H_wrf_d01 = helpers.plot_wrf_domain(ax1, 1, color='b')
H_wrf_d02_irene = helpers.plot_wrf_domain2(ax1, 2, color='b')
H_wrf_d02_lee = helpers.plot_wrf_domain_lee(ax1, 2, color='b')


## Add storm tracks.

## IBTracs Best Track data
## This was obtained from https://www.ncdc.noaa.gov/ibtracs/
##   --> Data Access --> CSV (Comma Separated Values) --> Choose which file you want.
## Current full link for North Atlantic:
##   https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.NA.list.v04r00.csv


best_track_csv_file = '/home/orca/data/best_track/IBTrACS/ibtracs.NA.list.v04r00.csv'
best_track_basin = 'NA'     # NA = North Atlantic; EP = East Pacific; WP = Western North Pacific
best_track_color = 'k'
best_track_year = 2011
best_track_linewidth = 0.7
print('Best Track: '+best_track_csv_file)

## Best Track: Irene
best_track_storm_name = 'Irene'
best_track = helpers.get_ibtracs_csv_storm_info(best_track_csv_file,best_track_basin
                                        , best_track_year, best_track_storm_name)
helpers.plot_best_track(best_track, ax1)
plt.text(-73, 14.0, 'H. Irene', fontsize=10, fontweight='bold')

## Best Track: Lee
best_track_storm_name = 'Lee'
best_track = helpers.get_ibtracs_csv_storm_info(best_track_csv_file,best_track_basin
                                        , best_track_year, best_track_storm_name)
helpers.plot_best_track(best_track, ax1)
plt.text(-94, 23.0, 'T.S. Lee', fontsize=10, fontweight='bold')

ax1.legend(handles=[H_wrf_d01, H_ocean, H_waves, H_wrf_d02_irene, H_wrf_d02_lee],
    labels=['WRFd01','HYCOM','UMWM', 'WRFd02, Irene', 'WRFd02, Lee'], fontsize=9,
    bbox_to_anchor=(1.05, 0.95), loc='lower right',ncol=2, frameon=False)

##
## Detailed Coastline Panel.
## (This was originally an inset, hence the "inset_area" labeling.)
##
inset_area = [-78.1,-72.01,34.8,40.8]

XBOX=[inset_area[0],inset_area[1],inset_area[1],inset_area[0],inset_area[0]]
YBOX=[inset_area[2],inset_area[2],inset_area[3],inset_area[3],inset_area[2]]


ax11 = fig.add_subplot(2,2,2)


## Define labels, track files, line colors, and line widths for model runs
##   in the dictionary below.

map11 = plot_map_background_inset(inset_area, ax=ax11
                           ,draw_meridians = np.arange(0,361,2)
                           ,draw_parallels = np.arange(-90,91,2))



## UWIN-CM Coast
DS = xr.open_dataset('/home/orca/bkerns/projects/coastal/model_fields/irene_and_lee_hourly_fields__2011082400.nc')
Z = DS['bathymetry_hycom'][0,:,:].copy().values
Z[~np.isfinite(Z)] = -10.0

## GOFS Coast
DS0 = xr.open_dataset('depth_GLBb0.08_09m11.nc')
Z0 = DS0['depth'][0,:,:].copy().values
Z0[~np.isfinite(Z0)] = -10.0

ax11.contour(DS0['Longitude'].values, DS0['Latitude'].values, Z0, levels=[0.0,],colors='k', linewidths=2.0)

ax11.contour(DS['longitude_hycom'][0,:,:].values,DS['latitude_hycom'][0,:,:].values, Z, levels=[0.0,],colors='r', linewidths=1.0)

## Proxy artists.
h11_0 = ax11.plot([-999, -998], [-999, -998], color='k', linewidth=2)[0]
h11_1 = ax11.plot([-999, -998], [-999, -998], color='r', linewidth=1)[0]
ax11.legend(handles=[h11_0, h11_1], labels=['GLobal HYCOM (GOFS)','UWIN-CM HYCOM'], loc='lower right', bbox_to_anchor=[1.05, 0.95], fontsize=9, frameon=False)



##
## Panels with station locations.
##

ax2 = fig.add_subplot(2,3,4)
map2 = plot_map_background(plot_area, ax = ax2)

stations_ndbc = pd.read_json('./stations_ndbc.json')
stations_ndbc.reset_index(inplace=True)

## Stations
H2 = ax2.scatter(stations_ndbc['LON'], stations_ndbc['LAT'], marker='x', s=70, c= 'r', linewidths=2.0, zorder=10000)
N2 = len(stations_ndbc)


## Tide Gauges
ax3 = fig.add_subplot(2,3,5)
map3 = plot_map_background(plot_area, ax = ax3, parallels_labels = [0,0,0,0])


## NOAA Tide Gauges
stations1 = pd.read_json('/home/orca/bkerns/projects/coastal/irene_lee_rainfall_paper/usgs_stream_flow/data.bak/processed/stations_noaa.json')
stations1.reset_index(inplace=True)

stations_with_navd88 = pd.read_csv('../datum_reference/noaa_tide_gauges_ssh_offset.csv')


H3 = ax3.scatter(stations1['lng'], stations1['lat'], marker='^', s=70, c= 'orange', linewidths=0.5, edgecolor='k')
N3 = len(stations1)

H33 = ax3.scatter(stations_with_navd88['lon'], stations_with_navd88['lat'], marker='^', s=70, c= 'orange', linewidths=1.5, edgecolor='k')


## USGS Estuary Sites
ax4 = fig.add_subplot(2,3,6)
map4 = plot_map_background(plot_area, ax = ax4, parallels_labels = [0,0,0,0])


## Estuary and tide gauge stations.
stations2 = pd.read_json('/home/orca/bkerns/projects/coastal/irene_lee_rainfall_paper/usgs_stream_flow/data.bak/processed/stations_estuary.json')
stations2.reset_index(inplace=True)
H4 = ax4.scatter(stations2['dec_long_va'], stations2['dec_lat_va'], s=50, c= 'c', linewidth=1.5, edgecolor='k')
N4 = len(stations2)

helpers.add_panel_label(ax1, 'a.')
helpers.add_panel_label(ax11, 'b.')
helpers.add_panel_label(ax2, 'c.')
helpers.add_panel_label(ax3, 'd.')
helpers.add_panel_label(ax4, 'e.')


## Legend
ax2.legend(handles=[H2,], labels=['NDBC Buoys',], bbox_to_anchor=[1.05, 0.95], loc='lower right', fontsize=9, frameon=False)
ax3.legend(handles=[H3,], labels=['NOAA Tide Gauges',], bbox_to_anchor=[1.05, 0.95], loc='lower right', fontsize=9, frameon=False)
ax4.legend(handles=[H4,], labels=['USGS Sites',], bbox_to_anchor=[1.05, 0.95], loc='lower right', fontsize=9, frameon=False)

## Formatting
plt.tight_layout()
fn_out_base = 'figure2_stations_map.revised'
#for ext in ['png','pdf']:
for ext in ['png',]:
    print(fn_out_base + '.' + ext)
    plt.savefig(fn_out_base + '.' + ext, dpi=150, bbox_inches='tight')


