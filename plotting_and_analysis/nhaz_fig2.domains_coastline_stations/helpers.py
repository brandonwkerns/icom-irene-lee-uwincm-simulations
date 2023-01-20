import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.basemap import Basemap
import json


def add_panel_label(ax, legend_text):
    ax.text(0.02, 1.02, legend_text, transform=ax.transAxes)


def plot_map_background0(plot_area = [-180, 180, -70, 70], ax = plt.gca()
                        , projection = 'cyl', resolution = 'i'
                        , draw_meridians = np.arange(0,360,2), meridians_line_width = 0.25
                        , meridians_labels = [0,0,0,1], meridians_fontsize=10
                        , meridians_dashes = [4,2]
                        , draw_parallels = np.arange(-90,90,2), parallels_line_width = 0.25
                        , parallels_labels = [1,0,0,0], parallels_fontsize=10
                        , parallels_dashes = [4,2]
                        , coastline_line_width = 0.7
                        , countries_line_width = 0.3
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
    # If no coastline in the entire plot area, drawcoastlines() would crash.
    try:
        map1.drawcoastlines(linewidth = coastline_line_width)
    except:
        pass
    # map1.drawstates(linewidth = countries_line_width,color='k')
    map1.drawcountries(linewidth = countries_line_width)

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




def plot_map_background(plot_area = [-180, 180, -70, 70], ax = plt.gca()
                        , projection = 'cyl', resolution = 'f'
                        , draw_meridians = np.arange(0,360,2), meridians_line_width = 0.25
                        , meridians_labels = [0,0,0,1], meridians_fontsize=10
                        , meridians_dashes = [4,2]
                        , draw_parallels = np.arange(-90,90,2), parallels_line_width = 0.25
                        , parallels_labels = [1,0,0,0], parallels_fontsize=10
                        , parallels_dashes = [4,2]
                        , coastline_line_width = 0.7
                        , countries_line_width = 0.3
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
    # If no coastline in the entire plot area, drawcoastlines() would crash.
    try:
        map1.drawcoastlines(linewidth = coastline_line_width)
    except:
        pass
    map1.drawstates(linewidth = countries_line_width,color='k')
    map1.drawcountries(linewidth = countries_line_width)

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



def plot_umwm_domain(ax, linewidth=1.0, color='k'):
    waves = xr.open_dataset(('/home/orca3/bkerns/models/' +
    'uwincm_hurricane_irene_umwm_fixes_new_coast_test_umwm3/' + 
    'run_aug26awo2_d02narrow_nudge0.0006uv_depth_limiter_new_coeff_fix_tides_big_umwm_v2_extend/' +
    'output/umwmout.grid'))

    lon = waves['lon']
    lat = waves['lat']

    H_waves=ax.plot(lon[:,0], lat[:,0], '--', linewidth=linewidth, color=color)[0]
    ax.plot(lon[0,:], lat[0,:], '--', linewidth=linewidth, color=color)
    ax.plot(lon[:,-1], lat[:,-1], '--', linewidth=linewidth, color=color)
    ax.plot(lon[-1,:], lat[-1,:], '--', linewidth=linewidth, color=color)

    return H_waves


def plot_hycom_domain(ax, linewidth=1.0, color='b'):
    hycom = xr.open_dataset(('/home/orca3/bkerns/models/' +
    'uwincm_hurricane_irene_umwm_fixes_new_coast_test_umwm3/' + 
    'run_aug26awo2_d02narrow_nudge0.0006uv_depth_limiter_new_coeff_fix_tides_big_umwm_v2_extend/' +
    'regional.depth.nc'))

    lon = hycom['longitude']
    lat = hycom['latitude']

    H_ocean=ax.plot(lon[:,0], lat[:,0], '-', linewidth=linewidth, color=color)[0]
    ax.plot(lon[0,:], lat[0,:], '-', linewidth=linewidth, color=color)
    ax.plot(lon[:,-1], lat[:,-1], '-', linewidth=linewidth, color=color)
    ax.plot(lon[-1,:], lat[-1,:], '-', linewidth=linewidth, color=color)

    return H_ocean


def plot_wrf_domain(ax, dom_num, linewidth=1.0, color='k'):
    wrf = xr.open_dataset(('/home/orca3/bkerns/models/' +
    'uwincm_hurricane_irene_umwm_fixes_new_coast_test_umwm3/' + 
    'run_aug26awo2_d02narrow_nudge0.0006uv_depth_limiter_new_coeff_fix_tides_big_umwm_v2_extend/' +
    'wrfout_d0{}_2011-08-24_00:00:00'.format(dom_num)))

    lon = wrf['XLONG'][0,:,:]
    lat = wrf['XLAT'][0,:,:]

    H_wrf=ax.plot(lon[:,0], lat[:,0], '-', linewidth=linewidth, color=color)[0]
    ax.plot(lon[0,:], lat[0,:], '-', linewidth=linewidth, color=color)
    ax.plot(lon[:,-1], lat[:,-1], '-', linewidth=linewidth, color=color)
    ax.plot(lon[-1,:], lat[-1,:], '-', linewidth=linewidth, color=color)

    return H_wrf


def plot_wrf_domain2(ax, dom_num, linewidth=1.3, color='k'):
    wrf = xr.open_dataset(('/home/orca3/bkerns/models/' +
    'uwincm_hurricane_irene_umwm_fixes_new_coast_test_umwm3/' + 
    'run_aug26awo2_d02narrow_nudge0.0006uv_depth_limiter_new_coeff_fix_tides_big_umwm_v2_extend/' +
    'wrfout_d0{}_2011-08-24_00:00:00'.format(dom_num)))

    lon = wrf['XLONG'][0,:,:]
    lat = wrf['XLAT'][0,:,:]

    ## Create path around the domain for patch.
    # X = lon[:,0]
    # X = np.append(X, lon[-1,:])
    # X = np.append(X, lon[::-1,-1])
    # X = np.append(X, lon[0,::-1])

    # Y = lat[:,0]
    # Y = np.append(Y, lat[-1,:])
    # Y = np.append(Y, lat[::-1,-1])
    # Y = np.append(Y, lat[0,::-1])

    # XY = np.dstack((X, Y))[0]
    # print(XY.shape)
    # ax.add_patch(patches.Polygon(XY, closed=True, alpha=0.5, hatch='//'))

    H_wrf=ax.plot(lon[:,0], lat[:,0], '--', linewidth=linewidth, color=color)[0]
    ax.plot(lon[0,:], lat[0,:], '--', linewidth=linewidth, color=color)
    ax.plot(lon[:,-1], lat[:,-1], '--', linewidth=linewidth, color=color)
    ax.plot(lon[-1,:], lat[-1,:], '--', linewidth=linewidth, color=color)

    return H_wrf


def plot_wrf_domain_lee(ax, dom_num, linewidth=1.3, color='k'):
    wrf = xr.open_dataset(('/home/orca3/bkerns/models/' +
    'uwincm_ts_lee/run_AWO1.3/' +
    'wrfout_d0{}_2011-09-02_00:00:00'.format(dom_num)))

    lon = wrf['XLONG'][0,:,:]
    lat = wrf['XLAT'][0,:,:]

    ## Create path around the domain for patch.
    # X = lon[:,0]
    # X = np.append(X, lon[-1,:])
    # X = np.append(X, lon[::-1,-1])
    # X = np.append(X, lon[0,::-1])

    # Y = lat[:,0]
    # Y = np.append(Y, lat[-1,:])
    # Y = np.append(Y, lat[::-1,-1])
    # Y = np.append(Y, lat[0,::-1])

    # XY = np.dstack((X, Y))[0]
    # print(XY.shape)
    # ax.add_patch(patches.Polygon(XY, closed=True, alpha=0.5, hatch='\\'))

    H_wrf=ax.plot(lon[:,0], lat[:,0], '-.', linewidth=linewidth, color=color)[0]
    ax.plot(lon[0,:], lat[0,:], '-.', linewidth=linewidth, color=color)
    ax.plot(lon[:,-1], lat[:,-1], '-.', linewidth=linewidth, color=color)
    ax.plot(lon[-1,:], lat[-1,:], '-.', linewidth=linewidth, color=color)

    return H_wrf


def get_ibtracs_csv_storm_info(csv_file, basin, year, name):

    """
    Extract storm data from IBTRACS csv file into a Pandas data frame.
    IBTRACS data can be obtained at https://www.ncdc.noaa.gov/ibtracs/index.php
    """
    data = pd.read_csv(csv_file, header=0, index_col=0, skiprows=[1], keep_default_na=False, na_values=['',' '])
    data_storm = data[(data['BASIN'].astype(str) == basin)
                          & (data['SEASON'].astype(str) == str(year))
                          & (data['NAME'].astype(str) == name.upper())]

    ## Add datetime column.
    data_storm.loc[:,'datetime'] = [dt.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')  for x in data_storm['ISO_TIME'].values]
    data_storm.loc[:,'timestamp'] = [(dt.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') - dt.datetime(1970,1,1,0,0,0)).total_seconds()  for x in data_storm['ISO_TIME'].values]

    return data_storm


def plot_best_track(best_track, ax):

    best_track_color='k'
    best_track_linewidth=1
    hour_markers_small=[] #(6,12,18)
    hour_markers_large=(0,)

    ax.plot(best_track['LON'].values, best_track['LAT'].values, '-'
              , color=best_track_color, label='Best Track'
              , linewidth=best_track_linewidth, zorder=1000)

    fmt = '%Y-%m-%d %H:%M:%S'
    datetime_list = [dt.datetime.strptime(x,fmt) for x in best_track['ISO_TIME'].values]

    ## Hour Markers.
    for hh in hour_markers_small:
        hour_idx = [x for x in range(len(best_track['datetime'].values)) if int(pd.DatetimeIndex(best_track['datetime'
]).hour[x]) == int(hh)]

        ax.scatter(best_track['LON'].values[hour_idx], best_track['LAT'].values[hour_idx]
                     , c='white', edgecolors=best_track_color, s=5, zorder=1010)

    for hh in hour_markers_large:
        hour_idx = [x for x in range(len(best_track['datetime'].values)) if int(pd.DatetimeIndex(best_track['datetime']).hour[x]) == int(hh)]

        ax.scatter(best_track['LON'].values[hour_idx], best_track['LAT'].values[hour_idx]
                     , c=best_track_color, s=10, zorder=1020)

        for this_hour_idx in hour_idx:
            DD = str(datetime_list[this_hour_idx].day)
            if int(DD) % 2 == 0:
                ax.text(best_track['LON'].values[this_hour_idx] + 1.0, best_track['LAT'].values[this_hour_idx], str(DD), fontsize=9, clip_on=True)


