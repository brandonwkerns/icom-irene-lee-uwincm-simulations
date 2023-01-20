import numpy as np
import pandas as pd
from netCDF4 import Dataset
import xarray as xr
import glob
import datetime as dt
import matplotlib.pyplot as plt

## Import Cartopy stuff.
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
import matplotlib.patheffects as path_effects
import matplotlib.ticker as mticker
import colorcet as cc
import cmocean

plt.close('all')

def plot_map_background(plot_area=[-180,180,-80,80], ax = None, crs=None):
    """
    Plot map background using the area specified as plot_area.
    NOTE: The projection MUST already have been defined when the axis was created!
    e.g., ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
    Projection information will be retrieved from the axis object using ax.projection.
    """

    # ax.set_extent(plot_area, ccrs.Geodetic())


    # crs = ccrs.Mercator()

    # Create a Stamen Terrain instance.
    if ax == None:
        stamen_terrain = cimgt.Stamen(style='terrain', cache=True)
        ax = plt.axes(projection=crs) #stamen_terrain.crs)
        ax.set_extent(plot_area)
    # ax.add_image(stamen_terrain, 8)

    res = '110m'  # Resolution/scale options: '10m','50m','110m' (higher scale --> higher spatial resolution)

    ## These pre-defined features are all at 110m resolution.
    # ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)
    # ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale=res,
        facecolor='none')
    ax.add_feature(states_provinces, linewidth=0.5)

    # Gridlines only work with ccrs.PlateCarree(). Plot still seems fine with axis using other projections.
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.5, color='k', alpha=0.7, linestyle='--')
    gl.top_labels = False
    # gl.left_labels = False
    # gl.bottom_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 2))
    gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, 2))

    return ax




def get_colormap():

    ## Colormaps
    CMAP={}
    CMAP['rainrate24h_global']={}
    CMAP['rainrate24h_global']['vmin'] = 0.0
    CMAP['rainrate24h_global']['vmax'] = 12.0
    CMAP['rainrate24h_global']['interval'] = 0.5

    cmap = plt.cm.jet
    N_keep = 12
    cmap2 = colors.LinearSegmentedColormap.from_list(name='custom', colors=cmap(range(cmap.N)), N=N_keep)
    color_list1 = cmap2(range(cmap2.N))
    ## Modify individual colors here.
    color_list11 = cmap2(range(cmap2.N))
    color_list11[0] = [1,1,1,1] # Make first color white.
    color_list11[1] = [0.7,0.7,0.7, 1] # Make first color gray.
    color_list11[2] = [0.0, 0.9, 1.0, 1.0]
    color_list11[3] = [0.0/255.0, 156.0/255.0, 255.0/255.0, 1.0]
    color_list11[4] = [0.0/255.0, 58.0/255.0, 255.0/255.0, 1.0]
    color_list11[5] = [0.3, 1.0, 0.3, 1.]
    color_list11[6] = [0.0, 0.8, 0.0, 1.]
    color_list11[9] = color_list1[10]  # Move top to colors of jet color map down.
    color_list11[10] = color_list1[11] # Move top to colors of jet color map down.
    color_list11[11] = [1,0,1,1] # Magenta on top.


    print(color_list11)
    ## Add color map to the dict.
    CMAP = colors.LinearSegmentedColormap.from_list(name='custom', colors=color_list11, N=N_keep)

    return CMAP


def set_up_plot(ax, plot_area=[-180,180,-80,80]):
    map = plot_map_background(plot_area, ax=ax, coastline_line_width = 1.0)
    return map


def get_atcf_storm_info(this_fn):

    cols = [2,6,7,8,9,19]
    names = ['yyyymmddhh','lat0','lon0','wspd','slp','rmw']
    this_data = pd.read_csv(this_fn, header=None, usecols=cols, names=names)
    this_data['datetime'] = [dt.datetime.strptime(x, '%Y%m%d%H')  for x in this_data['yyyymmddhh'].astype(str).values]
    this_data['timestamp'] = [(dt.datetime.strptime(x, '%Y%m%d%H') - dt.datetime(1970,1,1,0,0,0)).total_seconds()  for x in this_data['yyyymmddhh'].astype(str).values]

    this_data['lat'] = [np.double(x[:-1])/10.0 for x in this_data['lat0'].astype(str).values]
    this_data['lon'] = [np.double(x[:-1])/10.0 if x[-1] == 'E' else -1.0*np.double(x[:-1])/10.0 for x in this_data['lon0'].astype(str).values]

    ## Round.
    this_data['wspd'] = np.round(2.0*this_data['wspd'])/2.0
    this_data['slp'] = np.round(1.0*this_data['slp'])/1.0

    ## Get unique values. Use Pandas drop_duplicates()
    this_data.drop_duplicates(subset='yyyymmddhh', inplace=True)


    return this_data


def plot_best_track(map1, best_track, color='k'
        , hour_markers_large=[], hour_markers_small=[]):

    map1.plot(best_track['lon'].values, best_track['lat'].values, '-'
              , color='k', label='Best Track'
              , linewidth=3,transform=ccrs.PlateCarree())
    map1.plot(best_track['lon'].values, best_track['lat'].values, '--'
              , color=color, label='Best Track'
              , linewidth=1.5,transform=ccrs.PlateCarree())

    ## Hour Markers.
    for hh in hour_markers_small:
        hour_idx = [x for x in range(len(best_track['datetime'].values)) if int(pd.DatetimeIndex(best_track['datetime']).hour[x]) == int(hh)]

        map1.scatter(best_track['lon'].values[hour_idx], best_track['lat'].values[hour_idx]
                     , c='k', s=20,transform=ccrs.PlateCarree())

    for hh in hour_markers_large:
        hour_idx = [x for x in range(len(best_track['datetime'].values)) if int(pd.DatetimeIndex(best_track['datetime']).hour[x]) == int(hh)]

        map1.scatter(best_track['lon'].values[hour_idx], best_track['lat'].values[hour_idx]
                     , c='k', s=40,transform=ccrs.PlateCarree())


################################################################################
################################################################################
################################################################################

CMAP = "jet"
#CMAP_wspd = "jet" #cc.m_rainbow
CMAP_wspd = colors.ListedColormap(plt.cm.jet(range(0,255,15)))


#CMAP_swh = cc.m_linear_worb_100_25_c53

"""
#   30 colors from:        NCL  BlGrYeOrReVi200.rgb
colors_swh = [[0,0,255],        # line 3
              [1,31,213],     # 10
              [3,67,165],     # 18
              [6,139,69],    # 34
              [51,192,12],   # 50
              [150,223,6],   # 66
              [199,238,3],    # 74
              [249,253,0],    # 82
              [255,239,0],    # 90
              [255,220,0],    # 98
              [255 209   0].   #103
"""

colors_swh_all = np.loadtxt('BlGrYeOrReVi200.rgb', skiprows=2)
colors_swh_idx = [0,7,15,31,47,63,71,79,87,95,
                  100,105,110,115,120,125,130,135,140,145,
                  150,155,160,165,170,175,180,185,190,195]
colors_swh = [colors_swh_all[x] for x in colors_swh_idx]                  
colors_swh = np.array(colors_swh) / 255.0
CMAP_swh = colors.LinearSegmentedColormap.from_list(
    'Custom cmap', colors_swh, len(colors_swh[:,0]))

CMAP_ssh = 'viridis'

F={}
print('irene_and_lee__impact_fields.nc')
DS = xr.open_mfdataset('/home/orca/bkerns/projects/coastal/models/irene_and_lee_2011/impact_fields/relo_0826_v2_extend_irene/irene_and_lee_hourly_fields__*.nc')
DS_irene = DS.sel(time=slice('2011-08-24T00:00:00', '2011-09-02T00:00:00'))
DS_lee = DS.sel(time=slice('2011-09-02T00:00:00', '2011-09-12T00:00:00'))

X = DS['longitude_wrf'][0,:,:]
Y = DS['latitude_wrf'][0,:,:]

plot_area2 = [-80.0 + 360.0,-71.5 + 360.0,33,41.5]

## Map
fig=plt.figure(figsize=(7.5, 9.0))
# stamen_terrain = cimgt.Stamen(style='terrain-background', cache=True)

# crs = ccrs.Mercator()
# crs = stamen_terrain.crs
crs = ccrs.PlateCarree()

AX=[]
CAX=[]
# MAP=[]
for ii in range(6):
    ax = fig.add_subplot(3, 2, ii+1, projection=crs)
    # ax.add_image(stamen_terrain, 8)

    plot_map_background(plot_area2, ax=ax, crs=crs)
    ax.set_extent(plot_area2)
    # ax.add_image(stamen_terrain, 6)

    AX += [ax]

    # Add the colorbar axes anywhere in the figure. Its position will be
    # re-calculated at each figure resize. 
    CAX += [fig.add_axes([0, 0, 0.1, 0.1])]


fig.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)


def resize_colorbar(event):
    plt.draw()

    for nn, ax in enumerate(AX):
        posn = ax.get_position()
        CAX[nn].set_position([posn.x0 + posn.width + 0.01, posn.y0,
                            0.025, posn.height])

fig.canvas.mpl_connect('resize_event', resize_colorbar)

##
## Wind Speed
##

# Irene
plt.sca(AX[0])
Z = DS_irene['wspd'].max(dim='time')
cnt = AX[0].pcolormesh(X, Y, Z, vmin=0, vmax=34, cmap=CMAP_wspd, rasterized=True, transform=ccrs.PlateCarree())    
hcb = plt.colorbar(cnt, label='Max Wind Speed [m/s]', cax=CAX[0], extend='max')

AX[0].set_title('a. WSPD, Irene')

## Add model track.
model_track_irene = get_atcf_storm_info('/home/disk/orca/bkerns/projects/coastal/irene_lee_rainfall_paper/domains/awovispy/Irene_simulation_final.txt')
plot_best_track(AX[0], model_track_irene, color='w'
        , hour_markers_large=[], hour_markers_small=[]) # Model

# Lee
plt.sca(AX[1])
Z = DS_lee['wspd'].max(dim='time')
cnt = AX[1].pcolormesh(X, Y, Z, vmin=0, vmax=34, cmap=CMAP_wspd, rasterized=True, transform=ccrs.PlateCarree())
hcb = plt.colorbar(cnt, label='Max Wind Speed [m/s]', cax=CAX[1], extend='max')
AX[1].set_title('b. WSPD, Lee')

##
## Wave Height
##

X = DS['longitude_umwm'][0,:,:]
Y = DS['latitude_umwm'][0,:,:]

# Irene
plt.sca(AX[2])
Z = DS_irene['swh'].max(dim='time')
cnt = plt.contourf(X, Y, Z, np.arange(0,12.5,0.5), cmap=CMAP_swh, extend='max',transform=ccrs.PlateCarree())

# This is the fix for the white lines between contour levels
for c in cnt.collections:
    c.set_edgecolor("face")

hcb = plt.colorbar(cnt, label='Sig. Wave Hgt. [m]', cax=CAX[2])
AX[2].set_title('c. SWH, Irene')

## Add model track.
model_track_irene = get_atcf_storm_info('/home/disk/orca/bkerns/projects/coastal/irene_lee_rainfall_paper/domains/awovispy/Irene_simulation_final.txt')
plot_best_track(AX[2], model_track_irene, color='w'
        , hour_markers_large=[0], hour_markers_small=[6,12,18]) # Model

# Lee
plt.sca(AX[3])
Z = DS_lee['swh'].max(dim='time')
cnt = plt.contourf(X, Y, Z, np.arange(0,12.5,0.5), cmap=CMAP_swh, extend='max',transform=ccrs.PlateCarree())

# This is the fix for the white lines between contour levels
for c in cnt.collections:
    c.set_edgecolor("face")

hcb = plt.colorbar(cnt, label='Sig. Wave Hgt. [m]', cax=CAX[3])
AX[3].set_title('d. SWH, Lee')

##
## SST Change
##

DS_irene2 = DS.sel(time=slice('2011-08-24T00:00:00', '2011-08-29T00:00:00'))

X = DS['longitude_hycom'][0,:,:]
Y = DS['latitude_hycom'][0,:,:]

SST1 = DS_irene2['sst_hycom'][0,:,:]
SST2 = DS_irene2['sst_hycom'][-1,:,:]

dSST = SST2 - SST1


plt.sca(AX[4])
cnt = AX[4].pcolormesh(X, Y, SST2, vmin=20.0, vmax=29.0, cmap='cet_rainbow', rasterized=True, transform=ccrs.PlateCarree())
hcb = plt.colorbar(cnt, label='SST [$\degree$C]', cax=CAX[4], extend='both')
AX[4].set_title('e. Post-Storm SST, Irene')

## Add model track.
model_track_irene = get_atcf_storm_info('/home/disk/orca/bkerns/projects/coastal/irene_lee_rainfall_paper/domains/awovispy/Irene_simulation_final.txt')
plot_best_track(AX[4], model_track_irene, color='w'
        , hour_markers_large=[], hour_markers_small=[]) # Model

plt.sca(AX[5])
cnt = AX[5].pcolormesh(X, Y, dSST, vmin=-5, vmax=5, cmap=cc.cm['CET_D1A'], rasterized=True, transform=ccrs.PlateCarree())
hcb = plt.colorbar(cnt, label='$\Delta$SST [$\degree$C]', cax=CAX[5], extend='both')
AX[5].set_title('f. SST Change, Irene')

## Add model track.
model_track_irene = get_atcf_storm_info('/home/disk/orca/bkerns/projects/coastal/irene_lee_rainfall_paper/domains/awovispy/Irene_simulation_final.txt')
plot_best_track(AX[5], model_track_irene, color='w'
        , hour_markers_large=[], hour_markers_small=[]) # Model


# divider = make_axes_locatable(AX[3])
# cax     = divider.append_axes("right", size="5%", pad=0.10) 
# hcb = plt.colorbar(cax=cax, label='Max Sig. Wave Hgt. [m]')


for ax in AX:
    pe = [path_effects.Stroke(linewidth=2.5, foreground='white'),
                       path_effects.Normal()]
    ax.plot(-74.692, 38.460, '^', color='k', markersize=8, markerfacecolor='none', path_effects=pe, transform=ccrs.PlateCarree())
    ax.plot(-76.415, 38.556, '^', color='k', markersize=8, markerfacecolor='none', path_effects=pe, transform=ccrs.PlateCarree())

for ax in AX[4:6]:
    # ax.plot(-73.64, 38.708, 'x', color='k', markersize=8, markerfacecolor='none', path_effects=pe, transform=ccrs.PlateCarree())
    ax.plot(-74.122, 38.08, 'x', color='k', markersize=8, markerfacecolor='none', path_effects=pe, transform=ccrs.PlateCarree())


plt.tight_layout(w_pad = 4)

print('irene_and_lee__max_swaths2.revised.png')
resize_colorbar(None)
plt.savefig('irene_and_lee__max_swaths2.revised.png',dpi=300,bbox_inches='tight')

print('irene_and_lee__max_swaths2.revised.pdf')
plt.savefig('irene_and_lee__max_swaths2.revised.pdf',dpi=300,bbox_inches='tight')

print('irene_and_lee__max_swaths2.revised.eps')
plt.savefig('irene_and_lee__max_swaths2.revised.eps',dpi=300,bbox_inches='tight')

plt.close(fig)
