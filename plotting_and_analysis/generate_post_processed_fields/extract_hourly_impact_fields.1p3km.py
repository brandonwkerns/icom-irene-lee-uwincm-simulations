import numpy as np
from netCDF4 import Dataset
import glob
import datetime as dt

import pyhycom as hy
import os

## Stitch the model rainfall of Hurricane Irene and Tropical Storm Lee (2011)
## and subset the area to the eastern seaboard.
## 1.3 km rain will be used when d03 is present (e.g., for Irene).
## 4 km rain will fill in 9 1.3 km points surrounding the 4 km point otherwise.

irene_dir = ('/home/orca3/bkerns/models/uwincm_hurricane_irene_umwm_fixes_new_coast_test_umwm3'
            + '/run_aug26awo2_d02narrow_nudge0.0006uv_depth_limiter_new_coeff_fix_tides_big_umwm_v2_extend/')

#pre_lee_dir = ('/home/orca3/bkerns/models/uwincm_ts_lee/run_AWO4_pre_lee_v2')
lee_dir = ('/home/orca3/bkerns/models/uwincm_ts_lee/run_AWO4_v2')


BOX = [-100,-55,18,48]

#### Functions.
def get_nearest_ij(grid_lon, grid_lat, x, y):
    dist = np.sqrt(np.power(grid_lon-x,2) + np.power(grid_lat-y,2))
    return np.unravel_index(np.argmin(dist), (grid_lon.shape))

def get_time_index(datetime_list, this_datetime):
    try:
        return [x for x in range(len(datetime_list)) if datetime_list[x] == this_datetime][0]
    except:
        return None

def refine_grid(lon_orig,lat_orig,refine_factor=3):

    ny,nx = lon_orig.shape

    ## Initialize
    lon_refined = np.zeros([refine_factor*(ny-1)+1,refine_factor*(nx-1)+1])
    lat_refined = np.zeros([refine_factor*(ny-1)+1,refine_factor*(nx-1)+1])

    ## Fill in the points for the rows I have.
    for jj in range(ny):
        jj_refined = refine_factor*(jj)
        X1 = np.arange(nx)
        X2 = np.arange(refine_factor * (nx-1)+1) / (1.0*refine_factor)

        ## Lon
        Y1 = lon_orig[jj,:]
        Y2 = np.interp(X2,X1,Y1)
        lon_refined[jj_refined,:] = Y2

        ## Lat
        Y1 = lat_orig[jj,:]
        Y2 = np.interp(X2,X1,Y1)
        lat_refined[jj_refined,:] = Y2

    ## Fill in the rows inbetween.
    for ii_refined in range(refine_factor*(nx-1)+1):
        X1 = np.arange(ny)
        X2 = np.arange(refine_factor * (ny-1)+1) / (1.0*refine_factor)
    
        ## Lon
        Y1 = lon_refined[::refine_factor,ii_refined]
        Y2 = np.interp(X2,X1,Y1)
        lon_refined[:,ii_refined] = Y2

        ## Lat
        Y1 = lat_refined[::refine_factor,ii_refined]
        Y2 = np.interp(X2,X1,Y1)
        lat_refined[:,ii_refined] = Y2
   
    return (lon_refined, lat_refined)
    

def refine_field(data_orig,refine_factor=3):

    ny,nx = data_orig.shape

    ## Initialize
    data_refined = np.zeros([refine_factor*(ny-1)+1,refine_factor*(nx-1)+1])

    ## Fill in the points for the rows I have.
    for jj in range(ny):
        jj_refined = refine_factor*(jj)
        X1 = np.arange(nx)
        X2 = np.arange(refine_factor * (nx-1)+1) / (1.0*refine_factor)

        ## Data
        Y1 = data_orig[jj,:]
        Y2 = np.interp(X2,X1,Y1)
        data_refined[jj_refined,:] = Y2


    ## Fill in the rows inbetween.
    for ii_refined in range(refine_factor*(nx-1)+1):
        X1 = np.arange(ny)
        X2 = np.arange(refine_factor * (ny-1)+1) / (1.0*refine_factor)
    
        ## Data
        Y1 = data_refined[::refine_factor,ii_refined]
        Y2 = np.interp(X2,X1,Y1)
        data_refined[:,ii_refined] = Y2
   
    return data_refined



def get_wind_direction(u,v):
    ## u and v are numpy arrays. Should be the same length.

    phi_radians0 = np.arctan2(v,u) ## CCW from East, oceanography convention
    phi_deg0 = phi_radians0 * 180.0 / 3.14159
    phi_deg = 90 - phi_deg0   ## CW from North, oceanography convention

    ## Get in range 0-360.
    phi_deg[phi_deg < 0.0] = phi_deg[phi_deg < 0.0] + 360.0
    phi_deg[phi_deg > 360.0] = phi_deg[phi_deg > 360.0] - 360.0

    return phi_deg

    
    
#### 1. Set up.


##
## Irene 12 km grid information.
##
fn = sorted(glob.glob(irene_dir+'/wrfout_d01*'))[0]
DS = Dataset(fn)
lon12km = DS['XLONG'][:][0]
lat12km = DS['XLAT'][:][0]
DS.close()

## Refine to 1.3 km.
lon4km_entire_domain, lat4km_entire_domain = refine_grid(lon12km, lat12km)
lon1p3km_entire_domain, lat1p3km_entire_domain = refine_grid(lon4km_entire_domain, lat4km_entire_domain)



##
## Use this to get the master 1.3 km grid.
##

##
## Extract box.
##

jj0, ii0 = get_nearest_ij(lon1p3km_entire_domain,lat1p3km_entire_domain, BOX[0], BOX[2])
jj0, ii_dum = get_nearest_ij(lon1p3km_entire_domain,lat1p3km_entire_domain, 0.5*(BOX[0]+BOX[1]), BOX[2]) # Accounting for the projection of the WRF domain.
jj1, ii1 = get_nearest_ij(lon1p3km_entire_domain,lat1p3km_entire_domain, BOX[1], BOX[3])
jjdum, ii1 = get_nearest_ij(lon1p3km_entire_domain,lat1p3km_entire_domain, BOX[1], BOX[2])

print('Subset area: ',[ii0,ii1,jj0,jj1])
ny_collect = jj1-jj0+1
nx_collect = ii1-ii0+1

lon1p3km_entire_domain_keep = lon1p3km_entire_domain[jj0:jj1+1,ii0:ii1+1]
lat1p3km_entire_domain_keep = lat1p3km_entire_domain[jj0:jj1+1,ii0:ii1+1]

##
## Now set up times.
##
datetime_begin = dt.datetime(2011,8,24,0,0,0)
datetime_end = dt.datetime(2011,9,12,0,0,0)
time_res_hours = 1

total_hours = int((datetime_end - datetime_begin).total_seconds()/3600.0)
datetime_list = [datetime_begin + dt.timedelta(hours=x) for x in range(total_hours+1)]

##
## Initialize the big rain rate array.
##
rain_collect = np.nan*np.zeros((len(datetime_list),ny_collect,nx_collect),dtype=np.float32)
psfc_collect = np.nan*np.zeros((len(datetime_list),ny_collect,nx_collect),dtype=np.float32)
olr_collect = np.nan*np.zeros((len(datetime_list),ny_collect,nx_collect),dtype=np.float32)
u_collect = np.nan*np.zeros((len(datetime_list),ny_collect,nx_collect),dtype=np.float32)
v_collect = np.nan*np.zeros((len(datetime_list),ny_collect,nx_collect),dtype=np.float32)
wspd_collect = np.nan*np.zeros((len(datetime_list),ny_collect,nx_collect),dtype=np.float32)
wdir_collect = np.nan*np.zeros((len(datetime_list),ny_collect,nx_collect),dtype=np.float32)
swh = np.nan*np.zeros((len(datetime_list),ny_collect,nx_collect),dtype=np.float32)
mss = np.nan*np.zeros((len(datetime_list),ny_collect,nx_collect),dtype=np.float32)

waves_seamask = np.nan*np.zeros((ny_collect,nx_collect),dtype=np.float32)


##
## Irene 
##
tidx_list = []
tt = -1
tt_last = 0
for fn in sorted(glob.glob(irene_dir+'/wrfout_d01*')):
    if '.orig' in fn:
        continue
    tt += 1
    print('{0:03d}: '.format(tt)+fn)
    fn0 = fn.split('/')[-1]
    fmt = 'wrfout_d01_%Y-%m-%d_%H:00:00'
    this_datetime = dt.datetime.strptime(fn0,fmt)
    tidx = get_time_index(datetime_list, this_datetime)
    print('tidx = ',tidx)
    tidx_list += [tidx]
    
    with Dataset(fn) as DS:
        lon12km_this = DS['XLONG'][:][0]
        lat12km_this = DS['XLAT'][:][0]
        psfc_this    = DS['PSFC'][:][0]
        olr_this    = DS['OLR'][:][0]
        rainc_this  = DS['RAINC'][:][0]
        rainnc_this = DS['RAINNC'][:][0]
        u10_this = DS['U10'][:][0]
        v10_this = DS['V10'][:][0]
        terrain_this     = DS['HGT'][:][0]
        landmask_this    = DS['XLAND'][:][0]  

    rain_this = rainc_this + rainnc_this
    
    psfc_interp0 = refine_field(psfc_this)
    psfc_interp = refine_field(psfc_interp0)
    olr_interp0 = refine_field(olr_this)
    olr_interp = refine_field(olr_interp0)
    rain_interp0 = refine_field(rain_this)
    rain_interp = refine_field(rain_interp0)
    u10_interp0 = refine_field(u10_this)
    u10_interp = refine_field(u10_interp0)
    v10_interp0 = refine_field(v10_this)
    v10_interp = refine_field(v10_interp0)

    ## 12 km terrain and landmask.
    terrain0     = refine_field(terrain_this)
    terrain      = refine_field(terrain0)
    landmask0    = refine_field(landmask_this)
    landmask     = refine_field(landmask0)

    
    ## Refine with d02.      
    fn2 = fn.replace('d01','d02')
    print('{0:03d}: '.format(tt)+fn2)
    fn0 = fn2.split('/')[-1]
    fmt = 'wrfout_d02_%Y-%m-%d_%H:00:00'
    this_datetime = dt.datetime.strptime(fn0,fmt)
    with Dataset(fn2) as DS:
        lon4km_this = DS['XLONG'][:][0]
        lat4km_this = DS['XLAT'][:][0]
        psfc_this    = DS['PSFC'][:][0]
        olr_this    = DS['OLR'][:][0]
        rainc_this  = DS['RAINC'][:][0]
        rainnc_this = DS['RAINNC'][:][0]
        u10_this = DS['U10'][:][0]
        v10_this = DS['V10'][:][0]
    ny4km, nx4km = lon4km_this.shape
    
    idx1d_lower_left = np.argmin(np.sqrt(np.power(lon4km_this[0,0]-lon1p3km_entire_domain,2)+np.power(lat4km_this[0,0]-lat1p3km_entire_domain,2)).flatten())
    jj_ll4, ii_ll4 = np.unravel_index(idx1d_lower_left, lon1p3km_entire_domain.shape)
    print(jj_ll4,ii_ll4)
    jj_ur4 = jj_ll4 + 3 * (ny4km-1) + 1
    ii_ur4 = ii_ll4 + 3 * (nx4km-1) + 1
    psfc_interp[jj_ll4:jj_ur4, ii_ll4:ii_ur4] = refine_field(psfc_this)
    olr_interp[jj_ll4:jj_ur4, ii_ll4:ii_ur4] = refine_field(olr_this)
    rain_interp[jj_ll4:jj_ur4, ii_ll4:ii_ur4] = refine_field(rainc_this + rainnc_this)
    u10_interp[jj_ll4:jj_ur4, ii_ll4:ii_ur4] = refine_field(u10_this)
    v10_interp[jj_ll4:jj_ur4, ii_ll4:ii_ur4] = refine_field(v10_this)


    ## Refine with d03.
    fn3 = fn.replace('d01','d03')
    if os.path.exists(fn3):
        print('{0:03d}: '.format(tt)+fn3)
        fn0 = fn3.split('/')[-1]
        fmt = 'wrfout_d03_%Y-%m-%d_%H:00:00'
        this_datetime = dt.datetime.strptime(fn0,fmt)
        with Dataset(fn3) as DS:
            lon1p3km_this = DS['XLONG'][:][0]
            lat1p3km_this = DS['XLAT'][:][0]
            psfc_this    = DS['PSFC'][:][0]
            olr_this    = DS['OLR'][:][0]
            rainc_this  = DS['RAINC'][:][0]
            rainnc_this = DS['RAINNC'][:][0]
            u10_this = DS['U10'][:][0]
            v10_this = DS['V10'][:][0]
            
        ny1p3km, nx1p3km = lon1p3km_this.shape
        rain_1p3km_this = rainc_this + rainnc_this
        
        
        jstart_d03, istart_d03 = get_nearest_ij(lon1p3km_this,lat1p3km_this,lon1p3km_entire_domain[0,0],lat1p3km_entire_domain[0,0])        
        if jstart_d03 < ny1p3km-1: # If entire d03 is outside of the cut area, skip this step.

            jj_ll, ii_ll = get_nearest_ij(lon1p3km_entire_domain,lat1p3km_entire_domain,lon1p3km_this[0,0],lat1p3km_this[0,0])
            print(jj_ll,ii_ll)


            jj_ur = jj_ll + ny1p3km - jstart_d03
            ii_ur = ii_ll + nx1p3km - istart_d03
            psfc_interp[jj_ll:jj_ur, ii_ll:ii_ur] = psfc_this[jstart_d03:ny1p3km,istart_d03:nx1p3km]
            olr_interp[jj_ll:jj_ur, ii_ll:ii_ur] = olr_this[jstart_d03:ny1p3km,istart_d03:nx1p3km]
            #rain_interp[jj_ll+1:jj_ur-1, ii_ll+1:ii_ur-1] = rain_1p3km_this[jstart_d03+1:ny1p3km-1,istart_d03+1:nx1p3km-1] # For rain accumulation, don't use the edges to avoid moving nest artifacts.
            rain_interp[jj_ll:jj_ur, ii_ll:ii_ur] = rain_1p3km_this[jstart_d03:ny1p3km,istart_d03:nx1p3km]
            u10_interp[jj_ll:jj_ur, ii_ll:ii_ur] = u10_this[jstart_d03:ny1p3km,istart_d03:nx1p3km]
            v10_interp[jj_ll:jj_ur, ii_ll:ii_ur] = v10_this[jstart_d03:ny1p3km,istart_d03:nx1p3km]

        else:
            print('d03 completely outside of cut area.')

        #import matplotlib.pyplot as plt
        #plt.figure() ; plt.pcolormesh(rain_interp, cmap='jet') ; plt.show()

            
    wspd10_interp = np.sqrt(np.power(u10_interp,2) + np.power(v10_interp,2))
    wdir10_interp = get_wind_direction(u10_interp,v10_interp)

    rain_collect[tidx,:,:] = rain_interp[jj0:jj1+1,ii0:ii1+1]
    psfc_collect[tidx,:,:] = psfc_interp[jj0:jj1+1,ii0:ii1+1]
    olr_collect[tidx,:,:] = olr_interp[jj0:jj1+1,ii0:ii1+1]
    u_collect[tidx,:,:] = u10_interp[jj0:jj1+1,ii0:ii1+1]
    v_collect[tidx,:,:] = v10_interp[jj0:jj1+1,ii0:ii1+1]
    wspd_collect[tidx,:,:] = wspd10_interp[jj0:jj1+1,ii0:ii1+1]
    wdir_collect[tidx,:,:] = wdir10_interp[jj0:jj1+1,ii0:ii1+1]


    
    ##
    ## Waves stuff.
    ##
    fmt_waves = 'umwmout_%Y-%m-%d_%H:00:00.nc'
    fn_waves = (irene_dir+'/output/' + this_datetime.strftime(fmt_waves))
    DS = Dataset(fn_waves)
    lon_waves_this = DS['lon'][:][0]
    lat_waves_this = DS['lat'][:][0]
    swh_this = DS['swh'][:][0]
    mss_this = DS['mss'][:][0]
    waves_seamask_this = DS['seamask'][:][0]
    DS.close()
    ny_umwm, nx_umwm = lon_waves_this.shape

    ## Lower left corner of wave model is not necessarily in the cut out area.
    jstart, istart = get_nearest_ij(lon_waves_this,lat_waves_this,BOX[0],BOX[2])
    jend, iend = get_nearest_ij(lon_waves_this,lat_waves_this,BOX[1],BOX[3])

    if tt == 0:
        lon_waves_keep = lon_waves_this[jstart:jend+1,istart:iend+1]
        lat_waves_keep = lat_waves_this[jstart:jend+1,istart:iend+1]
        ny_waves ,nx_waves = lon_waves_keep.shape

        waves_seamask_keep = waves_seamask[jstart:jend+1,istart:iend+1]
        
        swh = np.zeros((len(datetime_list),ny_waves,nx_waves))
        mss = np.zeros((len(datetime_list),ny_waves,nx_waves))

    swh[tidx,:,:] = swh_this[jstart:jend+1,istart:iend+1]
    mss[tidx,:,:] = mss_this[jstart:jend+1,istart:iend+1]


    ## Ocean stuff.
    fmt_hy = 'archv.%Y_%j_%H.a'
    fn_hy = (irene_dir+'/' + this_datetime.strftime(fmt_hy))
    lon_hy = hy.getField('plon',irene_dir+'/regional.grid.a')
    lat_hy = hy.getField('plat',irene_dir+'/regional.grid.a')
    t_hy = hy.getField('temp', fn_hy, layers=[0,])
    s_hy = hy.getField('salin', fn_hy, layers=[0,])
    u_hy = hy.getField('u-vel', fn_hy, layers=[0,])
    v_hy = hy.getField('v-vel', fn_hy, layers=[0,])
    ub_hy = hy.getField('u_btrop', fn_hy)
    vb_hy = hy.getField('v_btrop', fn_hy)
    ssh_hy = hy.getField('srfhgt', fn_hy)
    print(lon_hy.shape)

    jstart, istart = get_nearest_ij(lon_hy,lat_hy,BOX[0],BOX[2])
    jend, iend = get_nearest_ij(lon_hy,lat_hy,BOX[1],BOX[3])
    if tt == 0:
        lon_hy_keep = lon_hy[jstart:jend+1,istart:iend+1]
        lat_hy_keep = lat_hy[jstart:jend+1,istart:iend+1]
        ny_hy ,nx_hy = lon_hy_keep.shape

        bathy = hy.getBathymetry(irene_dir+'/regional.depth.a')
        bathy_hy_keep = bathy[jstart:jend+1,istart:iend+1]
        
        sst_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        sss_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        u_sfc_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        v_sfc_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        u_baro_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        v_baro_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        ssh_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        
    sst_hy_keep[tidx,:,:] = t_hy[0,jstart:jend+1,istart:iend+1]
    sss_hy_keep[tidx,:,:] = s_hy[0,jstart:jend+1,istart:iend+1]
    u_sfc_hy_keep[tidx,:,:] = u_hy[0,jstart:jend+1,istart:iend+1]
    v_sfc_hy_keep[tidx,:,:] = v_hy[0,jstart:jend+1,istart:iend+1]
    u_baro_hy_keep[tidx,:,:] = ub_hy[jstart:jend+1,istart:iend+1]
    v_baro_hy_keep[tidx,:,:] = vb_hy[jstart:jend+1,istart:iend+1]
    ssh_hy_keep[tidx,:,:] = ssh_hy[jstart:jend+1,istart:iend+1]

    
## These are accumulated rain values. So subtract to get hourly rain.
rain_collect[0,:,:] = rain_collect[1,:,:]  # t=0 rain is zero, but set it to the first hour rain for continuity.
olr_collect[0,:,:]  = olr_collect[1,:,:]   # t=0 olr is zero, but set it to the first hour rain for continuity.

for tt0 in range(tt,tt_last+1,-1):
    if tt0 == 49:
        ## Vortex relocation at 26th 00 UTC. Do NOT subtract previous value for 26th 01 UTC.
        continue
    else:
        rain_collect[tt0,:,:] = rain_collect[tt0,:,:] - rain_collect[tt0-1,:,:]



##
## For Lee.
##

tt_last = tt+1
tt_last_lee = tt+1


for fn in sorted(glob.glob(lee_dir+'/wrfout_d01*'))[1:]: # I skip the init. time.
    tt += 1
    print('{0:03d}: '.format(tt)+fn)
    fn0 = fn.split('/')[-1]
    fmt = 'wrfout_d01_%Y-%m-%d_%H:00:00'
    this_datetime = dt.datetime.strptime(fn0,fmt)
    tidx = get_time_index(datetime_list, this_datetime)
    print('tidx = ',tidx)
    tidx_list += [tidx]
    
    DS = Dataset(fn)
    lon12km_this = DS['XLONG'][:][0]
    lat12km_this = DS['XLAT'][:][0]
    psfc_this    = DS['PSFC'][:][0]
    olr_this    = DS['OLR'][:][0]
    rainc_this  = DS['RAINC'][:][0]
    rainnc_this = DS['RAINNC'][:][0]
    u10_this = DS['U10'][:][0]
    v10_this = DS['V10'][:][0]
    terrain_this     = DS['HGT'][:][0]
    landmask_this    = DS['XLAND'][:][0]  
    DS.close()

    rain_this = rainc_this + rainnc_this
    
    psfc_interp0 = refine_field(psfc_this)
    psfc_interp = refine_field(psfc_interp0)
    olr_interp0 = refine_field(olr_this)
    olr_interp = refine_field(olr_interp0)
    rain_interp0 = refine_field(rain_this)
    rain_interp = refine_field(rain_interp0)
    u10_interp0 = refine_field(u10_this)
    u10_interp = refine_field(u10_interp0)
    v10_interp0 = refine_field(v10_this)
    v10_interp = refine_field(v10_interp0)

    ## 12 km terrain and landmask.
    terrain0     = refine_field(terrain_this)
    terrain      = refine_field(terrain0)
    landmask0    = refine_field(landmask_this)
    landmask     = refine_field(landmask0)

    
    ## Refine with d02.      
    fn2 = fn.replace('d01','d02')
    print('{0:03d}: '.format(tt)+fn2)
    fn0 = fn2.split('/')[-1]
    fmt = 'wrfout_d02_%Y-%m-%d_%H:00:00'
    this_datetime = dt.datetime.strptime(fn0,fmt)
    DS = Dataset(fn2)
    lon4km_this = DS['XLONG'][:][0]
    lat4km_this = DS['XLAT'][:][0]
    psfc_this    = DS['PSFC'][:][0]
    olr_this    = DS['OLR'][:][0]
    rainc_this  = DS['RAINC'][:][0]
    rainnc_this = DS['RAINNC'][:][0]
    u10_this = DS['U10'][:][0]
    v10_this = DS['V10'][:][0]
    DS.close()
    ny4km, nx4km = lon4km_this.shape
    rain_this = rainc_this + rainnc_this
    
    idx1d_lower_left = np.argmin(np.sqrt(np.power(lon4km_this[0,0]-lon1p3km_entire_domain,2)+np.power(lat4km_this[0,0]-lat1p3km_entire_domain,2)).flatten())
    jj_ll, ii_ll = np.unravel_index(idx1d_lower_left, lon1p3km_entire_domain.shape)
    print(jj_ll,ii_ll)
    jj_ur = jj_ll + 3 * (ny4km-1) + 1
    ii_ur = ii_ll + 3 * (nx4km-1) + 1
    psfc_interp[jj_ll:jj_ur, ii_ll:ii_ur] = refine_field(psfc_this)
    olr_interp[jj_ll:jj_ur, ii_ll:ii_ur] = refine_field(olr_this)
    rain_interp[jj_ll:jj_ur, ii_ll:ii_ur] = refine_field(rain_this)
    u10_interp[jj_ll:jj_ur, ii_ll:ii_ur] = refine_field(u10_this)
    v10_interp[jj_ll:jj_ur, ii_ll:ii_ur] = refine_field(v10_this)



    wspd10_interp = np.sqrt(np.power(u10_interp,2) + np.power(v10_interp,2))
    wdir10_interp = get_wind_direction(u10_interp,v10_interp)

    rain_collect[tidx,:,:] = rain_interp[jj0:jj1+1,ii0:ii1+1]
    psfc_collect[tidx,:,:] = psfc_interp[jj0:jj1+1,ii0:ii1+1]
    olr_collect[tidx,:,:] = olr_interp[jj0:jj1+1,ii0:ii1+1]
    u_collect[tidx,:,:] = u10_interp[jj0:jj1+1,ii0:ii1+1]
    v_collect[tidx,:,:] = v10_interp[jj0:jj1+1,ii0:ii1+1]
    wspd_collect[tidx,:,:] = wspd10_interp[jj0:jj1+1,ii0:ii1+1]
    wdir_collect[tidx,:,:] = wdir10_interp[jj0:jj1+1,ii0:ii1+1]


    
    ##
    ## Waves stuff.
    ##
    fmt_waves = 'umwmout_%Y-%m-%d_%H:00:00.nc'
    fn_waves = (lee_dir+'/output/' + this_datetime.strftime(fmt_waves))
    DS = Dataset(fn_waves)
    lon_waves_this = DS['lon'][:][0]
    lat_waves_this = DS['lat'][:][0]
    swh_this = DS['swh'][:][0]
    mss_this = DS['mss'][:][0]
    waves_seamask_this = DS['seamask'][:][0]
    DS.close()
    ny_umwm, nx_umwm = lon_waves_this.shape

    jstart, istart = get_nearest_ij(lon_waves_this,lat_waves_this,BOX[0],BOX[2])
    jend, iend = get_nearest_ij(lon_waves_this,lat_waves_this,BOX[1],BOX[3])

    if tt == 0:
        lon_waves_keep = lon_waves_this[jstart:jend+1,istart:iend+1]
        lat_waves_keep = lat_waves_this[jstart:jend+1,istart:iend+1]
        ny_waves ,nx_waves = lon_waves_keep.shape

        waves_seamask_keep = waves_seamask[jstart:jend+1,istart:iend+1]
        
        swh = np.zeros((len(datetime_list),ny_waves,nx_waves))
        mss = np.zeros((len(datetime_list),ny_waves,nx_waves))

    swh[tidx,:,:] = swh_this[jstart:jend+1,istart:iend+1]
    mss[tidx,:,:] = mss_this[jstart:jend+1,istart:iend+1]


    ## Ocean stuff.
    fmt_hy = 'archv.%Y_%j_%H.a'
    fn_hy = (lee_dir+'/' + this_datetime.strftime(fmt_hy))
    lon_hy = hy.getField('plon',irene_dir+'/regional.grid.a')
    lat_hy = hy.getField('plat',irene_dir+'/regional.grid.a')
    t_hy = hy.getField('temp', fn_hy, layers=[0,])
    s_hy = hy.getField('salin', fn_hy, layers=[0,])
    u_hy = hy.getField('u-vel', fn_hy, layers=[0,])
    v_hy = hy.getField('v-vel', fn_hy, layers=[0,])
    ub_hy = hy.getField('u_btrop', fn_hy)
    vb_hy = hy.getField('v_btrop', fn_hy)
    ssh_hy = hy.getField('srfhgt', fn_hy)
    print(lon_hy.shape)

    jstart, istart = get_nearest_ij(lon_hy,lat_hy,BOX[0],BOX[2])
    jend, iend = get_nearest_ij(lon_hy,lat_hy,BOX[1],BOX[3])
    if tt == 0:
        lon_hy_keep = lon_hy[jstart:jend+1,istart:iend+1]
        lat_hy_keep = lat_hy[jstart:jend+1,istart:iend+1]
        ny_hy ,nx_hy = lon_hy_keep.shape

        bathy = hy.getBathymetry(irene_dir+'/regional.depth.a')
        bathy_hy_keep = bathy[jstart:jend+1,istart:iend+1]
        
        sst_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        sss_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        u_sfc_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        v_sfc_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        u_baro_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        v_baro_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        ssh_hy_keep = np.zeros((len(datetime_list),ny_hy,nx_hy))
        
    sst_hy_keep[tidx,:,:] = t_hy[0,jstart:jend+1,istart:iend+1]
    sss_hy_keep[tidx,:,:] = s_hy[0,jstart:jend+1,istart:iend+1]
    u_sfc_hy_keep[tidx,:,:] = u_hy[0,jstart:jend+1,istart:iend+1]
    v_sfc_hy_keep[tidx,:,:] = v_hy[0,jstart:jend+1,istart:iend+1]
    u_baro_hy_keep[tidx,:,:] = ub_hy[jstart:jend+1,istart:iend+1]
    v_baro_hy_keep[tidx,:,:] = vb_hy[jstart:jend+1,istart:iend+1]
    ssh_hy_keep[tidx,:,:] = ssh_hy[jstart:jend+1,istart:iend+1]

    
## These are accumulated rain values. So subtract to get hourly rain.
for tt0 in range(tt,tt_last+1,-1):
    rain_collect[tt0,:,:] = rain_collect[tt0,:,:] - rain_collect[tt0-1,:,:]


terrain  = terrain[jj0:jj1+1,ii0:ii1+1]
landmask = landmask[jj0:jj1+1,ii0:ii1+1]

## Eliminate spurious negative rainfall caused by moving nest artifact.
rain_collect[rain_collect < 0.0] = 0.0
    

#### 3. Save hourly NetCDF files.
for tt0 in tidx_list:

    fn_out = 'irene_and_lee_hourly_fields__'+datetime_list[tt0].strftime('%Y%m%d%H')+'.nc'

    print('Writing: '+fn_out)
    with Dataset(fn_out,'w') as DS:

        DS.createDimension('time',1)
        DS.createDimension('y',ny_collect)
        DS.createDimension('x',nx_collect)
        DS.createDimension('y_umwm',ny_waves)
        DS.createDimension('x_umwm',nx_waves)
        DS.createDimension('y_hy',ny_hy)
        DS.createDimension('x_hy',nx_hy)
        DS.createVariable('time','f',('time',))
        DS.createVariable('longitude_wrf','f',('time','y','x',))
        DS.createVariable('latitude_wrf','f',('time','y','x',))
        DS.createVariable('terrain_wrf','f',('time','y','x',),zlib=True)
        DS.createVariable('landmask_wrf','f',('time','y','x',),zlib=True)
        DS.createVariable('rain_hourly','f',('time','y','x',),zlib=True)
        DS.createVariable('psfc','f',('time','y','x',),zlib=True)
        DS.createVariable('olr','f',('time','y','x',),zlib=True)
        DS.createVariable('u','f',('time','y','x',),zlib=True)
        DS.createVariable('v','f',('time','y','x',),zlib=True)
        DS.createVariable('wspd','f',('time','y','x',),zlib=True)
        DS.createVariable('wdir','f',('time','y','x',),zlib=True)

        DS.createVariable('longitude_umwm','f',('time','y_umwm','x_umwm',))
        DS.createVariable('latitude_umwm','f',('time','y_umwm','x_umwm',))
        DS.createVariable('waves_seamask','f',('time','y_umwm','x_umwm',),zlib=True)
        DS.createVariable('swh','f',('time','y_umwm','x_umwm',),zlib=True)
        DS.createVariable('mss','f',('time','y_umwm','x_umwm',),zlib=True)
        
        DS.createVariable('longitude_hycom','f',('time','y_hy','x_hy',))
        DS.createVariable('latitude_hycom','f',('time','y_hy','x_hy',))
        DS.createVariable('bathymetry_hycom','f',('time','y_hy','x_hy',))
        DS.createVariable('sst_hycom','f',('time','y_hy','x_hy',),zlib=True)
        DS.createVariable('sss_hycom','f',('time','y_hy','x_hy',),zlib=True)
        
        DS.createVariable('u_sfc_hycom','f',('time','y_hy','x_hy',),zlib=True)
        DS.createVariable('u_baro_hycom','f',('time','y_hy','x_hy',),zlib=True)
        DS.createVariable('v_sfc_hycom','f',('time','y_hy','x_hy',),zlib=True)
        DS.createVariable('v_baro_hycom','f',('time','y_hy','x_hy',),zlib=True)
        DS.createVariable('spd_sfc_hycom','f',('time','y_hy','x_hy',),zlib=True)
        DS.createVariable('spd_baro_hycom','f',('time','y_hy','x_hy',),zlib=True)
        DS.createVariable('spd_sfc_plus_baro_hycom','f',('time','y_hy','x_hy',),zlib=True)
        
        DS.createVariable('ssh_hycom','f',('time','y_hy','x_hy',),zlib=True)


        ## Write Data
        DS['time'][:] = np.array([(x - dt.datetime(2011,1,1,0,0,0)).total_seconds()/3600.0 for x in datetime_list[tt0:tt0+1]])
        DS['time'].units = 'hours since 2011-01-01 0000 UTC'
        DS['longitude_wrf'][:] = lon1p3km_entire_domain_keep
        DS['latitude_wrf'][:]  = lat1p3km_entire_domain_keep
        DS['terrain_wrf'][:]   = terrain
        DS['landmask_wrf'][:]  = landmask
        DS['rain_hourly'][:] = rain_collect[tt0,:,:]
        DS['olr'][:] = olr_collect[tt0,:,:]
        DS['psfc'][:] = psfc_collect[tt0,:,:]
        DS['u'][:] = u_collect[tt0,:,:]
        DS['v'][:] = v_collect[tt0,:,:]
        DS['wspd'][:] = wspd_collect[tt0,:,:]
        DS['wdir'][:] = wdir_collect[tt0,:,:]

        DS['longitude_umwm'][:] = lon_waves_keep
        DS['latitude_umwm'][:]  = lat_waves_keep
        DS['swh'][:] = swh[tt0,:,:]
        DS['mss'][:] = mss[tt0,:,:]
        DS['waves_seamask'][:]  = waves_seamask_keep
        
        DS['longitude_hycom'][:] = lon_hy_keep
        DS['latitude_hycom'][:]  = lat_hy_keep
        DS['bathymetry_hycom'][:]   = bathy_hy_keep
        
        DS['sst_hycom'][:] = sst_hy_keep[tt0,:,:]
        DS['sss_hycom'][:] = sss_hy_keep[tt0,:,:]
        
        DS['u_sfc_hycom'][:] = u_sfc_hy_keep[tt0,:,:]
        DS['v_sfc_hycom'][:] = v_sfc_hy_keep[tt0,:,:]
        DS['spd_sfc_hycom'][:] = np.sqrt(np.power(u_sfc_hy_keep[tt0,:,:],2) + np.power(v_sfc_hy_keep[tt0,:,:],2))
        DS['u_baro_hycom'][:] = u_baro_hy_keep[tt0,:,:]
        DS['v_baro_hycom'][:] = v_baro_hy_keep[tt0,:,:]
        DS['spd_baro_hycom'][:] = np.sqrt(np.power(u_baro_hy_keep[tt0,:,:],2) + np.power(v_baro_hy_keep[tt0,:,:],2))
        DS['spd_sfc_plus_baro_hycom'][:] = np.sqrt(np.power(u_sfc_hy_keep+u_baro_hy_keep,2) + np.power(v_sfc_hy_keep+v_baro_hy_keep,2))[tt0,:,:]
        
        DS['ssh_hycom'][:] = ssh_hy_keep[tt0,:,:] / 9.8

    
