#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates numpy files based on data processed from WRF output files
"""

import xarray as xr
import cfgrib
import numpy as np
from netCDF4 import Dataset
from wrf import getvar, vinterp
from scipy.interpolate import griddata

### Define functions for reading grib files

def get_slp(filepath, min_lat=-90, max_lat=90, min_lon=-180, max_lon=180):
    
    file = xr.open_dataset(filepath, engine='cfgrib',
                backend_kwargs={'filter_by_keys': {'typeOfLevel': 'meanSea','shortName': 'prmsl'}})
    var = file.get('prmsl').to_dataframe()
    latitudes = var.index.get_level_values('latitude')
    longitudes = var.index.get_level_values('longitude')
    map_function = lambda lon: (lon - 360) if (lon > 180) else lon
    remapped_longitudes = longitudes.map(map_function)
    var['longitude'] = remapped_longitudes
    var['latitude'] = latitudes
    lat_filter = (var['latitude'] >= min_lat) & (var['latitude'] <= max_lat)
    lon_filter = (var['longitude'] >= min_lon) & (var['longitude'] <= max_lon)
    var = var.loc[lat_filter & lon_filter]
    var = var.set_index(['latitude', 'longitude']).to_xarray()
    
    return var

### file paths

nov15_12z = '/p/work1/lloveras/real/nov2018/gfs_files/analysis/gfs.0p25.2018111512.f000.grib2'
nov16_00z = '/p/work1/lloveras/real/nov2018/gfs_files/analysis/gfs.0p25.2018111600.f000.grib2'
nov16_12z = '/p/work1/lloveras/real/nov2018/gfs_files/analysis/gfs.0p25.2018111612.f000.grib2'

slp1 = get_slp(nov15_12z)
slp2 = get_slp(nov16_00z)
slp3 = get_slp(nov16_12z)

### Lat-lon values

lons_4km = np.load('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/lons.npy')
lats_4km = np.load('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/lats.npy')

lons_gfs, lats_gfs = np.meshgrid(np.asarray(slp1['longitude']), np.asarray(slp1['latitude']))
lons_gfs[lons_gfs > 0] -= 360

### Interpolate

slp1_4km = griddata((lons_gfs.ravel(), lats_gfs.ravel()), 
                     np.asarray(slp1['prmsl']).ravel(), (lons_4km, lats_4km), method='linear')

slp2_4km = griddata((lons_gfs.ravel(), lats_gfs.ravel()), 
                     np.asarray(slp2['prmsl']).ravel(), (lons_4km, lats_4km), method='linear')

slp3_4km = griddata((lons_gfs.ravel(), lats_gfs.ravel()), 
                     np.asarray(slp3['prmsl']).ravel(), (lons_4km, lats_4km), method='linear')

### Save the output

np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/gfs_anl/slp_11-15_12',slp1_4km)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/gfs_anl/slp_11-16_00',slp2_4km)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/gfs_anl/slp_11-16_12',slp3_4km)
