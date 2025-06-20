#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates numpy files based on data processed from WRF output files
"""

import numpy as np
from netCDF4 import Dataset
from wrf import getvar, vinterp
from scipy.interpolate import griddata

### netCDF files

ctl = Dataset('/p/work1/lloveras/real/nov2018/4km_files/ctl/wrfout_d01_2018-11-13_12_00_00')
gfs = Dataset('/p/work1/lloveras/real/nov2018/30km_files/gfs/wrfin_d01_2018-11-15_12_00_00')

### Lat-lon
lats_4km = np.asarray(getvar(ctl,'lat'))
lons_4km = np.asarray(getvar(ctl,'lon'))
lons_4km[lons_4km > 0] -= 360

lats_30km = np.asarray(getvar(gfs,'lat'))
lons_30km = np.asarray(getvar(gfs,'lon'))
lons_30km[lons_30km > 0] -= 360

### Settings

ti = 0
plevs = np.arange(300,510,10)
plev = [400]

### Read in the data

pv_gfs = np.asarray(getvar(gfs,'pvo',timeidx=ti))
z_gfs = np.asarray(getvar(gfs,'z',timeidx=ti))

### Interpolate onto levels

pva_gfs = np.mean(np.asarray(vinterp(gfs, pv_gfs, 'pressure', plevs, timeidx=0)),axis=0)
pvl_gfs = np.asarray(vinterp(gfs, pv_gfs, 'pressure', plev, timeidx=0))
zl_gfs = np.asarray(vinterp(gfs, z_gfs, 'pressure', plev, timeidx=0))

### Horizontally interpolate onto 4-km mesh

pva_gfs_4km = griddata((lons_30km.ravel(), lats_30km.ravel()), 
                             pva_gfs.ravel(), (lons_4km, lats_4km), method='linear')

pvl_gfs_4km = griddata((lons_30km.ravel(), lats_30km.ravel()), 
                             pvl_gfs.ravel(), (lons_4km, lats_4km), method='linear')

zl_gfs_4km = griddata((lons_30km.ravel(), lats_30km.ravel()), 
                             zl_gfs.ravel(), (lons_4km, lats_4km), method='linear')

### Save the output

np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/gfs_anl/pv300-500_11-15_12',pva_gfs_4km)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/gfs_anl/pv400_11-15_12',pvl_gfs_4km)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/gfs_anl/z400_11-15_12',zl_gfs_4km)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/lats',lats_4km)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/lons',lons_4km)
