#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates numpy files based on data processed from WRF output files
"""

import numpy as np
from netCDF4 import Dataset
from wrf import getvar, vinterp

### netCDF files

file = Dataset('/p/work1/lloveras/real/nov2018/4km_files/adj_hole/wrfout_d01_2018-11-13_12_00_00')

### Height

z = np.asarray(getvar(file,'z',timeidx=16))
zl = np.asarray(vinterp(file, z, 'pressure', [500], timeidx=16))

### PV

pvo = np.asarray(getvar(file,'pvo',timeidx=16))
zl2 = np.asarray(vinterp(file, z, 'pressure', [400], timeidx=16))
pvl = np.asarray(vinterp(file, pvo, 'pressure', [400], timeidx=16))
pva = np.mean(np.asarray(vinterp(file, pvo, 'pressure', np.arange(300,510,10), timeidx=16)),axis=0)

### SLP

slp = np.asarray(getvar(file,'slp',timeidx=20))
slp2 = np.asarray(getvar(file,'slp',timeidx=24))

### dBZ

dbz = np.asarray(getvar(file,'mdbz',timeidx=20))
dbz2 = np.asarray(getvar(file,'mdbz',timeidx=24))

### Omega

w = np.asarray(getvar(file,'omg',timeidx=20))
wl = np.asarray(vinterp(file, w, 'pressure', [500], timeidx=20))

w2 = np.asarray(getvar(file,'omg',timeidx=24))
wl2 = np.asarray(vinterp(file, w2, 'pressure', [500], timeidx=24))

### Precip

pre = (np.asarray(getvar(file,'RAINC',timeidx=24)) +\
                 np.asarray(getvar(file,'RAINNC',timeidx=24)) +\
                 np.asarray(getvar(file,'SNOWNC',timeidx=24)) +\
                 np.asarray(getvar(file,'HAILNC',timeidx=24)) +\
                 np.asarray(getvar(file,'GRAUPELNC',timeidx=24)))

### Save output

np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/z500_11-15_12',zl[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/z400_11-15_12',zl2[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/pv400_11-15_12',pvl[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/pv300-500_11-15_12',pva)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/slp_11-16_00',slp)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/dbz_11-16_00',dbz)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/omg500_11-16_00',wl[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/slp_11-16_12',slp2)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/dbz_11-16_12',dbz2)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/omg500_11-16_12',wl2[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_hole/precip-accu_11-16_12',pre)
