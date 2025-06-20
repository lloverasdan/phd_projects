#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates numpy files based on data processed from WRF output files
"""

import numpy as np
from netCDF4 import Dataset
from wrf import getvar, vinterp

### netCDF files

file = Dataset('/p/work1/lloveras/real/nov2018/4km_files/wave/wrfout_d01_2018-11-13_12_00_00')

### Height

z1 = np.asarray(getvar(file,'z',timeidx=4))
zl1 = np.asarray(vinterp(file, z1, 'pressure', [500], timeidx=4))

z2 = np.asarray(getvar(file,'z',timeidx=8))
zl2 = np.asarray(vinterp(file, z2, 'pressure', [500], timeidx=8))

z3 = np.asarray(getvar(file,'z',timeidx=12))
zl3 = np.asarray(vinterp(file, z3, 'pressure', [500], timeidx=12))

z4 = np.asarray(getvar(file,'z',timeidx=16))
zl4 = np.asarray(vinterp(file, z4, 'pressure', [500], timeidx=16))

### Wind

v1 = np.asarray(getvar(file,'va',timeidx=4))
vl1 = np.asarray(vinterp(file, v1, 'pressure', [500], timeidx=4))

v2 = np.asarray(getvar(file,'va',timeidx=8))
vl2 = np.asarray(vinterp(file, v2, 'pressure', [500], timeidx=8))

v3 = np.asarray(getvar(file,'va',timeidx=12))
vl3 = np.asarray(vinterp(file, v3, 'pressure', [500], timeidx=12))

v4 = np.asarray(getvar(file,'va',timeidx=16))
vl4 = np.asarray(vinterp(file, v4, 'pressure', [500], timeidx=16))

### PV

pvo = np.asarray(getvar(file,'pvo',timeidx=16))
pvl = np.asarray(vinterp(file, pvo, 'pressure', [400], timeidx=16))
zl = np.asarray(vinterp(file, z4, 'pressure', [400], timeidx=16))
pva = np.mean(np.asarray(vinterp(file, pvo, 'pressure', np.arange(300,510,10), timeidx=16)),axis=0)

### SLP

slp = np.asarray(getvar(file,'slp',timeidx=20))
slp2 = np.asarray(getvar(file,'slp',timeidx=24))

### dBZ

dbz1 = np.asarray(getvar(file,'mdbz',timeidx=4))
dbz2 = np.asarray(getvar(file,'mdbz',timeidx=8))
dbz3 = np.asarray(getvar(file,'mdbz',timeidx=12))
dbz4 = np.asarray(getvar(file,'mdbz',timeidx=16))
dbz5 = np.asarray(getvar(file,'mdbz',timeidx=20))
dbz6 = np.asarray(getvar(file,'mdbz',timeidx=24))

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

np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/z500_11-14_00',zl1[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/z500_11-14_12',zl2[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/z500_11-15_00',zl3[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/z500_11-15_12',zl4[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/v500_11-14_00',vl1[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/v500_11-14_12',vl2[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/v500_11-15_00',vl3[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/v500_11-15_12',vl4[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/dbz_11-14_00',dbz1)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/dbz_11-14_12',dbz2)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/dbz_11-15_00',dbz3)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/dbz_11-15_12',dbz4)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/dbz_11-16_00',dbz5)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/dbz_11-16_12',dbz6)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/z400_11-15_12',zl[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/pv400_11-15_12',pvl[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/pv300-500_11-15_12',pva)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/slp_11-16_00',slp)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/slp_11-16_12',slp2)
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/omg500_11-16_00',wl[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/omg500_11-16_12',wl2[0])
np.save('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/wave/precip-accu_11-16_12',pre)
