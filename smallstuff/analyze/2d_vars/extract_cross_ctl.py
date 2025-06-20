#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

from netCDF4 import Dataset
import numpy as np
from wrf import getvar, interplevel

### Times
ti = np.asarray([0,1,2,3,4,5,6]) # time indices in half hours
nx = 2000
nz = 100
yval = 765

### Load the netCDF files
nc1 = Dataset('/p/work1/lloveras/adj_4km/short_runs/ctl_48h/wrfout_d01_2021-01-03_00_00_00')

### Read in the data
dbz1 = np.zeros((len(ti),nz,nx))
t1 = np.zeros((len(ti),nz,nx))
z1 = np.zeros((len(ti),nz,nx))

for i in range(len(ti)):
    dbz1[i,:,:] = np.asarray(getvar(nc1,'dbz',timeidx=int(ti[i])))[:,yval,:]
    t1[i,:,:] = np.asarray(getvar(nc1,'T',timeidx=int(ti[i])))[:,yval,:]
    z1[i,:,:] = np.asarray(getvar(nc1,'z',timeidx=int(ti[i])))[:,yval,:]
    
### Save the output
np.save('/p/work1/lloveras/adj_4km/proc_short/dbz_ctl',dbz1)
np.save('/p/work1/lloveras/adj_4km/proc_short/th_ctl',t1)
np.save('/p/work1/lloveras/adj_4km/proc_short/z_ctl',z1)