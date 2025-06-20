#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

from netCDF4 import Dataset
import numpy as np
from wrf import getvar, destagger

CP = 1005.7
RD = 287.04
P0 = 1000.
TR = 300.
LV = 2.501e6
EPS = 1.

### Input
nc1 = Dataset('/p/work1/lloveras/real/nov2018/4km_files/ctl/wrfout_d01_2018-11-13_12_00_00')
nc2 = Dataset('/p/work1/lloveras/real/nov2018/4km_files/adj_cape_box/wrfout_d01_2018-11-13_12_00_00')
lats = np.load('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/lats.npy')
lons = np.load('/p/work1/lloveras/real/nov2018/proc_oct/4km_files/lons.npy')

min_lat = 20
max_lat = 50
min_lon = -110
max_lon = -50

### Output
dte_file = '/p/work1/lloveras/real/nov2018/proc_oct/4km_files/adj_cape_box/dte_loc'

### Analysis
nt = 25
dte_vals = np.zeros(nt)
for ti in range(nt):
    du = np.asarray(getvar(nc1,'U',timeidx=ti)) - np.asarray(getvar(nc2,'U',timeidx=ti))
    dv = np.asarray(getvar(nc1,'V',timeidx=ti)) - np.asarray(getvar(nc2,'V',timeidx=ti))
    dt = np.asarray(getvar(nc1,'THM',timeidx=ti)) - np.asarray(getvar(nc2,'THM',timeidx=ti))
    dte = destagger(du,stagger_dim=2)**2 + destagger(dv,stagger_dim=1)**2 + (CP/TR)*dt**2
    
    dte[:,lons > max_lon] = np.NaN
    dte[:,lons < min_lon] = np.NaN
    dte[:,lats > max_lat] = np.NaN
    dte[:,lats < min_lat] = np.NaN
    
    dte_vals[ti] = np.nansum(dte)

np.save(dte_file,dte_vals)
