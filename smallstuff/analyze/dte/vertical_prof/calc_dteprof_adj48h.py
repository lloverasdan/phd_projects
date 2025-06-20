#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

from netCDF4 import Dataset
import numpy as np
from wrf import getvar, vinterp

### Constants
TR = 300.
CP = 1005.7

### Time step
ti = 16

### Load the netCDF files
nc_ctl = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_48h/wrfout_d01_2021-01-03_00_00_00')
nc_pert = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/adj_48h/wrfout_d01_2021-01-03_00_00_00')

### Read in the data
u_ctl = np.asarray(getvar(nc_ctl,'ua',timeidx=ti))
u_pert = np.asarray(getvar(nc_pert,'ua',timeidx=ti))

v_ctl = np.asarray(getvar(nc_ctl,'va',timeidx=ti))
v_pert = np.asarray(getvar(nc_pert,'va',timeidx=ti))

tk_ctl = np.asarray(getvar(nc_ctl,'tk',timeidx=ti))
tk_pert = np.asarray(getvar(nc_pert,'tk',timeidx=ti))

### Interpolate to height levels
levs = np.arange(950,40,-10)

u_ctl_p = vinterp(nc_ctl,u_ctl,'pressure',levs,timeidx=ti)
u_pert_p = vinterp(nc_pert,u_pert,'pressure',levs,timeidx=ti)

v_ctl_p = vinterp(nc_ctl,v_ctl,'pressure',levs,timeidx=ti)
v_pert_p = vinterp(nc_pert,v_pert,'pressure',levs,timeidx=ti)

tk_ctl_p = vinterp(nc_ctl,tk_ctl,'pressure',levs,timeidx=ti)
tk_pert_p = vinterp(nc_pert,tk_pert,'pressure',levs,timeidx=ti)

du = u_pert_p - u_ctl_p
dv = v_pert_p - v_ctl_p
dtk = tk_pert_p - tk_ctl_p

dte = np.nansum(du**2 + dv**2 + (CP/TR)*dtk**2,axis=(1,2))

np.save('/p/work1/lloveras/adj_4km/processed/dte/vertical_profiles/dte_adj48h',dte)
