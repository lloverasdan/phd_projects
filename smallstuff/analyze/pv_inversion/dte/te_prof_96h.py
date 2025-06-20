#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

from netCDF4 import Dataset
import numpy as np
from wrf import getvar, vinterp
from inversion_tools import qgpv_invert

### Constants
G = 9.81  # gravitational acceleration, in m/s^2
T0 = 300.  # reference potential temperature, in K
TR = 300.
P0 = 1.0e5  # reference pressure, in Pa
CP = 1004.  # specific heat at constant pressure, in J/(K*kg)       
CV = 717.  # specific heat at constant volume, in J/(K*kg)
RD = 287.  # ideal gas constant for dry air, in J/(K*kg)
RV = 461.6 # ideal gas constant for water vapor, in J/(K*kg)
F0 = 1.0e-4  # Coriolis parameter, in s^-1
SVPT0 = 273.15
GAMMA = CP/CV
KAPPA = RD/CP

### Time index
ti = 16

### Grid params
nx = 2000
ny = 1800
nz = 100
hres = 4 # km
zl = 20 # km
cutoff = 30 # z grid point cutoff for inversion
xl, yl, x, y, xg, yg, dz, zp, facz = qgpv_invert.cartesian_mesh(nx, ny, nz, hres, zl)

### Load the netCDF files
nc_in = Dataset('/p/work1/lloveras/adj_4km/wrf_output/wrfin_ctl')
nc_out = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_96h/wrfout_d01_2021-01-05_00_00_00')

### Input
u_in = np.asarray(vinterp(nc_in,getvar(nc_in,'ua'),'ght_msl',zp/1000.))
v_in = np.asarray(vinterp(nc_in,getvar(nc_in,'va'),'ght_msl',zp/1000.))
th_in = np.asarray(vinterp(nc_in,getvar(nc_in,'th'),'ght_msl',zp/1000.))

u_in = np.tile(np.expand_dims(u_in[:,:,0],-1),nx)
v_in = np.tile(np.expand_dims(v_in[:,:,0],-1),nx)
th_in = np.tile(np.expand_dims(th_in[:,:,0],-1),nx)

### Output
u_out = np.asarray(vinterp(nc_out,getvar(nc_out,'ua',timeidx=ti),'ght_msl',zp/1000.,timeidx=ti))
v_out = np.asarray(vinterp(nc_out,getvar(nc_out,'va',timeidx=ti),'ght_msl',zp/1000.,timeidx=ti))
th_out = np.asarray(vinterp(nc_out,getvar(nc_out,'th',timeidx=ti),'ght_msl',zp/1000.,timeidx=ti))

### Compute DTE
du = u_out - u_in
dv = v_out - v_in
dth = th_out - th_in

dte = 0.5*(np.sum(du**2,axis=(1,2)) + np.sum(dv**2,axis=(1,2)) + CP/TR*np.sum(dth**2,axis=(1,2)))

### Save file
np.save('/p/work1/lloveras/adj_4km/processed/pv_inversion/te_96h',dte)
