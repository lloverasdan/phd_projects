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

### Grid params
nx = 2000
ny = 1800
nz = 100
hres = 4 # km
zl = 20 # km
cutoff = 30 # z grid point cutoff for inversion
xl, yl, x, y, xg, yg, dz, zp, facz = qgpv_invert.cartesian_mesh(nx, ny, nz, hres, zl)

### Load the netCDF file
nc_in = Dataset('/p/work1/lloveras/adj_4km/wrf_output/wrfin_ctl')

### Read in the data and interpolate
th_r_z = np.asarray(vinterp(nc_in,getvar(nc_in,'th'),'ght_msl',zp/1000.))[:,0,0]
dt_r_dz = np.ones(nz)*np.nanmean(np.gradient(th_r_z,dz))
n_r = np.sqrt(G/np.nanmean(th_r_z)*dt_r_dz)

### Save output
np.save('/p/work1/lloveras/adj_4km/processed/pv_inversion/dt_r_dz',dt_r_dz)
np.save('/p/work1/lloveras/adj_4km/processed/pv_inversion/n_r',n_r)