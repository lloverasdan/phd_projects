#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""
import numpy as np
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
wavelength = 1000 # km
zl = 20 # km
cutoff = 30 # z grid point cutoff for inversion
xl, yl, x, y, xg, yg, dz, zp, facz = qgpv_invert.cartesian_mesh(nx, ny, nz, hres, zl)

### File paths
file_path_up = '/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_72h_10/dte_up'
file_path_low = '/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_72h_10/dte_low'

### Load the numpy files
u_up = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_72h_10/u_up.npy')
v_up = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_72h_10/v_up.npy')
th_up = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_72h_10/th_up.npy')

u_low = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_72h_10/u_low.npy')
v_low = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_72h_10/v_low.npy')
th_low = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_72h_10/th_low.npy')

### Compute DTE
dte_up = 0.5*(np.sum(u_up**2,axis=(1,2)) + np.sum(v_up**2,axis=(1,2))\
              + CP/TR*np.sum(th_up**2,axis=(1,2)))
dte_low = 0.5*(np.sum(u_low**2,axis=(1,2)) + np.sum(v_low**2,axis=(1,2))\
              + CP/TR*np.sum(th_low**2,axis=(1,2)))

### Save files
np.save(file_path_up,dte_up)
np.save(file_path_low,dte_low)