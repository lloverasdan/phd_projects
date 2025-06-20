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

### Time step
ti = 16

### Grid params
nx = 2000
ny = 1800
nz = 100
hres = 4 # km
wavelength = 1000 # km
zl = 20 # km
cutoff = 30 # z grid point cutoff for inversion
xl, yl, x, y, xg, yg, dz, zp, facz = qgpv_invert.cartesian_mesh(nx, ny, nz, hres, zl)

### Load the netCDF files
nc_ctl = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_96h/wrfout_d01_2021-01-05_00_00_00')
nc_pert = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/wave_96h/wrfout_d01_2021-01-05_00_00_00')

### Path for output
file_path_up = '/p/work1/lloveras/adj_4km/processed/pv_inversion/wave_96h/dqgpv_up'
file_path_low = '/p/work1/lloveras/adj_4km/processed/pv_inversion/wave_96h/dqgpv_low'
file_path_ts = '/p/work1/lloveras/adj_4km/processed/pv_inversion/wave_96h/dts'

### Read in the reference profile
dt_r_dz = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/dt_r_dz.npy')

### Read in the data and interpolate
p_ctl = np.asarray(vinterp(nc_ctl,getvar(nc_ctl,'p',timeidx=ti),'ght_msl',zp/1000.,timeidx=ti))
p_pert = np.asarray(vinterp(nc_pert,getvar(nc_pert,'p',timeidx=ti),'ght_msl',zp/1000.,timeidx=ti))

th_ctl = np.asarray(vinterp(nc_ctl,getvar(nc_ctl,'th',timeidx=ti),'ght_msl',zp/1000.,timeidx=ti))
th_pert = np.asarray(vinterp(nc_pert,getvar(nc_pert,'th',timeidx=ti),'ght_msl',zp/1000.,timeidx=ti))

### Filter
p_ctl = qgpv_invert.horizontal_filter(p_ctl, nx, ny, hres*1000., wavelength*1000.)
p_pert = qgpv_invert.horizontal_filter(p_pert, nx, ny, hres*1000., wavelength*1000.)

th_ctl = qgpv_invert.horizontal_filter(th_ctl, nx, ny, hres*1000., wavelength*1000.)
th_pert = qgpv_invert.horizontal_filter(th_pert, nx, ny, hres*1000., wavelength*1000.)

ts_ctl = th_ctl[0,:,:]
ts_pert = th_pert[0,:,:]

### Balance
rho_ctl = P0/(RD*th_ctl)*(p_ctl/P0)**(CV/CP)
rho_pert = P0/(RD*th_pert)*(p_pert/P0)**(CV/CP)

dpdy_ctl, dpdx_ctl = np.gradient(p_ctl, hres*1000., axis=(1,2))
dpdy_pert, dpdx_pert = np.gradient(p_pert, hres*1000., axis=(1,2))

u_ctl = -dpdy_ctl/(rho_ctl*F0)
u_pert = -dpdy_pert/(rho_pert*F0)

v_ctl = dpdx_ctl/(rho_ctl*F0)
v_pert = dpdx_pert/(rho_pert*F0)

### Compute differences and QGPV
dth = th_pert - th_ctl
dts = ts_pert - ts_ctl

du = u_pert - u_ctl
dv = v_pert - v_ctl

dth_dz = np.gradient(dth, dz, axis=0)
du_dy = np.gradient(du, hres*1000., axis=1)
dv_dx = np.gradient(dv, hres*1000., axis=2)

dqgpv_up = np.zeros((nz,ny,nx))
dqgpv_low = np.zeros((nz,ny,nx))
for k in range(nz):
    if k <= cutoff:
        dqgpv_low[k,:,:] = dv_dx[k,:,:] - du_dy[k,:,:] + (F0/dt_r_dz[k])*dth_dz[k,:,:]
        dqgpv_up[k,:,:] = np.zeros((ny,nx))

    else:
        dqgpv_low[k,:,:] = np.zeros((ny,nx))
        dqgpv_up[k,:,:] = dv_dx[k,:,:] - du_dy[k,:,:] + (F0/dt_r_dz[k])*dth_dz[k,:,:]

### Save output
np.save(file_path_up,dqgpv_up)
np.save(file_path_low,dqgpv_low)
np.save(file_path_ts,dts)
