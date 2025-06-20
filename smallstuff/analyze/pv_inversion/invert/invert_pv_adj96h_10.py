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
kmax, lmax, facx, facy, dxg, dyg = qgpv_invert.spectral_mesh(nx, ny, xl, yl)

### File paths
file_path_u_up = '/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_96h_10/u_up'
file_path_v_up = '/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_96h_10/v_up'
file_path_th_up = '/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_96h_10/th_up'

file_path_u_low = '/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_96h_10/u_low'
file_path_v_low = '/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_96h_10/v_low'
file_path_th_low = '/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_96h_10/th_low'

### Load the numpy files
dqgpv_up = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_96h_10/dqgpv_up.npy')
dqgpv_low = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_96h_10/dqgpv_low.npy')
dts = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/adj_96h_10/dts.npy')
dt_r_dz = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/dt_r_dz.npy')
n_r = np.load('/p/work1/lloveras/adj_4km/processed/pv_inversion/n_r.npy')
bu_fac = (n_r[:-1]/F0)**2

### Initialize upper anomaly

### Boundary conditions (theta perturbation at model bottom and top)
ubcxy = np.zeros((ny,nx))
lbcxy = np.zeros((ny,nx))    

### Transform PV anomaly to spectral space in horizontal
lbcsp = np.fft.fft2(lbcxy)
ubcsp = np.fft.fft2(ubcxy)
pvsp = np.zeros((nz,ny,nx),dtype=np.dtype(np.complex64))
for k in range(0,nz):
    pvsp[k,:,:] = np.fft.fft2(np.squeeze(dqgpv_up[k,:,:]))

### Apply Neumann boundary conditions    
pvsp[0,:,:] = np.squeeze(pvsp[0,:,:]) + F0*(lbcsp/(facz[0]*dz*dt_r_dz[0]))
pvsp[nz-1,:,:] = np.squeeze(pvsp[nz-1,:,:]) - F0*(ubcsp/(facz[nz-1]*dz*dt_r_dz[nz-2]))

### Initialize lower anomaly

### Boundary conditions (theta perturbation at model bottom and top)
ubcs = np.zeros((ny,nx))
lbcs = np.zeros((ny,nx)) #dts

### Transform PV anomaly to spectral space in horizontal
lbcsps = np.fft.fft2(lbcs)
ubcsps = np.fft.fft2(ubcs)
pvsps = np.zeros((nz,ny,nx),dtype=np.dtype(np.complex64))
for k in range(0,nz):
    pvsps[k,:,:] = np.fft.fft2(np.squeeze(dqgpv_low[k,:,:]))

### Apply Neumann boundary conditions    
pvsps[0,:,:] = np.squeeze(pvsps[0,:,:]) + F0*(lbcsps/(facz[0]*dz*dt_r_dz[0]))
pvsps[nz-1,:,:] = np.squeeze(pvsps[nz-1,:,:]) - F0*(ubcsps/(facz[nz-1]*dz*dt_r_dz[nz-2]))

### Upper level
fbsp, ftsp, fzbsp, fztsp, fsp = qgpv_invert.qgpv_inversion(nx, ny, nz, bu_fac, \
        facx, facy, facz, kmax, lmax, pvsp, ubcsp, lbcsp, dz)

u_up_pert, v_up_pert, theta_up_pert, rho_up_pert, fxy_up = qgpv_invert.qgpv_solver(fsp, nx, ny, nz, \
        bu_fac, dxg, dyg, dz, lbcxy, ubcxy, dt_r_dz[:-1])

p_up_pert = fxy_up*F0

### Lower level
fbsps, ftsps, fzbsps, fztsps, fsps = qgpv_invert.qgpv_inversion(nx, ny, nz, bu_fac, \
    facx, facy, facz, kmax, lmax, pvsps, ubcsps, lbcsps, dz)

u_surf_pert, v_surf_pert, theta_surf_pert, rho_surf_pert, fxy_surf = qgpv_invert.qgpv_solver(fsps, nx, ny, nz, \
    bu_fac, dxg, dyg, dz, lbcs, ubcs, dt_r_dz[:-1])

p_surf_pert = fxy_surf*F0

### Save files
np.save(file_path_u_up,u_up_pert)
np.save(file_path_v_up,v_up_pert)
np.save(file_path_th_up,theta_up_pert)

np.save(file_path_u_low,u_surf_pert)
np.save(file_path_v_low,v_surf_pert)
np.save(file_path_th_low,theta_surf_pert)
