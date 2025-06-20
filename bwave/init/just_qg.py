#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  18 09:18:04 2022

@author: lloverasdan
"""

import numpy as np
from bwave_ideal_wrf import qgpv_pert
from ambiance import Atmosphere

### Constants
G = 9.81  # gravitational acceleration, in m/s^2
T0 = 300.  # reference potential temperature, in K
P0 = 1.0e5  # reference pressure, in Pa
CP = 1004.  # specific heat at constant pressure, in J/(K*kg)       
CV = 717.  # specific heat at constant volume, in J/(K*kg)
RD = 287.  # ideal gas constant for dry air, in J/(K*kg)
RV = 461.6 # ideal gas constant for water vapor, in J/(K*kg)
F0 = 1.0e-4  # Coriolis parameter, in s^-1
SVPT0 = 273.15
GAMMA = CP/CV
KAPPA = RD/CP

### Grid parameters
nx = 400 # number of grid points in x direction
ny = 360 # number of grid points in y direction
nz = 100 # number of grid points in z direction
hres = 20. # horizontal grid resolution in km
zl = 20. # model top in km
dz = 200. # z grid spacing in m
ly = hres*ny*1000. # length in y direction in m

### Tropopause anomaly parameters
up_pert = True # tropopause anomaly option
up_pert_mag = 1.25e-4 # magnitude of tropopause anomaly in s^-1
x_up_pert = 140 # x center gridpoint for tropopause anomaly
y_up_pert = 165 # y center gridpoint for tropopause anomaly
z_up_pert = 40 # z center gridpoint for tropopause anomaly
ax_up_pert = 200. # x decay scale of tropopause anomaly in km
ay_up_pert = 600. # y decay scale of tropopause anomaly in km
az_up_pert = 1.5 # z decay scale of tropopause anomaly in km

### Surface anomaly parameters
surf_pert = True # surface anomaly option
surf_pert_mag = 4.0 # magnitude of surface theta anomaly in K
x_surf_pert = 200 # x center gridpoint for surface anomaly
y_surf_pert = 135 # y center gridpoint for surface anomaly
ax_surf_pert = 600. # x decay scale for surface anomaly in km
ay_surf_pert = 200. # y decay scale for surface anomaly in km

### Set up grids
xl, yl, x, y, xg, yg, dz, zp = qgpv_pert.cartesian_mesh(nx, ny, nz, hres, zl)
kmax, lmax, facx, facy, dxg, dyg = qgpv_pert.spectral_mesh(nx, ny, xl, yl)

### Compute reference static stability using standard atmosphere
atm_ref = Atmosphere(zp)
t_ref = atm_ref.temperature
p_ref = atm_ref.pressure
th_ref = (P0/p_ref)**(RD/CP)*t_ref

n_ref = np.zeros(nz-1)
dtdz_ref = np.zeros(nz-1)
for k in range(1,nz-1):
    dtdz_ref[k] = (th_ref[k+1] - th_ref[k-1])/(zp[k+1] - zp[k-1])
    n_ref[k] = np.sqrt(G/th_ref[k]*dtdz_ref[k])
    
dtdz_ref[0] = dtdz_ref[1]
n_ref[0] = n_ref[1]

### Tropopause anomaly
if up_pert:
    
    ### Initialize
    pvxy, bcxy, pvsp, b_fac = qgpv_pert.trop_anom(up_pert_mag,\
        x_up_pert, y_up_pert, z_up_pert, az_up_pert, ax_up_pert, ay_up_pert, x, y, zp, xg, yg, nx, ny, \
        nz, n_ref, dz, dtdz_ref)
    
    ### Invert
    fsp = qgpv_pert.qgpv_inversion(nx, ny, nz, b_fac, \
            facx, facy, kmax, lmax, pvsp, dz)

    ### Compute perturbations
    u_up_pert, v_up_pert, theta_up_pert, rho_up_pert, fxy_up = qgpv_pert.qgpv_solver(fsp, nx, ny, nz, \
            b_fac, dxg, dyg, dz, bcxy, dtdz_ref)

    p_up_pert = fxy_up*F0
    
else:
    u_up_pert = np.zeros((nz,ny,nx))
    v_up_pert = np.zeros((nz,ny,nx))
    theta_up_pert = np.zeros((nz,ny,nx))
    rho_up_pert = np.zeros((nz,ny,nx))
    p_up_pert = np.zeros((nz,ny,nx))

### Surface anomaly
if surf_pert:
    
    ### Initialize
    pvs, bcs, pvsps, b_facs = qgpv_pert.surf_anom(surf_pert_mag,\
        x_surf_pert, y_surf_pert, ax_surf_pert, ay_surf_pert, x, y, zp, xg, yg, nx, ny, nz,\
        n_ref, dz, dtdz_ref)
    
    ### Invert
    fsps = qgpv_pert.qgpv_inversion(nx, ny, nz, b_facs, \
        facx, facy, kmax, lmax, pvsps, dz)
    
    ### Compute perturbations
    u_surf_pert, v_surf_pert, theta_surf_pert, rho_surf_pert, fxy_surf = qgpv_pert.qgpv_solver(fsps, nx, ny, nz, \
        b_facs, dxg, dyg, dz, bcs, dtdz_ref)

    p_surf_pert = fxy_surf*F0

else:
    u_surf_pert = np.zeros((nz,ny,nx))
    v_surf_pert = np.zeros((nz,ny,nx))
    theta_surf_pert = np.zeros((nz,ny,nx))
    rho_surf_pert = np.zeros((nz,ny,nx))
    p_surf_pert = np.zeros((nz,ny,nx))

### Save the output    
np.save('/p/work1/lloveras/bwave/processed/qg_perts/u_up',u_up_pert)
np.save('/p/work1/lloveras/bwave/processed/qg_perts/v_up',v_up_pert)
np.save('/p/work1/lloveras/bwave/processed/qg_perts/th_up',theta_up_pert)
np.save('/p/work1/lloveras/bwave/processed/qg_perts/p_up',p_up_pert)

np.save('/p/work1/lloveras/bwave/processed/qg_perts/u_surf',u_surf_pert)
np.save('/p/work1/lloveras/bwave/processed/qg_perts/v_surf',v_surf_pert)
np.save('/p/work1/lloveras/bwave/processed/qg_perts/th_surf',theta_surf_pert)
np.save('/p/work1/lloveras/bwave/processed/qg_perts/p_surf',p_surf_pert)

np.save('/p/work1/lloveras/bwave/processed/qg_perts/z',zp)
np.save('/p/work1/lloveras/bwave/processed/qg_perts/n',n_ref)
np.save('/p/work1/lloveras/bwave/processed/qg_perts/dtdz',dtdz_ref)
