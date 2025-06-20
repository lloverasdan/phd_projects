#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  18 09:18:04 2022

@author: lloverasdan
"""

import numpy as np
from netCDF4 import Dataset
from wrf import vinterp, getvar
import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from bwave_ideal_wrf import qgpv_pert

### Define constants
nx = 2000
ny = 1800
nz = 90
ti = 8
hres = 4
G = 9.81
F0 = 1.0e-4
RHO = 1.0
zbot = 0.2
ztop = 18
z = np.linspace(zbot,ztop,nz)

### Read in WRF data and interpolate
nc1 = Dataset('/p/work1/lloveras/adj_4km/long_runs/ctl_48h/wrfout_d01_2021-01-05_00_00_00')
nc2 = Dataset('/p/work1/lloveras/adj_4km/long_runs/adj_48h_10/wrfout_d01_2021-01-05_00_00_00')

pv1 = np.asarray(getvar(nc1,'pvo',timeidx=ti))*1e-6
pv2 = np.asarray(getvar(nc2,'pvo',timeidx=ti))*1e-6

u1 = np.asarray(getvar(nc1,'ua',timeidx=ti))
u2 = np.asarray(getvar(nc2,'ua',timeidx=ti))

v1 = np.asarray(getvar(nc1,'va',timeidx=ti))
v2 = np.asarray(getvar(nc2,'va',timeidx=ti))

th1 = np.asarray(getvar(nc1,'th',timeidx=ti))

pv1 = np.asarray(vinterp(nc1,pv1,'ght_msl',z,timeidx=ti))
pv2 = np.asarray(vinterp(nc2,pv2,'ght_msl',z,timeidx=ti))

u1 = np.asarray(vinterp(nc1,u1,'ght_msl',z,timeidx=ti))
u2 = np.asarray(vinterp(nc2,u2,'ght_msl',z,timeidx=ti))

v1 = np.asarray(vinterp(nc1,v1,'ght_msl',z,timeidx=ti))
v2 = np.asarray(vinterp(nc2,v2,'ght_msl',z,timeidx=ti))

th1 = np.asarray(vinterp(nc1,th1,'ght_msl',z,timeidx=ti))

### Compute reference state
th_avg = np.nanmean(th1,axis=(1,2))
dtdz_avg = np.zeros(nz)
n_avg = np.zeros(nz-1)
for k in range(1,nz-1):
    dtdz_avg[k] = (th_avg[k+1] - th_avg[k-1])/((z[k+1] - z[k-1])*1000.)
    n_avg[k] = np.sqrt(G/th_avg[k]*dtdz_avg[k])

dtdz_avg[0] = dtdz_avg[1]
dtdz_avg[-1] = dtdz_avg[-2]
n_avg[0] = np.sqrt(G/th_avg[0]*dtdz_avg[0])

bu = n_avg/F0
bu_fac = bu**2

### Compute QGPV
qgpv1 = np.zeros((nz,ny,nx))
qgpv2 = np.zeros((nz,ny,nx))
for k in range(nz):
    qgpv1[k,:,:] = RHO*pv1[k,:,:]/dtdz_avg[k]
    qgpv2[k,:,:] = RHO*pv2[k,:,:]/dtdz_avg[k]
    
dtdz_avg = dtdz_avg[:-1]

### Compute difference fields and transform into spectral space
du = u2 - u1
dv = v2 - v1
dqgpv = qgpv2 - qgpv1
pvsp = np.zeros((nz,ny,nx),dtype=np.dtype(np.complex64))
for k in range(nz):
    pvsp[k,:,:] = np.fft.fft2(np.squeeze(dqgpv[k,:,:]))

### Invert
lbcxy = np.zeros((ny,nx))
ubcxy = np.zeros((ny,nx))
lbcsp = np.fft.fft2(lbcxy)
ubcsp = np.fft.fft2(ubcxy)
xl, yl, x, y, xg, yg, dz, zp, facz = qgpv_pert.cartesian_mesh(nx, ny, nz, hres, ztop)
kmax, lmax, facx, facy, dxg, dyg = qgpv_pert.spectral_mesh(nx, ny, xl, yl)

fbsp, ftsp, fzbsp, fztsp, fsp = qgpv_pert.qgpv_inversion(nx, ny, nz, bu_fac, \
        facx, facy, facz, kmax, lmax, pvsp, ubcsp, lbcsp, dz)

du_qg, dv_qg, dth_qg, dr_qg, dphi_qg = qgpv_pert.qgpv_solver(fsp, nx, ny, nz, \
        bu_fac, dxg, dyg, dz, lbcxy, ubcxy, dtdz_avg)

### Save files
np.save('/p/work1/lloveras/adj_4km/processed_data/qgpv_96h/dqgpv_adj_10',dqgpv)
np.save('/p/work1/lloveras/adj_4km/processed_data/qgpv_96h/du_adj_10',du)
np.save('/p/work1/lloveras/adj_4km/processed_data/qgpv_96h/dv_adj_10',dv)
np.save('/p/work1/lloveras/adj_4km/processed_data/qgpv_96h/du_adj_qg_10',du_qg)
np.save('/p/work1/lloveras/adj_4km/processed_data/qgpv_96h/dv_adj_qg_10',dv_qg)
