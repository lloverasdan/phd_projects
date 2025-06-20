#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  18 09:18:04 2022

@author: lloverasdan
"""

import numpy as np
from netCDF4 import Dataset
from bwave_ideal_wrf import epv_jet, qgpv_pert
from wrf import interplevel, destagger, getvar

CP = 1005.7
RD = 287.04
P0 = 1000.
TR = 300.
LV = 2.501e6
EPS = 1.

dx = 32000.
dy = 32000.
nx = 250
ny = 225
nz_coamps = 100
nz_wrf = 100
shape_coamps = [nz_coamps, ny, nx]
shape_wrf = [nz_wrf, ny, nx]
ti = 16 # noshear: 2 days = 8, anticyc: 4 days = 16
rol = -125  # noshear: 2560 km = -80, anticyc: 4000 km = -125

rst_dir = '/p/work/lloveras/bwave/32km_files/anticyc/restart_moist_adj/wrfrst_d01_2021-01-05_00_00_00'
out_dir = '/p/work/lloveras/bwave/32km_files/anticyc/output/wrfout_anticyc_moist'
coamps_dir = '/p/work/lloveras/bwave/32km_files/jim_files/files_anticyc/'

ncfile = Dataset(rst_dir,'r+')
outfile = Dataset(out_dir)
wrf_u = np.roll(np.asarray(getvar(ncfile,'U_2')),rol,axis=-1)
wrf_v = np.roll(np.asarray(getvar(ncfile,'V_2')),rol,axis=-1)
wrf_th = np.roll(np.asarray(getvar(ncfile,'T_2')),rol,axis=-1)
wrf_p = np.roll(np.asarray(getvar(ncfile,'P')),rol,axis=-1)
wrf_qv = np.roll(np.asarray(getvar(ncfile,'QVAPOR')),rol,axis=-1)
wrf_full_p = np.roll(np.asarray(getvar(outfile,'pressure',timeidx=ti)),rol,axis=-1)

coamps_u = np.fromfile(coamps_dir + 'aaauu1_sig_020007_000050_1a0250x0225_2023040400_00000000_fcstfld', dtype='>f4')
coamps_v = np.fromfile(coamps_dir + 'aaavv1_sig_020007_000050_1a0250x0225_2023040400_00000000_fcstfld', dtype='>f4')
coamps_ex = np.fromfile(coamps_dir + 'aaapp1_sig_020007_000050_1a0250x0225_2023040400_00000000_fcstfld', dtype='>f4')
coamps_th = np.fromfile(coamps_dir + 'aaath1_sig_020007_000050_1a0250x0225_2023040400_00000000_fcstfld', dtype='>f4')
coamps_qv = np.fromfile(coamps_dir + 'aaaqv1_sig_020007_000050_1a0250x0225_2023040400_00000000_fcstfld', dtype='>f4')
coamps_full_p = np.fromfile(coamps_dir + 'ttlprs_sig_020007_000050_1a0250x0225_2023040400_00000000_fcstfld', dtype='>f4')

coamps_u = np.flip(np.reshape(coamps_u, shape_coamps),axis=0)
coamps_v = np.flip(np.reshape(coamps_v, shape_coamps),axis=0)
coamps_ex = np.flip(np.reshape(coamps_ex, shape_coamps),axis=0)
coamps_th = np.flip(np.reshape(coamps_th, shape_coamps),axis=0)
coamps_qv = np.flip(np.reshape(coamps_qv, shape_coamps),axis=0)
coamps_full_p = np.flip(np.reshape(coamps_full_p, shape_coamps),axis=0)

coamps_pert_ex = (coamps_full_p/P0)**(RD/CP) + coamps_ex
coamps_pert_p = P0*coamps_pert_ex**(CP/RD)
coamps_p = (coamps_pert_p - coamps_full_p)*100.

u_pert = np.zeros([nz_wrf,ny,nx+1])
v_pert = np.zeros([nz_wrf,ny+1,nx])
p_pert = np.zeros(shape_wrf)
th_pert = np.zeros(shape_wrf)
qv_pert = np.zeros(shape_wrf)
for k in range(nz_wrf):
    u_pert[k,:,:-1] = interplevel(coamps_u, coamps_full_p, wrf_full_p[k,:,:])
    v_pert[k,:-1,:] = interplevel(coamps_v, coamps_full_p, wrf_full_p[k,:,:])
    th_pert[k,:,:] = interplevel(coamps_th, coamps_full_p, wrf_full_p[k,:,:])
    p_pert[k,:,:] = interplevel(coamps_p, coamps_full_p, wrf_full_p[k,:,:])
    qv_pert[k,:,:] = interplevel(coamps_qv, coamps_full_p, wrf_full_p[k,:,:])
    
u_pert[:,:,-1] = u_pert[:,:,0]
v_pert[:,-1,:] = v_pert[:,-2,:]

u_pert = np.nan_to_num(u_pert)
v_pert = np.nan_to_num(v_pert)
th_pert = np.nan_to_num(th_pert)
p_pert = np.nan_to_num(p_pert)
qv_pert = np.nan_to_num(qv_pert)

ncfile.variables['U_2'][0,:,:,:] = np.roll(wrf_u + u_pert,-rol,axis=-1)
ncfile.variables['V_2'][0,:,:,:] = np.roll(wrf_v + v_pert,-rol,axis=-1)
ncfile.variables['T_2'][0,:,:,:] = np.roll(wrf_th + th_pert,-rol,axis=-1)
ncfile.variables['P'][0,:,:,:] = np.roll(wrf_p + p_pert,-rol,axis=-1)
ncfile.variables['QVAPOR'][0,:,:,:] = np.roll(wrf_qv + qv_pert,-rol,axis=-1)

ncfile.close()