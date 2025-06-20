#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  18 09:18:04 2022

@author: lloverasdan
"""

import numpy as np
import netCDF4 as nc
import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from bwave_ideal_wrf import wrf_fields, write_wrfinputd01
from wrf import getvar, destagger, interplevel

CP = 1005.7
RD = 287.04
P0 = 1000.
TR = 300.
LV = 2.501e6
EPS = 1.

nx_coamps = 250
ny_coamps = 225
dx_coamps = 32000.
nz_coamps = 100
nx_wrf = 2000
ny_wrf = 1800
nz_wrf = 100
dx_wrf = 4000.
cent_length = 2500000.
shape_coamps = [nz_coamps, ny_coamps, nx_coamps]
shape_wrf = [nz_wrf, ny_wrf, nx_wrf]

coamps_dir = '/p/work1/lloveras/adj_4km/input_data/files_from_jim/adj_72h/'
adj_dir = '/p/work1/lloveras/adj_4km/input_data/rst_files/adj_72h_100/wrfrst_d01_2021-01-04_00_00_00'

coamps_u = np.fromfile(coamps_dir + 'aaauu1_sig_020007_000050_1a0250x0225_2022070500_00000000_fcstfld', dtype='>f4')
coamps_v = np.fromfile(coamps_dir + 'aaavv1_sig_020007_000050_1a0250x0225_2022070500_00000000_fcstfld', dtype='>f4')
coamps_ex = np.fromfile(coamps_dir + 'aaapp1_sig_020007_000050_1a0250x0225_2022070500_00000000_fcstfld', dtype='>f4')
coamps_th = np.fromfile(coamps_dir + 'aaath1_sig_020007_000050_1a0250x0225_2022070500_00000000_fcstfld', dtype='>f4')
coamps_qv = np.fromfile(coamps_dir + 'aaaqv1_sig_020007_000050_1a0250x0225_2022070500_00000000_fcstfld', dtype='>f4')
coamps_full_p = np.fromfile(coamps_dir + 'ttlprs_sig_020007_000050_1a0250x0225_2022070500_00000000_fcstfld', dtype='>f4')

coamps_u = np.flip(np.reshape(coamps_u, shape_coamps),axis=0).astype(float)
coamps_v = np.flip(np.reshape(coamps_v, shape_coamps),axis=0).astype(float)
coamps_ex = np.flip(np.reshape(coamps_ex, shape_coamps),axis=0).astype(float)
coamps_th = np.flip(np.reshape(coamps_th, shape_coamps),axis=0).astype(float)
coamps_qv = np.flip(np.reshape(coamps_qv, shape_coamps),axis=0).astype(float)
coamps_full_p = np.flip(np.reshape(coamps_full_p, shape_coamps),axis=0).astype(float)

coamps_pert_ex = (coamps_full_p/P0)**(RD/CP) + coamps_ex
coamps_pert_p = P0*coamps_pert_ex**(CP/RD)
coamps_p = (coamps_pert_p - coamps_full_p)*100.

x_coamps = np.linspace(int(dx_coamps),int(nx_coamps*dx_coamps),int(nx_coamps))
x_wrf = np.linspace(int(dx_wrf),int(nx_wrf*dx_wrf),int(nx_wrf))

coamps_u_temp = np.zeros([nz_coamps, ny_coamps, nx_wrf])
coamps_v_temp = np.zeros([nz_coamps, ny_coamps, nx_wrf])
coamps_p_temp = np.zeros([nz_coamps, ny_coamps, nx_wrf])
coamps_th_temp = np.zeros([nz_coamps, ny_coamps, nx_wrf])
coamps_qv_temp = np.zeros([nz_coamps, ny_coamps, nx_wrf])
coamps_full_p_temp = np.zeros([nz_coamps, ny_coamps, nx_wrf])
for k in range(nz_coamps):
    for j in range(ny_coamps):
        for i in range(nx_wrf):
            coamps_u_temp[k,j,i] = wrf_fields.interp_0(coamps_u[k,j,:], x_coamps, x_wrf[i], nx_coamps)
            coamps_v_temp[k,j,i] = wrf_fields.interp_0(coamps_v[k,j,:], x_coamps, x_wrf[i], nx_coamps)
            coamps_p_temp[k,j,i] = wrf_fields.interp_0(coamps_p[k,j,:], x_coamps, x_wrf[i], nx_coamps)
            coamps_th_temp[k,j,i] = wrf_fields.interp_0(coamps_th[k,j,:], x_coamps, x_wrf[i], nx_coamps)
            coamps_qv_temp[k,j,i] = wrf_fields.interp_0(coamps_qv[k,j,:], x_coamps, x_wrf[i], nx_coamps)
            coamps_full_p_temp[k,j,i] = wrf_fields.interp_0(coamps_full_p[k,j,:], x_coamps, x_wrf[i], nx_coamps)

y_coamps = np.linspace(int(dx_coamps),int(ny_coamps*dx_coamps),int(ny_coamps))
y_wrf = np.linspace(int(dx_wrf),int(ny_wrf*dx_wrf),int(ny_wrf))

coamps_u_4km = np.zeros([nz_coamps, ny_wrf, nx_wrf])
coamps_v_4km = np.zeros([nz_coamps, ny_wrf, nx_wrf])
coamps_p_4km = np.zeros([nz_coamps, ny_wrf, nx_wrf])
coamps_th_4km = np.zeros([nz_coamps, ny_wrf, nx_wrf])
coamps_qv_4km = np.zeros([nz_coamps, ny_wrf, nx_wrf])
coamps_full_p_4km = np.zeros([nz_coamps, ny_wrf, nx_wrf])
for k in range(nz_coamps):
    for j in range(ny_wrf):
        for i in range(nx_wrf):
            coamps_u_4km[k,j,i] = wrf_fields.interp_0(coamps_u_temp[k,:,i], y_coamps, y_wrf[j], ny_coamps)
            coamps_v_4km[k,j,i] = wrf_fields.interp_0(coamps_v_temp[k,:,i], y_coamps, y_wrf[j], ny_coamps)
            coamps_p_4km[k,j,i] = wrf_fields.interp_0(coamps_p_temp[k,:,i], y_coamps, y_wrf[j], ny_coamps)
            coamps_th_4km[k,j,i] = wrf_fields.interp_0(coamps_th_temp[k,:,i], y_coamps, y_wrf[j], ny_coamps)
            coamps_qv_4km[k,j,i] = wrf_fields.interp_0(coamps_qv_temp[k,:,i], y_coamps, y_wrf[j], ny_coamps)
            coamps_full_p_4km[k,j,i] = wrf_fields.interp_0(coamps_full_p_temp[k,:,i], y_coamps, y_wrf[j], ny_coamps)

ncfile = nc.Dataset(adj_dir,'r+')
wrf_u = getvar(ncfile,'U_2')
wrf_v = getvar(ncfile,'V_2')
wrf_th = getvar(ncfile,'T_2')
wrf_p = getvar(ncfile,'P')
wrf_pb = getvar(ncfile,'PB')
wrf_qv = getvar(ncfile,'QVAPOR')
wrf_full_p = (np.asarray(wrf_p) + np.asarray(wrf_pb))/100.

psurf = wrf_p[0,:,:]
minp_ind = np.unravel_index(psurf.argmin(), psurf.shape)
cent_ind =  int(cent_length/dx_wrf)
rol_val = -(cent_ind - minp_ind[1])

coamps_u_4km_roll = np.roll(coamps_u_4km, rol_val, axis=-1)
coamps_v_4km_roll = np.roll(coamps_v_4km, rol_val, axis=-1)
coamps_p_4km_roll = np.roll(coamps_p_4km, rol_val, axis=-1)
coamps_th_4km_roll = np.roll(coamps_th_4km, rol_val, axis=-1)
coamps_qv_4km_roll = np.roll(coamps_qv_4km, rol_val, axis=-1)
coamps_full_p_4km_roll = np.roll(coamps_full_p_4km, rol_val, axis=-1)

u_pert = np.zeros([nz_wrf,ny_wrf,nx_wrf+1])
v_pert = np.zeros([nz_wrf,ny_wrf+1,nx_wrf])
p_pert = np.zeros(shape_wrf)
th_pert = np.zeros(shape_wrf)
qv_pert = np.zeros(shape_wrf)
for k in range(nz_wrf):
    u_pert[k,:,:-1] = interplevel(coamps_u_4km_roll, coamps_full_p_4km_roll, wrf_full_p[k,:,:])
    v_pert[k,:-1,:] = interplevel(coamps_v_4km_roll, coamps_full_p_4km_roll, wrf_full_p[k,:,:])
    th_pert[k,:,:] = interplevel(coamps_th_4km_roll, coamps_full_p_4km_roll, wrf_full_p[k,:,:])
    p_pert[k,:,:] = interplevel(coamps_p_4km_roll, coamps_full_p_4km_roll, wrf_full_p[k,:,:])
    qv_pert[k,:,:] = interplevel(coamps_qv_4km_roll, coamps_full_p_4km_roll, wrf_full_p[k,:,:])
    
u_pert[:,:,-1] = u_pert[:,:,0]
v_pert[:,-1,:] = v_pert[:,-2,:]

u_pert = np.nan_to_num(u_pert)
v_pert = np.nan_to_num(v_pert)
th_pert = np.nan_to_num(th_pert)
p_pert = np.nan_to_num(p_pert)
qv_pert = np.nan_to_num(qv_pert)

ncfile.variables['U_1'][0,:,:,:] = wrf_u + u_pert/100.
ncfile.variables['V_1'][0,:,:,:] = wrf_v + v_pert/100.
ncfile.variables['T_1'][0,:,:,:] = wrf_th + th_pert/100.
ncfile.variables['P'][0,:,:,:] = wrf_p + p_pert/100.
ncfile.variables['QVAPOR'][0,:,:,:] = wrf_qv + qv_pert/100.

ncfile.variables['U_2'][0,:,:,:] = wrf_u + u_pert/100.
ncfile.variables['V_2'][0,:,:,:] = wrf_v + v_pert/100.
ncfile.variables['T_2'][0,:,:,:] = wrf_th + th_pert/100.

ncfile.close()
