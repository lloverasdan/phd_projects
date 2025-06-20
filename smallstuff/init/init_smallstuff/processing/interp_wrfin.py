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
from wrf import getvar, destagger

# DEFINE

nc_in = nc.Dataset('/p/work1/lloveras/48h_adj_30km/in_files/wrfin_ctl')
title_str = "OUTPUT FROM IDEAL V3.6.1 PREPROCESSOR"
time_str = "2021-01-01_00:00:00"

nx_in = 500
nx_in_trunc = 267
ny_in = 240
dx_in = 30.
hres_m_in = dx_in*1000.

nx_out = 2000
ny_out = 1800
dx_out = 4.
hres_m_out = dx_out*1000.
rdx_out = 1/hres_m_out
rdy_out = 1/hres_m_out

dt = 24.
mp_physics = 17
ra_lw_physics = 0
ra_sw_physics = 0
sf_sfclay_physics = 1
sf_surface_physics = 0
bl_pbl_physics = 1
cu_physics = 0
num_soil_layers = 5
diff_opt = 1
km_opt = 4
damp_opt = 3
dampcoef = 0.4
khdif = 0.
kvdif = 0.
spec_bdy_width = 8
sf_lake_physics = 0
surface_input_source = 1
hypsometric_opt = 0
num_land_cat = 24
num_soil_cat = 16
ptop = 5000.
nz = 100

# READ

u_field_in = destagger(np.array(getvar(nc_in,'U')),stagger_dim=2)[:,:,:nx_in_trunc]
v_field_in = destagger(np.array(getvar(nc_in,'V')),stagger_dim=1)[:,:,:nx_in_trunc]
t_in = np.array(getvar(nc_in,'T'))[:,:,:nx_in_trunc]
ph_in = np.array(getvar(nc_in,'PH'))[:,:,:nx_in_trunc]
phb_in = np.array(getvar(nc_in,'PHB'))[:,:,:nx_in_trunc]
p_in = np.array(getvar(nc_in,'P'))[:,:,:nx_in_trunc]
pb_in = np.array(getvar(nc_in,'PB'))[:,:,:nx_in_trunc]
t_init_in = np.array(getvar(nc_in,'T_INIT'))[:,:,:nx_in_trunc]
mu_in = np.array(getvar(nc_in,'MU'))[:,:nx_in_trunc]
mub_in = np.array(getvar(nc_in,'MUB'))[:,:nx_in_trunc]
moist_in = np.array(getvar(nc_in,'QVAPOR'))[:,:,:nx_in_trunc]
tsk_in = np.array(getvar(nc_in,'TSK'))[:,:nx_in_trunc]
tmn_in = np.array(getvar(nc_in,'TMN'))[:,:nx_in_trunc]

fnm = np.array(getvar(nc_in,'FNM'))
fnp = np.array(getvar(nc_in,'FNP'))
rdnw = np.array(getvar(nc_in,'RDNW'))
rdn = np.array(getvar(nc_in,'RDN'))
dnw = np.array(getvar(nc_in,'DNW'))
dn = np.array(getvar(nc_in,'DN'))
t_base = np.array(getvar(nc_in,'T_BASE'))
cfn = float(getvar(nc_in,'CFN'))
cfn1 = float(getvar(nc_in,'CFN1'))
cf1 = float(getvar(nc_in,'CF1'))
cf2 = float(getvar(nc_in,'CF2'))
cf3 = float(getvar(nc_in,'CF3'))
u_base = np.array(getvar(nc_in,'U_BASE'))
v_base = np.array(getvar(nc_in,'V_BASE'))
qv_base = np.array(getvar(nc_in,'QV_BASE'))
znw = np.array(getvar(nc_in,'ZNW'))
znu = np.array(getvar(nc_in,'ZNU'))

# INTERPOLATE

x_in = np.linspace(int(dx_in),int(nx_in_trunc*dx_in),int(nx_in_trunc))
x_out = np.linspace(int(dx_out),int(nx_out*dx_out),int(nx_out))

u_field_temp = np.zeros([nz, ny_in, nx_out+1])
v_field_temp = np.zeros([nz, ny_in+1, nx_out])
t_temp = np.zeros([nz, ny_in, nx_out])
ph_temp = np.zeros([nz+1, ny_in, nx_out])
phb_temp = np.zeros([nz+1, ny_in, nx_out])
p_temp = np.zeros([nz, ny_in, nx_out])
pb_temp = np.zeros([nz, ny_in, nx_out])
t_init_temp = np.zeros([nz, ny_in, nx_out])
mu_temp = np.zeros([ny_in, nx_out])
mub_temp = np.zeros([ny_in, nx_out])
moist_temp = np.zeros([nz, ny_in, nx_out])
tsk_temp = np.zeros([ny_in, nx_out])
tmn_temp = np.zeros([ny_in, nx_out])

for j in range(ny_in):
    for i in range(nx_out):
        ph_temp[nz,j,i] = wrf_fields.interp_0(ph_in[nz,j,:], x_in, x_out[i], nx_in_trunc)
        phb_temp[nz,j,i] = wrf_fields.interp_0(phb_in[nz,j,:], x_in, x_out[i], nx_in_trunc)
        for k in range(nz):
            u_field_temp[k,j,i] = wrf_fields.interp_0(u_field_in[k,j,:], x_in, x_out[i], nx_in_trunc)
            v_field_temp[k,j,i] = wrf_fields.interp_0(v_field_in[k,j,:], x_in, x_out[i], nx_in_trunc)
            t_temp[k,j,i] = wrf_fields.interp_0(t_in[k,j,:], x_in, x_out[i], nx_in_trunc)
            ph_temp[k,j,i] = wrf_fields.interp_0(ph_in[k,j,:], x_in, x_out[i], nx_in_trunc)
            phb_temp[k,j,i] = wrf_fields.interp_0(phb_in[k,j,:], x_in, x_out[i], nx_in_trunc)
            p_temp[k,j,i] = wrf_fields.interp_0(p_in[k,j,:], x_in, x_out[i], nx_in_trunc)
            pb_temp[k,j,i] = wrf_fields.interp_0(pb_in[k,j,:], x_in, x_out[i], nx_in_trunc)
            t_init_temp[k,j,i] = wrf_fields.interp_0(t_init_in[k,j,:], x_in, x_out[i], nx_in_trunc)
            mu_temp[j,i] = wrf_fields.interp_0(mu_in[j,:], x_in, x_out[i], nx_in_trunc)
            mub_temp[j,i] = wrf_fields.interp_0(mub_in[j,:], x_in, x_out[i], nx_in_trunc)
            moist_temp[k,j,i] = wrf_fields.interp_0(moist_in[k,j,:], x_in, x_out[i], nx_in_trunc)
            tsk_temp[j,i] = wrf_fields.interp_0(tsk_in[j,:], x_in, x_out[i], nx_in_trunc)
            tmn_temp[j,i] = wrf_fields.interp_0(tmn_in[j,:], x_in, x_out[i], nx_in_trunc)
            
u_field_temp[:,:,-1] = u_field_temp[:,:,-2]
v_field_temp[:,-1,:] = v_field_temp[:,-1,:]

y_in = np.linspace(int(dx_in),int(ny_in*dx_in),int(ny_in))
y_out = np.linspace(int(dx_out),int(ny_out*dx_out),int(ny_out))

u_field_out = np.zeros([nz, ny_out, nx_out+1])
v_field_out = np.zeros([nz, ny_out+1, nx_out])
t_out = np.zeros([nz, ny_out, nx_out])
ph_out = np.zeros([nz+1, ny_out, nx_out])
phb_out = np.zeros([nz+1, ny_out, nx_out])
p_out = np.zeros([nz, ny_out, nx_out])
pb_out = np.zeros([nz, ny_out, nx_out])
t_init_out = np.zeros([nz, ny_out, nx_out])
mu_out = np.zeros([ny_out, nx_out])
mub_out = np.zeros([ny_out, nx_out])
moist_out = np.zeros([nz, ny_out, nx_out])
tsk_out = np.zeros([ny_out, nx_out])
tmn_out = np.zeros([ny_out, nx_out])

for j in range(ny_out):
    for i in range(nx_out):
        ph_out[nz,j,i] = wrf_fields.interp_0(ph_temp[nz,:,i], y_in, y_out[j], ny_in)
        phb_out[nz,j,i] = wrf_fields.interp_0(phb_temp[nz,:,i], y_in, y_out[j], ny_in)
        for k in range(nz):
            u_field_out[k,j,i] = wrf_fields.interp_0(u_field_temp[k,:,i], y_in, y_out[j], ny_in)
            v_field_out[k,j,i] = wrf_fields.interp_0(v_field_temp[k,:,i], y_in, y_out[j], ny_in)
            t_out[k,j,i] = wrf_fields.interp_0(t_temp[k,:,i], y_in, y_out[j], ny_in)
            ph_out[k,j,i] = wrf_fields.interp_0(ph_temp[k,:,i], y_in, y_out[j], ny_in)
            phb_out[k,j,i] = wrf_fields.interp_0(phb_temp[k,:,i], y_in, y_out[j], ny_in)
            p_out[k,j,i] = wrf_fields.interp_0(p_temp[k,:,i], y_in, y_out[j], ny_in)
            pb_out[k,j,i] = wrf_fields.interp_0(pb_temp[k,:,i], y_in, y_out[j], ny_in)
            t_init_out[k,j,i] = wrf_fields.interp_0(t_init_temp[k,:,i], y_in, x_out[j], ny_in)
            mu_out[j,i] = wrf_fields.interp_0(mu_temp[:,i], y_in, y_out[j], ny_in)
            mub_out[j,i] = wrf_fields.interp_0(mub_temp[:,i], y_in, y_out[j], ny_in)
            moist_out[k,j,i] = wrf_fields.interp_0(moist_temp[k,:,i], y_in, y_out[j], ny_in)
            tsk_out[j,i] = wrf_fields.interp_0(tsk_temp[:,i], y_in, y_out[j], ny_in)
            tmn_out[j,i] = wrf_fields.interp_0(tmn_temp[:,i], y_in, y_out[j], ny_in)
            
u_field_out[:,:,-1] = u_field_out[:,:,-2]
v_field_out[:,-1,:] = v_field_out[:,-1,:]

# WRITE

nc_out = nc.Dataset('/p/work1/lloveras/48h_adj_30km/in_files/wrfin_ctl_4km','w',format='NETCDF3_64BIT_OFFSET')
nc_out = write_wrfinputd01.write(nc_out, nx_out, ny_out, nz, hres_m_out, title_str, \
        time_str, u_field_out, v_field_out, t_out, ph_out, phb_out, t_init_out, mu_out, mub_out, p_out, pb_out, \
        fnm, fnp, rdnw, rdn, dnw, dn, t_base, cfn, cfn1, rdx_out, rdy_out, cf1, \
        cf2, cf3, tsk_out, u_base, v_base, qv_base, tmn_out, moist_out, znw, znu, \
        diff_opt, km_opt, damp_opt, dampcoef, khdif, kvdif, mp_physics, \
        ra_lw_physics, ra_sw_physics, sf_sfclay_physics, sf_surface_physics, \
        bl_pbl_physics, cu_physics, sf_lake_physics, surface_input_source, \
        hypsometric_opt, dt, num_land_cat, num_soil_layers, num_soil_cat, \
        spec_bdy_width, ptop)
nc_out.close()
