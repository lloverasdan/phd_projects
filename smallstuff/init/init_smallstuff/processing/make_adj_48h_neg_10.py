#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  18 09:18:04 2022

@author: lloverasdan
"""

import numpy as np
import netCDF4 as nc
from wrf import getvar

ctl_file = '/p/work1/lloveras/adj_4km/rst_input/ctl_48h/wrfrst_d01_2021-01-03_00_00_00'
adj_file = '/p/work1/lloveras/adj_4km/rst_input/adj_48h/wrfrst_d01_2021-01-03_00_00_00'
new_dir = '/p/work1/lloveras/adj_4km/rst_input/adj_48h_neg_10/wrfrst_d01_2021-01-03_00_00_00'

nc_ctl = nc.Dataset(ctl_file)
nc_adj = nc.Dataset(adj_file)

u_ctl = np.asarray(getvar(nc_ctl,'U_2'))
v_ctl = np.asarray(getvar(nc_ctl,'V_2'))
t_ctl = np.asarray(getvar(nc_ctl,'T_2'))
p_ctl = np.asarray(getvar(nc_ctl,'P'))
qv_ctl = np.asarray(getvar(nc_ctl,'QVAPOR'))

u_adj = np.asarray(getvar(nc_adj,'U_2'))
v_adj = np.asarray(getvar(nc_adj,'V_2'))
t_adj = np.asarray(getvar(nc_adj,'T_2'))
p_adj = np.asarray(getvar(nc_adj,'P'))
qv_adj = np.asarray(getvar(nc_adj,'QVAPOR'))

u_pert = u_adj - u_ctl
v_pert = v_adj - v_ctl
t_pert = t_adj - t_ctl
p_pert = p_adj - p_ctl
qv_pert = qv_adj - qv_ctl

ncfile = nc.Dataset(new_dir,'r+')
ncfile.variables['U_1'][0,:,:,:] = u_ctl - u_pert/10.
ncfile.variables['V_1'][0,:,:,:] = v_ctl - v_pert/10.
ncfile.variables['T_1'][0,:,:,:] = t_ctl - t_pert/10.
ncfile.variables['P'][0,:,:,:] = p_ctl - p_pert/10.
ncfile.variables['QVAPOR'][0,:,:,:] = qv_ctl - qv_pert/10.
ncfile.variables['U_2'][0,:,:,:] = u_ctl - u_pert/10.
ncfile.variables['V_2'][0,:,:,:] = v_ctl - v_pert/10.
ncfile.variables['T_2'][0,:,:,:] = t_ctl - t_pert/10.

ncfile.close()
