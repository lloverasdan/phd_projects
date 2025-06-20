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
from bwave_ideal_wrf import wrf_fields
from wrf import getvar, destagger, interplevel

CP = 1005.7
RD = 287.04
P0 = 1000.
TR = 300.
LV = 2.501e6
EPS = 1.
G = 9.81

nx = 2000
ny = 1800
nz = 100
dx = 4000.
lxy = 28000.
dz = 200.
h = 3000.
xnot = int(nx/2)
rol = 975
ynot = 925
ax = 350000.
ay = 450000.

ctl_dir = '/p/work1/lloveras/adj_4km/input_data/rst_files/ctl_96h/wrfrst_d01_2021-01-05_00_00_00'
adj_dir = '/p/work1/lloveras/adj_4km/input_data/rst_files/adj_96h_100/wrfrst_d01_2021-01-05_00_00_00'
wave_dir = '/p/work1/lloveras/adj_4km/input_data/rst_files/wave_96h_100/wrfrst_d01_2021-01-05_00_00_00'

ncfile_ctl = nc.Dataset(ctl_dir)
ncfile_adj = nc.Dataset(adj_dir)

u_ctl = getvar(ncfile_ctl,'U_2')
v_ctl = getvar(ncfile_ctl,'V_2')
th_ctl = getvar(ncfile_ctl,'T_2')
p_ctl = getvar(ncfile_ctl,'P')
pb_ctl = getvar(ncfile_ctl,'PB')
qv_ctl = getvar(ncfile_ctl,'QVAPOR')
full_th_ctl = np.asarray(th_ctl) + TR
full_p_ctl = (np.asarray(p_ctl) + np.asarray(pb_ctl))/100.
full_tk_ctl = full_th_ctl*(full_p_ctl/P0)**(RD/CP)
z_ctl = (destagger(np.asarray(getvar(ncfile_ctl,'PH_2')) + np.asarray(getvar(ncfile_ctl,'PHB')),-3))/G

u_adj = getvar(ncfile_adj,'U_2')
v_adj = getvar(ncfile_adj,'V_2')
th_adj = getvar(ncfile_adj,'T_2')
p_adj = getvar(ncfile_adj,'P')
pb_adj = getvar(ncfile_adj,'PB')
qv_adj = getvar(ncfile_adj,'QVAPOR')
full_th_adj = np.asarray(th_adj) + TR
full_p_adj = (np.asarray(p_adj) + np.asarray(pb_adj))/100.
full_tk_adj = full_th_adj*(full_p_adj/P0)**(RD/CP) 

du_adj = np.asarray(u_adj) - np.asarray(u_ctl)
dv_adj = np.asarray(v_adj) - np.asarray(v_ctl)
dtk_adj = full_tk_adj - full_tk_ctl

dte1_adj = np.sum(np.sum(np.sum(np.squeeze(du_adj[:,:,:]**2.0),1),1))
dte2_adj = np.sum(np.sum(np.sum(np.squeeze(dv_adj[:,:,:]**2.0),1),1))
dte3_adj = np.sum(np.sum(np.sum(np.squeeze(CP/TR*(dtk_adj[:,:,:]**2.0)),1),1))
dte_adj = 0.5*(dte1_adj + dte2_adj + dte3_adj)

dqv_adj = np.asarray(qv_adj) - np.asarray(qv_ctl)
dth_adj = full_th_adj - full_th_ctl
dqv_max = np.amax(np.abs(dqv_adj))
dth_max = np.amax(np.abs(dth_adj))
aqv = dqv_max/dth_max

### Wave perturbation
ncfile_wave = nc.Dataset(wave_dir,'r+')
th_pert_wave = np.zeros((nz,ny,nx))
qv_pert_wave = np.zeros((nz,ny,nx))
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            r = np.sqrt(((i - xnot)*dx/ax)**2 + ((j - ynot)*dx/ay)**2)
            if r > np.pi/2: r = np.pi/2.
            th_pert_wave[k,j,i] = np.cos(np.pi*k*dz/(2*h))*np.cos(-r)*np.sin(2*np.pi*i*dx/lxy)*np.sin(2*np.pi*j*dx/lxy)
            qv_pert_wave[k,j,i] = aqv*np.cos(np.pi*k*dz/(2*h))*np.cos(-r)*np.sin(2*np.pi*i*dx/lxy)*np.sin(2*np.pi*j*dx/lxy)
            if k*dz > h:
                th_pert_wave[k,j,i] = 0.
                qv_pert_wave[k,j,i] = 0.
                break

### Balance hydrostatically
p_pert_wave = np.zeros((nz,ny,nx))
for k in range(nz-1,0,-1):
    p_pert_wave[k-1,:,:] = p_pert_wave[k,:,:] - (z_ctl[k,:,:] - z_ctl[k-1,:,:])*(th_pert_wave[k,:,:] + th_pert_wave[k-1,:,:])*G/(2*TR)

### Roll
p_pert_wave = np.roll(p_pert_wave, rol, axis=-1)
th_pert_wave = np.roll(th_pert_wave, rol, axis=-1)
qv_pert_wave = np.roll(qv_pert_wave, rol, axis=-1)

### Make sure DTE is the same
full_th_wave = full_th_ctl + th_pert_wave
full_p_wave = full_p_ctl + p_pert_wave/100.
full_tk_wave = full_th_wave*(full_p_wave/P0)**(RD/CP) 
dtk_wave = full_tk_wave - full_tk_ctl

dte3_wave = np.sum(np.sum(np.sum(np.squeeze(CP/TR*(dtk_wave[:,:,:]**2.0)),1),1))
dte_wave = 0.5*(dte3_wave)

fac = dte_adj/dte_wave
sca = np.sqrt(fac)
p_pert_wave = sca*p_pert_wave
th_pert_wave = sca*th_pert_wave
qv_pert_wave = sca*qv_pert_wave

ncfile_wave.variables['P'][0,:,:,:] = p_ctl + p_pert_wave
ncfile_wave.variables['T_1'][0,:,:,:] = th_ctl + th_pert_wave
ncfile_wave.variables['T_2'][0,:,:,:] = th_ctl + th_pert_wave
ncfile_wave.variables['QVAPOR'][0,:,:,:] = qv_ctl + qv_pert_wave

ncfile_wave.close()
