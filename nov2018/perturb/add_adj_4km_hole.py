#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  18 09:18:04 2022

@author: lloverasdan
"""

### Import libraries

import xarray as xr
import numpy as np
from netCDF4 import Dataset
from wrf import interplevel, getvar, destagger
from scipy.interpolate import griddata

### Constants

CP = 1005.7
RD = 287.04
P0 = 1000.
TR = 300.
LV = 2.501e6
EPS = 1.

### COAMPS parameters

nx_c = 301
ny_c = 201
nz_c = 45
coamps_dir = '/p/work2/doyle/coamps/data/cycnov_pverror/'

### WRF parameters

nx_w = 3000
ny_w = 1875
nz_w = 59
pert_dir = '/p/work1/lloveras/real/nov2018/4km_files/adj_hole/wrfrst_d01_2018-11-13_12_00_00'
ctl_dir = '/p/work1/lloveras/real/nov2018/4km_files/ctl/wrfout_d01_2018-11-13_12_00_00'
ncfile_ctl = Dataset(ctl_dir)

### Interpolation settings

rescale_unity = False

largescale_filter = False
lam = 1000.
dx = 4.

latlon_box = True
min_lat = 25
max_lat = 40
min_lon = -107
max_lon = -92

### Read in the COAMPS data

### Latitude and longitude
lats_c = np.fromfile(coamps_dir + 'latitu_sfc_000000_000000_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4')
lons_c = np.fromfile(coamps_dir + 'longit_sfc_000000_000000_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4')
lats_c = np.reshape(lats_c, [ny_c, nx_c])
lons_c = np.reshape(lons_c, [ny_c, nx_c])
lons_c[lons_c > 0] -= 360

### Full pressure in hPa
p_c = (np.fromfile(coamps_dir + '9x/ttlprs_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4') \
      + np.fromfile(coamps_dir + 'ttlprs_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4'))/2
p_c = np.flip(np.reshape(p_c, [nz_c, ny_c, nx_c]), axis=0)

### Theta
t_c = np.fromfile(coamps_dir + '9x/pottmp_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4') \
        - np.fromfile(coamps_dir + 'pottmp_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4')
t_c = np.flip(np.reshape(t_c, [nz_c, ny_c, nx_c]),axis=0)

### Wind
u_c = np.fromfile(coamps_dir + '9x/uuwind_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4') \
        - np.fromfile(coamps_dir + 'uuwind_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4')
v_c = np.fromfile(coamps_dir + '9x/vvwind_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4') \
        - np.fromfile(coamps_dir + 'vvwind_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4')
w_c = np.fromfile(coamps_dir + '9x/wwwind_sig_029735_000000_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4') \
        - np.fromfile(coamps_dir + 'wwwind_sig_029735_000000_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4')
u_c = np.flip(np.reshape(u_c, [nz_c, ny_c, nx_c]),axis=0)
v_c = np.flip(np.reshape(v_c, [nz_c, ny_c, nx_c]),axis=0)
w_c = np.flip(np.reshape(w_c, [nz_c+1, ny_c, nx_c]),axis=0)
w_c_ds = destagger(w_c,stagger_dim=0)

### Water vapor
qv_c = np.fromfile(coamps_dir + '9x/wvapor_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4') \
        - np.fromfile(coamps_dir + 'wvapor_sig_028485_000010_1a0301x0201_2018111312_00000000_fcstfld', dtype='>f4')
qv_c = np.flip(np.reshape(qv_c, [nz_c, ny_c, nx_c]),axis=0)

### Read in the WRF data

### Latitude and longitude
lats_w = np.asarray(getvar(ncfile_ctl,'lat',timeidx=0))
lons_w = np.asarray(getvar(ncfile_ctl,'lon',timeidx=0))
lons_w[lons_w > 0] -= 360

### Full pressure in hPa
p_w = np.asarray(getvar(ncfile_ctl,'pressure',timeidx=0))

### Heights in m
z_w = np.asarray(getvar(ncfile_ctl,'z',timeidx=0))

### INTERPOLATION STARTS HERE

### Start by interpolating WRF pressure field onto COAMPS horizontal grid

p_w_cxy = np.zeros((nz_w,ny_c,nx_c))
for k in range(nz_w):
    p_w_cxy[k,:,:] = griddata((lons_w.ravel(), lats_w.ravel()),
                              p_w[k,:,:].ravel(), (lons_c, lats_c), method='linear')

### Now interpolate COAMPS perturbations onto WRF pressure levels

t_c_wz = np.zeros([nz_w, ny_c, nx_c])
u_c_wz = np.zeros((nz_w, ny_c, nx_c))
v_c_wz = np.zeros((nz_w, ny_c, nx_c))
w_c_wz = np.zeros((nz_w, ny_c, nx_c))
qv_c_wz = np.zeros((nz_w, ny_c, nx_c))
for k in range(nz_w):
    t_c_wz[k,:,:] = np.nan_to_num(interplevel(t_c, p_c, p_w_cxy[k,:,:]))
    u_c_wz[k,:,:] = np.nan_to_num(interplevel(u_c, p_c, p_w_cxy[k,:,:]))
    v_c_wz[k,:,:] = np.nan_to_num(interplevel(v_c, p_c, p_w_cxy[k,:,:]))
    w_c_wz[k,:,:] = np.nan_to_num(interplevel(w_c_ds, p_c, p_w_cxy[k,:,:]))
    qv_c_wz[k,:,:] = np.nan_to_num(interplevel(qv_c, p_c, p_w_cxy[k,:,:]))

### Finally interpolate perturbations onto WRF horizontal grid

t_c_w = np.zeros((nz_w, ny_w, nx_w))
u_c_w = np.zeros((nz_w, ny_w, nx_w))
v_c_w = np.zeros((nz_w, ny_w, nx_w))
w_c_w = np.zeros((nz_w, ny_w, nx_w))
qv_c_w = np.zeros((nz_w, ny_w, nx_w))
for k in range(nz_w):
    t_c_w[k,:,:] = np.nan_to_num(griddata((lons_c.ravel(), lats_c.ravel()), 
                                          t_c_wz[k,:,:].ravel(), (lons_w, lats_w), method='linear'))
    u_c_w[k,:,:] = np.nan_to_num(griddata((lons_c.ravel(), lats_c.ravel()), 
                                          u_c_wz[k,:,:].ravel(), (lons_w, lats_w), method='linear'))
    v_c_w[k,:,:] = np.nan_to_num(griddata((lons_c.ravel(), lats_c.ravel()), 
                                          v_c_wz[k,:,:].ravel(), (lons_w, lats_w), method='linear'))
    w_c_w[k,:,:] = np.nan_to_num(griddata((lons_c.ravel(), lats_c.ravel()), 
                                          w_c_wz[k,:,:].ravel(), (lons_w, lats_w), method='linear'))
    qv_c_w[k,:,:] = np.nan_to_num(griddata((lons_c.ravel(), lats_c.ravel()), 
                                          qv_c_wz[k,:,:].ravel(), (lons_w, lats_w), method='linear'))
    
### Rescale perturbations to retain magnitude prior to interpolation

t_max_c = np.amax(np.abs(t_c))
u_max_c = np.amax(np.abs(u_c))
v_max_c = np.amax(np.abs(v_c))
w_max_c = np.amax(np.abs(w_c))
qv_max_c = np.amax(np.abs(qv_c))
 
t_max_w = np.amax(np.abs(t_c_w))
u_max_w = np.amax(np.abs(u_c_w))
v_max_w = np.amax(np.abs(v_c_w))
w_max_w = np.amax(np.abs(w_c_w))
qv_max_w = np.amax(np.abs(qv_c_w))    
    
t_r = t_max_c/t_max_w
u_r = u_max_c/u_max_w
v_r = v_max_c/v_max_w
w_r = w_max_c/w_max_w
qv_r = qv_max_c/qv_max_w

t_full = t_c_w*t_r
u_full = u_c_w*u_r
v_full = v_c_w*v_r
w_full = w_c_w*w_r
qv_full = qv_c_w*qv_r

### Option to rescale to unity

if rescale_unity:
    
    s = max(t_max_c,u_max_c,qv_max_c*1000.)
    t_full = t_full/s
    u_full = u_full/s
    v_full = v_full/s
    w_full = w_full/s
    qv_full = qv_full/s

### Option to filter out small scales + rescale to have equal DTE

if largescale_filter:
    
    kx = np.fft.fftfreq(nx_w,dx)
    ky = np.fft.fftfreq(ny_w,dx)
    peak = 1/lam

    t_fft = np.fft.fft2(t_full)
    t_fft[:,np.abs(ky) > peak, :] = 0.
    t_fft[:,:,np.abs(kx) > peak] = 0.
    t_filt = np.real(np.fft.ifft2(t_fft))

    u_fft = np.fft.fft2(u_full)
    u_fft[:,np.abs(ky) > peak, :] = 0.
    u_fft[:,:,np.abs(kx) > peak] = 0.
    u_filt = np.real(np.fft.ifft2(u_fft))

    v_fft = np.fft.fft2(v_full)
    v_fft[:,np.abs(ky) > peak, :] = 0.
    v_fft[:,:,np.abs(kx) > peak] = 0.
    v_filt = np.real(np.fft.ifft2(v_fft))

    w_fft = np.fft.fft2(w_full)
    w_fft[:,np.abs(ky) > peak, :] = 0.
    w_fft[:,:,np.abs(kx) > peak] = 0.
    w_filt = np.real(np.fft.ifft2(w_fft))

    qv_fft = np.fft.fft2(qv_full)
    qv_fft[:,np.abs(ky) > peak, :] = 0.
    qv_fft[:,:,np.abs(kx) > peak] = 0.
    qv_filt = np.real(np.fft.ifft2(qv_fft))
    
    dte_full = np.sum(u_full**2) + np.sum(v_full**2) + (CP/TR)*np.sum(t_full**2)
    dte_filt = np.sum(u_filt**2) + np.sum(v_filt**2) + (CP/TR)*np.sum(t_filt**2)
    
    r_filt = np.sqrt(dte_full/dte_filt)
    t_filt = t_filt*r_filt
    u_filt = u_filt*r_filt
    v_filt = v_filt*r_filt
    w_filt = w_filt*r_filt
    qv_filt = qv_filt*r_filt

### Option to isolate perturbations to box + rescale to have equal DTE

if latlon_box:
    
    t_trunc = np.copy(t_full)
    u_trunc = np.copy(u_full)
    v_trunc = np.copy(v_full)
    w_trunc = np.copy(w_full)
    qv_trunc = np.copy(qv_full)

    t_trunc[:,(lons_w > min_lon) & (lons_w < max_lon)
               & (lats_w > min_lat) & (lats_w < max_lat)] = np.NaN
    
    u_trunc[:,(lons_w > min_lon) & (lons_w < max_lon)
               & (lats_w > min_lat) & (lats_w < max_lat)] = np.NaN
    
    v_trunc[:,(lons_w > min_lon) & (lons_w < max_lon)
               & (lats_w > min_lat) & (lats_w < max_lat)] = np.NaN
    
    w_trunc[:,(lons_w > min_lon) & (lons_w < max_lon)
               & (lats_w > min_lat) & (lats_w < max_lat)] = np.NaN
    
    qv_trunc[:,(lons_w > min_lon) & (lons_w < max_lon)
               & (lats_w > min_lat) & (lats_w < max_lat)] = np.NaN

    t_trunc = np.nan_to_num(t_trunc)
    u_trunc = np.nan_to_num(u_trunc)
    v_trunc = np.nan_to_num(v_trunc)
    w_trunc = np.nan_to_num(w_trunc)
    qv_trunc = np.nan_to_num(qv_trunc)
    
    dte_full = np.sum(u_full**2) + np.sum(v_full**2) + (CP/TR)*np.sum(t_full**2)
    dte_trunc = np.sum(u_trunc**2) + np.sum(v_trunc**2) + (CP/TR)*np.sum(t_trunc**2)
    
    r_trunc = np.sqrt(dte_full/dte_trunc)
    t_trunc = t_trunc*r_trunc
    u_trunc = u_trunc*r_trunc
    v_trunc = v_trunc*r_trunc
    w_trunc = w_trunc*r_trunc
    qv_trunc = qv_trunc*r_trunc

### Collect perturbations    
    
t_pert = t_full
u_pert = u_full
v_pert = v_full
w_pert = w_full
qv_pert = qv_full

if largescale_filter:
    
    t_pert = t_filt
    u_pert = u_filt
    v_pert = v_filt
    w_pert = w_filt
    qv_pert = qv_filt

if latlon_box:
    
    t_pert = t_trunc
    u_pert = u_trunc
    v_pert = v_trunc
    w_pert = w_trunc
    qv_pert = qv_trunc

### FINALLY ADD PERTURBATIONS TO FILE HERE

### Open the netCDF file for writing
ncfile_pert = Dataset(pert_dir,'r+')

### Theta
t1 = np.asarray(ncfile_pert.variables['THM_1'][0,:,:,:])
ncfile_pert.variables['THM_1'][0,:,:,:] = t1 + t_pert

t2 = np.asarray(ncfile_pert.variables['THM_2'][0,:,:,:])
ncfile_pert.variables['THM_2'][0,:,:,:] = t2 + t_pert

### Zonal wind
u1 = np.asarray(ncfile_pert.variables['U_1'][0,:,:,:])
ncfile_pert.variables['U_1'][0,:,:,:-1] = u1[:,:,:-1] + u_pert
ncfile_pert.variables['U_1'][0,:,:,-1] = u1[:,:,-1]

u2 = np.asarray(ncfile_pert.variables['U_2'][0,:,:,:])
ncfile_pert.variables['U_2'][0,:,:,:-1] = u2[:,:,:-1] + u_pert
ncfile_pert.variables['U_2'][0,:,:,-1] = u2[:,:,-1]

### Meridional wind

v1 = np.asarray(ncfile_pert.variables['V_1'][0,:,:,:])
ncfile_pert.variables['V_1'][0,:,:-1,:] = v1[:,:-1,:] + v_pert
ncfile_pert.variables['V_1'][0,:,-1,:] = v1[:,-1,:]

v2 = np.asarray(ncfile_pert.variables['V_2'][0,:,:,:])
ncfile_pert.variables['V_2'][0,:,:-1,:] = v2[:,:-1,:] + v_pert
ncfile_pert.variables['V_2'][0,:,-1,:] = v2[:,-1,:]

### Vertical wind

w1 = np.asarray(ncfile_pert.variables['W_1'][0,:,:,:])
ncfile_pert.variables['W_1'][0,:-1,:,:] = w1[:-1,:,:] + w_pert
ncfile_pert.variables['W_1'][0,-1,:,:] = w1[-1,:,:]

w2 = np.asarray(ncfile_pert.variables['W_2'][0,:,:,:])
ncfile_pert.variables['W_2'][0,:-1,:,:] = w2[:-1,:,:] + w_pert
ncfile_pert.variables['W_2'][0,-1,:,:] = w2[-1,:,:]

### Water vapor
qv = np.asarray(ncfile_pert.variables['QVAPOR'][0,:,:,:])
ncfile_pert.variables['QVAPOR'][0,:,:,:] = qv + qv_pert

### Close the perturbed netCDF file
ncfile_pert.close()
