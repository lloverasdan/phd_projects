#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

from netCDF4 import Dataset
import numpy as np
from wrf import getvar, vinterp

### Time steps
ti = np.asarray([1,2,3,4,5,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,54,60,66,72,78,84,90,96])

### Grid parameters
nx = 2000
ny = 1800
dx = 4.0e3
lam = 1000.0e3

### Constants
F0 = 1e-4
G = 9.81

### Load the netCDF files
nc_ctl_0 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/short_runs/ctl_48h/wrfout_d01_2021-01-03_00_00_00')
nc_pert_0 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/short_runs/adj_48h_100/wrfout_d01_2021-01-03_00_00_00')

nc_ctl_48 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_48h/wrfout_d01_2021-01-03_00_00_00')
nc_pert_48 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/adj_48h_100/wrfout_d01_2021-01-03_00_00_00')

nc_ctl_96 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_48h/wrfout_d01_2021-01-05_00_00_00')
nc_pert_96 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/adj_48h_100/wrfout_d01_2021-01-05_00_00_00')

### Output file path
file_path = '/p/work1/lloveras/adj_4km/processed/balance/bal_adj_48h_100'

### Define filters
def horizontal_filter(var, nx, ny, dx, lam):
    
    var = np.append(var,np.flip(var,axis=1),axis=1)
    
    kx = np.fft.fftfreq(nx,dx)
    ky = np.fft.fftfreq(int(ny*2),dx)
    peak = 1/lam
    
    var_fft = np.fft.fft2(var)
    var_fft[:,np.abs(ky) > peak, :] = 0.
    var_fft[:,:,np.abs(kx) > peak] = 0.
    var_filt = np.real(np.fft.ifft2(var_fft))
    
    var_filt = var_filt[:,:ny,:]
    
    return var_filt

def horizontal_filter_v(var, nx, ny, dx, lam):
    
    var = np.append(var,-np.flip(var,axis=1),axis=1)
    
    kx = np.fft.fftfreq(nx,dx)
    ky = np.fft.fftfreq(int(ny*2),dx)
    peak = 1/lam
    
    var_fft = np.fft.fft2(var)
    var_fft[:,np.abs(ky) > peak, :] = 0.
    var_fft[:,:,np.abs(kx) > peak] = 0.
    var_filt = np.real(np.fft.ifft2(var_fft))
    
    var_filt = var_filt[:,:ny,:]
    
    return var_filt

### Define computation of balance number
def bal(u_ctl, v_ctl, ph_ctl, u_pert, v_pert, ph_pert, nx, ny, dx, lam):
    
    du = horizontal_filter(np.nan_to_num(u_pert - u_ctl), nx, ny, dx, lam)
    dv = horizontal_filter_v(np.nan_to_num(v_pert - v_ctl), nx, ny, dx, lam)
    dph = horizontal_filter(np.nan_to_num(ph_pert - ph_ctl), nx, ny, dx, lam)

    dph_dy, dph_dx = np.gradient(dph, dx, axis=(1,2))
    
    du_g = -dph_dy/F0
    dv_g = dph_dx/F0
    
    du_ag = du - du_g
    dv_ag = dv - dv_g
    
    num = np.sum(np.sqrt(du_ag**2 + dv_ag**2))/np.sum(np.sqrt(du_g**2 + dv_g**2))
    
    return num

### Read in data and compute
bal_pert = np.zeros(len(ti))
levs = np.arange(950,40,-10)
for i in range(len(ti)):
    if ti[i] <= 6:
        rol = 0
        
        u_ctl = np.roll(np.asarray(vinterp(nc_ctl_0,getvar(nc_ctl_0,'ua',timeidx=int(ti[i]*2)),\
                        'pressure',levs,timeidx=int(ti[i]*2))),rol,axis=-1)
        u_pert = np.roll(np.asarray(vinterp(nc_pert_0,getvar(nc_pert_0,'ua',timeidx=int(ti[i]*2)),\
                        'pressure',levs,timeidx=int(ti[i]*2))),rol,axis=-1)
        
        v_ctl = np.roll(np.asarray(vinterp(nc_ctl_0,getvar(nc_ctl_0,'va',timeidx=int(ti[i]*2)),\
                        'pressure',levs,timeidx=int(ti[i]*2))),rol,axis=-1)
        v_pert = np.roll(np.asarray(vinterp(nc_pert_0,getvar(nc_pert_0,'va',timeidx=int(ti[i]*2)),\
                        'pressure',levs,timeidx=int(ti[i]*2))),rol,axis=-1)
        
        ph_ctl = np.roll(np.asarray(vinterp(nc_ctl_0,getvar(nc_ctl_0,'geopt',timeidx=int(ti[i]*2)),\
                        'pressure',levs,timeidx=int(ti[i]*2))),rol,axis=-1)
        ph_pert = np.roll(np.asarray(vinterp(nc_pert_0,getvar(nc_pert_0,'geopt',timeidx=int(ti[i]*2)),\
                        'pressure',levs,timeidx=int(ti[i]*2))),rol,axis=-1)
        
        bal_pert[i] = bal(u_ctl,v_ctl,ph_ctl,u_pert,v_pert,ph_pert,nx,ny,dx,lam)
    
    elif ti[i] > 6 and ti[i] <= 48:
        rol = -1000
        
        u_ctl = np.roll(np.asarray(vinterp(nc_ctl_48,getvar(nc_ctl_48,'ua',timeidx=int(ti[i]/3)),\
                        'pressure',levs,timeidx=int(ti[i]/3))),rol,axis=-1)
        u_pert = np.roll(np.asarray(vinterp(nc_pert_48,getvar(nc_pert_48,'ua',timeidx=int(ti[i]/3)),\
                        'pressure',levs,timeidx=int(ti[i]/3))),rol,axis=-1)
        
        v_ctl = np.roll(np.asarray(vinterp(nc_ctl_48,getvar(nc_ctl_48,'va',timeidx=int(ti[i]/3)),\
                        'pressure',levs,timeidx=int(ti[i]/3))),rol,axis=-1)
        v_pert = np.roll(np.asarray(vinterp(nc_pert_48,getvar(nc_pert_48,'va',timeidx=int(ti[i]/3)),\
                        'pressure',levs,timeidx=int(ti[i]/3))),rol,axis=-1)
        
        ph_ctl = np.roll(np.asarray(vinterp(nc_ctl_48,getvar(nc_ctl_48,'geopt',timeidx=int(ti[i]/3)),\
                        'pressure',levs,timeidx=int(ti[i]/3))),rol,axis=-1)
        ph_pert = np.roll(np.asarray(vinterp(nc_pert_48,getvar(nc_pert_48,'geopt',timeidx=int(ti[i]/3)),\
                        'pressure',levs,timeidx=int(ti[i]/3))),rol,axis=-1)
        
        bal_pert[i] = bal(u_ctl,v_ctl,ph_ctl,u_pert,v_pert,ph_pert,nx,ny,dx,lam)
        
    else:
        rol = -1500
        
        u_ctl = np.roll(np.asarray(vinterp(nc_ctl_96,getvar(nc_ctl_96,'ua',timeidx=int(ti[i]/6 - 8)),\
                        'pressure',levs,timeidx=int(ti[i]/6 - 8))),rol,axis=-1)
        u_pert = np.roll(np.asarray(vinterp(nc_pert_96,getvar(nc_pert_96,'ua',timeidx=int(ti[i]/6 - 8)),\
                        'pressure',levs,timeidx=int(ti[i]/6 - 8))),rol,axis=-1)
        
        v_ctl = np.roll(np.asarray(vinterp(nc_ctl_96,getvar(nc_ctl_96,'va',timeidx=int(ti[i]/6 - 8)),\
                        'pressure',levs,timeidx=int(ti[i]/6 - 8))),rol,axis=-1)
        v_pert = np.roll(np.asarray(vinterp(nc_pert_96,getvar(nc_pert_96,'va',timeidx=int(ti[i]/6 - 8)),\
                        'pressure',levs,timeidx=int(ti[i]/6 - 8))),rol,axis=-1)      
        
        ph_ctl = np.roll(np.asarray(vinterp(nc_ctl_96,getvar(nc_ctl_96,'geopt',timeidx=int(ti[i]/6 - 8)),\
                        'pressure',levs,timeidx=int(ti[i]/6 - 8))),rol,axis=-1)
        ph_pert = np.roll(np.asarray(vinterp(nc_pert_96,getvar(nc_pert_96,'geopt',timeidx=int(ti[i]/6 - 8)),\
                        'pressure',levs,timeidx=int(ti[i]/6 - 8))),rol,axis=-1)
        
        bal_pert[i] = bal(u_ctl,v_ctl,ph_ctl,u_pert,v_pert,ph_pert,nx,ny,dx,lam)
    
### Save the output
np.save(file_path,bal_pert)
