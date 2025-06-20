#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:54:47 2020

@author: lloverasdan
"""

import numpy as np
import xarray as xr
from numba import njit

@njit
def low_pert(nx, ny, nz, dx, dz, a, h, l):
    
    pert = np.zeros((nz,ny,nx))
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                pert[z,y,x] = a*np.exp(-z*dz/h)*np.sin(2*np.pi*x*dx/l)*np.sin(2*np.pi*y*dx/l)
                
    return pert

@njit
def gaussian_pert(nx, ny, nz, dx, dz, z0, a, h, l):
    
    pert = np.zeros((nz,ny,nx))
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                pert[z,y,x] = a*np.exp(-(z*dz - z0)**2/(2*h**2))*np.sin(2*np.pi*x*dx/l)*np.sin(2*np.pi*y*dx/l)
                
    return pert

def pert_wrfrst(file_in, file_out, dx, dz, a, h, l, pert_type='low', z0=10):
    
    data = xr.open_dataset(file_in,chunks={'south_north':10,'Time':1})
    T1 = data.get('T_1')
    T2 = data.get('T_2')
    x = np.asarray(T1)
    nz = x.shape[1]
    ny = x.shape[2]
    nx = x.shape[3]
    if pert_type=='low':
        pert = low_pert(nx,ny,nz,dx,dz,a,h,l)
    elif pert_type=='gaussian':
        pert = gaussian_pert(nx,ny,nz,dx,dz,z0,a,h,l)
    T1 = T1+pert
    T2 = T2+pert
    data.T_1.values = T1
    data.T_2.values = T2
    out = data.to_netcdf(file_out)
    
    return out
