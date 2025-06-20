#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:54:47 2020

@author: lloverasdan
"""

import numpy as np
from wrf import getvar, vinterp
from scipy.fftpack import fft, fft2
from numba import njit

TR = 300.
P0 = 1000.
LV = 2.501e6
EPS = 1.
CP = 1005.7
RD = 287.04

#-------------------- DTE --------------------#


def dte(base_ncfile, pert_ncfile, in_ncfile, timeid=0):
    """
    Computes the domain-integrated difference total energy at a given time
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param pert_ncfile: netCDF Dataset: perturbed WRF output file opened using netCDF4.Dataset()
    :param time: int: time index (default 0)
    Returns a float
    """
    u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=timeid))
    v_base = np.asarray(getvar(base_ncfile,'va',timeidx=timeid))
    tk_base = np.asarray(getvar(base_ncfile,'tk',timeidx=timeid))
    u_pert = np.asarray(getvar(pert_ncfile,'ua',timeidx=timeid))
    v_pert = np.asarray(getvar(pert_ncfile,'va',timeidx=timeid))
    tk_pert = np.asarray(getvar(pert_ncfile,'tk',timeidx=timeid))
    du = np.abs(u_pert - u_base)
    dv = np.abs(v_pert - v_base)
    dt = np.abs(tk_pert - tk_base)

    [nz, ny, nx] = np.shape(u_base)
    u_in = np.tile(np.expand_dims(np.asarray(getvar(in_ncfile,'ua'))[:,:,0],-1),nx)
    v_in = np.tile(np.expand_dims(np.asarray(getvar(in_ncfile,'va'))[:,:,0],-1),nx)
    tk_in = np.tile(np.expand_dims(np.asarray(getvar(in_ncfile,'tk'))[:,:,0],-1),nx)
    du_norm = np.abs(u_base - u_in)
    dv_norm = np.abs(u_base - u_in)
    dt_norm = np.abs(tk_base - tk_in)

    du = np.where(du<0.1, np.nan, du)
    dv = np.where(dv<0.1, np.nan, dv)
    dt = np.where(dt<0.1, np.nan, dt)
    
    du_norm = np.where(du_norm<1.0, np.nan, du_norm)
    dv_norm = np.where(dv_norm<1.0, np.nan, dv_norm)
    dt_norm = np.where(dt_norm<1.0, np.nan, dt_norm)
    
    #ratio = (du/du_norm + dv/dv_norm + dt/dt_norm)/3
    ratio = (du/du_norm + dt/dt_norm)/2
    ratio[ratio>1] = 1
    bins = np.linspace(0, 1.0, 11)
    hist, edges = np.histogram(ratio,bins)

    #dems = du_norm + dv_norm + dt_norm
    dems = du_norm + dt_norm
    percent = np.count_nonzero(~np.isnan(ratio))/np.count_nonzero(~np.isnan(dems))
    
    return hist, percent

# def dte(base_ncfile, pert_ncfile, in_ncfile, timeid=0):
#     """
#     Computes the domain-integrated difference total energy at a given time
#     :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
#     :param pert_ncfile: netCDF Dataset: perturbed WRF output file opened using netCDF4.Dataset()
#     :param time: int: time index (default 0)
#     Returns a float
#     """
#     u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=timeid))
#     v_base = np.asarray(getvar(base_ncfile,'va',timeidx=timeid))
#     tk_base = np.asarray(getvar(base_ncfile,'tk',timeidx=timeid))
#     u_pert = np.asarray(getvar(pert_ncfile,'ua',timeidx=timeid))
#     v_pert = np.asarray(getvar(pert_ncfile,'va',timeidx=timeid))
#     tk_pert = np.asarray(getvar(pert_ncfile,'tk',timeidx=timeid))
#     du = u_pert - u_base
#     dv = v_pert - v_base
#     dt = tk_pert - tk_base
#     dte = np.sum(du**2) + np.sum(dv**2) + CP/TR*np.sum(dt**2)

# #     [nz, ny, nx] = np.shape(u_base)
# #     u_in = np.tile(np.expand_dims(np.asarray(getvar(in_ncfile,'ua'))[:,:,0],-1),nx)
# #     v_in = np.tile(np.expand_dims(np.asarray(getvar(in_ncfile,'va'))[:,:,0],-1),nx)
# #     tk_in = np.tile(np.expand_dims(np.asarray(getvar(in_ncfile,'tk'))[:,:,0],-1),nx)
# #     du_norm = u_base - u_in
# #     dv_norm = u_base - u_in
# #     dt_norm = tk_base - tk_in
# #     dte_norm = np.sum(du_norm**2) + np.sum(dv_norm**2) + CP/TR*np.sum(dt_norm**2)
    
# #     dte_tot = 0.5*dte/dte_norm
    
#     return dte

def moist_dte(base_ncfile, pert_ncfile, timeid=0):
    """
    Computes the difference moist total energy at a given time
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param pert_ncfile: netCDF Dataset: perturbed WRF output file opened using netCDF4.Dataset()
    :param time: int: time index (default 0)
    Returns a float
    """
    u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=timeid))
    v_base = np.asarray(getvar(base_ncfile,'va',timeidx=timeid))
    tk_base = np.asarray(getvar(base_ncfile,'tk',timeidx=timeid))
    p_base = np.asarray(getvar(base_ncfile,'p',timeidx=timeid))
    qv_base = np.asarray(getvar(base_ncfile,'QVAPOR',timeidx=timeid))
    u_pert = np.asarray(getvar(pert_ncfile,'ua',timeidx=timeid))
    v_pert = np.asarray(getvar(pert_ncfile,'va',timeidx=timeid))
    tk_pert = np.asarray(getvar(pert_ncfile,'tk',timeidx=timeid))
    p_pert = np.asarray(getvar(pert_ncfile,'p',timeidx=timeid))
    qv_pert = np.asarray(getvar(pert_ncfile,'QVAPOR',timeidx=timeid))
    du = u_base - u_pert
    dv = v_base - v_pert
    dt = tk_base - tk_pert
    dp = p_base - p_pert
    dq = qv_base - qv_pert
    dte1 = np.sum(np.sum(np.sum(np.squeeze(du[:,:,:]**2.0),1),1))
    dte2 = np.sum(np.sum(np.sum(np.squeeze(dv[:,:,:]**2.0),1),1))
    dte3 = np.sum(np.sum(np.sum(np.squeeze(CP/TR*(dt[:,:,:]**2.0)),1),1))
    dte4 = np.sum(np.sum(np.sum(np.squeeze((RD*TR/((P0*100)**2))*(dp[:,:,:]**2.0)),1),1))
    dte5 = np.sum(np.sum(np.sum(np.squeeze(EPS*((LV**2)/(CP*TR))*(dq[:,:,:]**2.0)),1),1))
    dte = 0.5*(dte1 + dte2 + dte3 + dte4 + dte5)
    
    return dte

def dte_time_shift(base_ncfile, pert_ncfile, timeid=0):
    """
    Computes the domain-integrated difference total energy at times in which the output are offset
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param pert_ncfile: netCDF Dataset: perturbed WRF output file opened using netCDF4.Dataset()
    :param time: int: time index (default 0)
    Returns a float
    """
    u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=timeid))
    v_base = np.asarray(getvar(base_ncfile,'va',timeidx=timeid))
    tk_base = np.asarray(getvar(base_ncfile,'tk',timeidx=timeid))
    u_pert = np.asarray(getvar(pert_ncfile,'ua',timeidx=int(timeid*2)))
    v_pert = np.asarray(getvar(pert_ncfile,'va',timeidx=int(timeid*2)))
    tk_pert = np.asarray(getvar(pert_ncfile,'tk',timeidx=int(timeid*2)))
    du = u_base - u_pert
    dv = v_base - v_pert
    dt = tk_base - tk_pert
    dte1 = np.sum(np.sum(np.sum(np.squeeze(du[:,:,:]**2.0),1),1))
    dte2 = np.sum(np.sum(np.sum(np.squeeze(dv[:,:,:]**2.0),1),1))
    dte3 = np.sum(np.sum(np.sum(np.squeeze(CP/TR*(dt[:,:,:]**2.0)),1),1))
    dte = 0.5*(dte1 + dte2 + dte3)
    
    return dte

def dte_phys_shift(base_ncfile, pert_ncfile, timeid=0):
    """
    Computes the domain-integrated difference total energy at a given time in which spatial fields are offset
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param pert_ncfile: netCDF Dataset: perturbed WRF output file opened using netCDF4.Dataset()
    :param time: int: time index (default 0)
    Returns a float
    """
    u_base = np.roll(np.asarray(getvar(base_ncfile,'ua',timeidx=timeid)),-20)
    v_base = np.roll(np.asarray(getvar(base_ncfile,'va',timeidx=timeid)),-20)
    tk_base = np.roll(np.asarray(getvar(base_ncfile,'tk',timeidx=timeid)),-20)
    u_pert = np.asarray(getvar(pert_ncfile,'ua',timeidx=timeid))
    v_pert = np.asarray(getvar(pert_ncfile,'va',timeidx=timeid))
    tk_pert = np.asarray(getvar(pert_ncfile,'tk',timeidx=timeid))
    du = u_base - u_pert
    dv = v_base - v_pert
    dt = tk_base - tk_pert
    dte1 = np.sum(np.sum(np.sum(np.squeeze(du[:,:,:]**2.0),1),1))
    dte2 = np.sum(np.sum(np.sum(np.squeeze(dv[:,:,:]**2.0),1),1))
    dte3 = np.sum(np.sum(np.sum(np.squeeze(CP/TR*(dt[:,:,:]**2.0)),1),1))
    dte = 0.5*(dte1 + dte2 + dte3)
    
    return dte

def te_back_1d(base_ncfile, dx, timeid=0):
    """
    Computes the one-dimensional zonal background total energy spectrum at a given time
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param dx: float: grid spacing in x direction in meters
    :param time: int: time index (default 0)
    Returns a three-dimensional array
    """
    u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=timeid))
    v_base = np.asarray(getvar(base_ncfile,'va',timeidx=timeid))
    tk_base = np.asarray(getvar(base_ncfile,'tk',timeidx=timeid))
    nz = u_base.shape[0]
    ny = u_base.shape[1]
    nx = u_base.shape[2]
    te_spec = np.zeros((nz,ny,nx))
    for k in range(nz):
        for j in range(ny):
            u_fft = np.abs(fft(u_base[k,j,:]))
            v_fft = np.abs(fft(v_base[k,j,:]))
            tk_fft = np.abs(fft(tk_base[k,j,:]))
            for i in range(nx):
                te_spec[k,j,i] = (dx/(2*np.pi*nx))*(np.abs(u_fft[i]*u_fft[i]) + \
                                np.abs(v_fft[i]*v_fft[i]) + \
                                np.abs(CP/TR*tk_fft[i]*tk_fft[i]))
                    
    return te_spec                
    
def dte_1d(base_ncfile, pert_ncfile, dx, timeid=0):
    """
    Computes the one-dimensional zonal difference total energy spectrum at a given time
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param pert_ncfile: netCDF Dataset: perturbed WRF restart file opened using netCDF4.Dataset()
    :param dx: float: grid spacing in x direction in meters
    :param time: int: time index (default 0)
    Returns a three-dimensional array
    """    
    u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=timeid))
    v_base = np.asarray(getvar(base_ncfile,'va',timeidx=timeid))
    tk_base = np.asarray(getvar(base_ncfile,'tk',timeidx=timeid))
    u_pert = np.asarray(getvar(pert_ncfile,'ua',timeidx=timeid))
    v_pert = np.asarray(getvar(pert_ncfile,'va',timeidx=timeid))
    tk_pert = np.asarray(getvar(pert_ncfile,'tk',timeidx=timeid))   
    nz = u_base.shape[0]
    ny = u_base.shape[1]
    nx = u_base.shape[2]
    du = u_base - u_pert
    dv = v_base - v_pert
    dt = tk_base - tk_pert
    dte_spec = np.zeros((nz,ny,nx))
    for k in range(nz):
        for j in range(ny):
            du_fft = np.abs(fft(du[k,j,:]))
            dv_fft = np.abs(fft(dv[k,j,:]))
            dt_fft = np.abs(fft(dt[k,j,:]))
            for i in range(nx):
                dte_spec[k,j,i] = (dx/(2*np.pi*nx))*(np.abs(du_fft[i]*du_fft[i]) + \
                                np.abs(dv_fft[i]*dv_fft[i]) + \
                                np.abs(CP/TR*dt_fft[i]*dt_fft[i]))
                
    return dte_spec

@njit
def te_bin(u_fft, v_fft, tk_fft, nx, ny, nmax, dkx, dky, dkh, dx, kx, ky):
    
    e_p = np.zeros(nmax)
    for p in range(nmax):
        kp = (p+1)*dkh
        c = 1
        for k in range(nx):
            for l in range(ny):
                mag = (kx[k]**2 + ky[l]**2)**(1/2)
                if mag < kp + dkh/2 and mag >= kp - dkh/2:
                    val = np.real(u_fft[l,k]*np.conj(u_fft[l,k]) + v_fft[l,k]*np.conj(v_fft[l,k]) + (CP/TR)*tk_fft[l,k]*np.conj(tk_fft[l,k]))
                    e_p[p] = e_p[p] + val
                    c = c + 1

        e_p[p] = e_p[p]*dx*dx*np.minimum(dkx,dky)*2*\
                    np.pi*kp/(np.minimum(dkx,dky)*c*8*nx*ny*np.pi**2)
            
    return e_p

def te_back_2d(ncfile_ctl, dx, zbot=0.2, ztop=18, nz=90, timeid=0):

    ### Read in the data
    u_ctl = np.asarray(getvar(ncfile_ctl,'ua',timeidx=timeid))
    v_ctl = np.asarray(getvar(ncfile_ctl,'va',timeidx=timeid))
    tk_ctl = np.asarray(getvar(ncfile_ctl,'tk',timeidx=timeid))
    ### Interpolate onto height surface
    z = np.linspace(zbot,ztop,nz)
    u_ctl = vinterp(ncfile_ctl,u_ctl,'ght_msl',z,timeidx=timeid)
    v_ctl = vinterp(ncfile_ctl,v_ctl,'ght_msl',z,timeidx=timeid)
    tk_ctl = vinterp(ncfile_ctl,tk_ctl,'ght_msl',z,timeidx=timeid)
    ### Extend periodically
    u_ctl = np.append(u_ctl,np.flip(u_ctl,axis=1),axis=1)
    v_ctl = np.append(v_ctl,-np.flip(v_ctl,axis=1),axis=1)
    tk_ctl = np.append(tk_ctl,np.flip(tk_ctl,axis=1),axis=1)
    ### Compute spectra
    ny = u_ctl.shape[1]
    nx = u_ctl.shape[2]
    lx = nx*dx
    ly = ny*dx
    dkx = 2*np.pi/lx
    dky = 2*np.pi/ly
    dkh = np.maximum(dkx,dky)
    kx = 2*np.pi*np.fft.fftfreq(nx,dx)
    ky = 2*np.pi*np.fft.fftfreq(ny,dx)
    nmax = int(np.ceil(np.sqrt(2)*np.maximum(nx/2,ny/2)))
    te = np.zeros((nz,nmax))
    for m in range(nz):
        u_fft = np.fft.fft2(u_ctl[m,:,:])
        v_fft = np.fft.fft2(v_ctl[m,:,:])
        tk_fft = np.fft.fft2(tk_ctl[m,:,:])
        te[m,:] = te_bin(u_fft, v_fft, tk_fft, nx, ny, nmax, dkx, dky, dkh, dx, kx, ky)
    
    return te

def dte_2d(ncfile_ctl, ncfile_pert, dx, zbot=0.2, ztop=18, nz=90, timeid=0):

    ### Read in the data
    u_ctl = np.asarray(getvar(ncfile_ctl,'ua',timeidx=timeid))
    v_ctl = np.asarray(getvar(ncfile_ctl,'va',timeidx=timeid))
    tk_ctl = np.asarray(getvar(ncfile_ctl,'tk',timeidx=timeid))
    u_pert = np.asarray(getvar(ncfile_pert,'ua',timeidx=timeid))
    v_pert = np.asarray(getvar(ncfile_pert,'va',timeidx=timeid))
    tk_pert = np.asarray(getvar(ncfile_pert,'tk',timeidx=timeid))
    ### Interpolate onto height surface
    z = np.linspace(zbot,ztop,nz)
    u_ctl = vinterp(ncfile_ctl,u_ctl,'ght_msl',z,timeidx=timeid)
    v_ctl = vinterp(ncfile_ctl,v_ctl,'ght_msl',z,timeidx=timeid)
    tk_ctl = vinterp(ncfile_ctl,tk_ctl,'ght_msl',z,timeidx=timeid)
    u_pert = vinterp(ncfile_pert,u_pert,'ght_msl',z,timeidx=timeid)
    v_pert = vinterp(ncfile_pert,v_pert,'ght_msl',z,timeidx=timeid)
    tk_pert = vinterp(ncfile_pert,tk_pert,'ght_msl',z,timeidx=timeid)
    ### Extend periodically
    u_ctl = np.append(u_ctl,np.flip(u_ctl,axis=1),axis=1)
    v_ctl = np.append(v_ctl,-np.flip(v_ctl,axis=1),axis=1)
    tk_ctl = np.append(tk_ctl,np.flip(tk_ctl,axis=1),axis=1)
    u_pert = np.append(u_pert,np.flip(u_pert,axis=1),axis=1)
    v_pert = np.append(v_pert,-np.flip(v_pert,axis=1),axis=1)
    tk_pert = np.append(tk_pert,np.flip(tk_pert,axis=1),axis=1)
    ### Compute spectra
    ny = u_ctl.shape[1]
    nx = u_ctl.shape[2]
    lx = nx*dx
    ly = ny*dx
    dkx = 2*np.pi/lx
    dky = 2*np.pi/ly
    dkh = np.maximum(dkx,dky)
    kx = 2*np.pi*np.fft.fftfreq(nx,dx)
    ky = 2*np.pi*np.fft.fftfreq(ny,dx)
    nmax = int(np.ceil(np.sqrt(2)*np.maximum(nx/2,ny/2)))
    du = u_ctl - u_pert
    dv = v_ctl - v_pert
    dtk = tk_ctl - tk_pert
    dte = np.zeros((nz,nmax))
    for m in range(nz):
        du_fft = np.fft.fft2(du[m,:,:])
        dv_fft = np.fft.fft2(dv[m,:,:])
        dtk_fft = np.fft.fft2(dtk[m,:,:])
        dte[m,:] = te_bin(du_fft, dv_fft, dtk_fft, nx, ny, nmax, dkx, dky, dkh, dx, kx, ky)
    
    return dte

#-------------------- KE --------------------#

def ke_back_1d(base_ncfile, dx, timeid=0):
    """
    Computes the one-dimensional zonal background kinetic energy spectrum at a given time
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param dx: float: grid spacing in x direction in meters
    :param time: int: time index (default 0)
    Returns a three-dimensional array
    """
    u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=timeid))
    v_base = np.asarray(getvar(base_ncfile,'va',timeidx=timeid))
    nz = u_base.shape[0]
    ny = u_base.shape[1]
    nx = u_base.shape[2]
    ke = np.zeros((nz,ny,nx))
    for k in range(nz):
        for j in range(ny):
            u_fft = np.fft.fft(u_base[k,j,:])
            v_fft = np.fft.fft(v_base[k,j,:])
            for i in range(nx):
                ke[k,j,i] = (dx/(2*np.pi*nx))*(np.abs(u_fft[i]*np.conj(u_fft[i])) + \
                                np.abs(v_fft[i]*np.conj(v_fft[i])))
                
    return ke

def ke_back_1dy(base_ncfile, dx, timeid=0):
    """
    Computes the one-dimensional zonal background kinetic energy spectrum at a given time
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param dx: float: grid spacing in x direction in meters
    :param time: int: time index (default 0)
    Returns a three-dimensional array
    """
    u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=timeid))
    v_base = np.asarray(getvar(base_ncfile,'va',timeidx=timeid))
    nz = u_base.shape[0]
    ny = u_base.shape[1]
    nx = u_base.shape[2]
    ke = np.zeros((nz,ny,nx))
    for k in range(nz):
        for j in range(nx):
            u_fft = np.fft.fft(u_base[k,:,j])
            v_fft = np.fft.fft(v_base[k,:,j])
            for i in range(ny):
                ke[k,i,j] = (dx/(2*np.pi*nx))*(np.abs(u_fft[i]*np.conj(u_fft[i])) + \
                                np.abs(v_fft[i]*np.conj(v_fft[i])))
                
    return ke

def ke_pert_1d(base_ncfile, pert_ncfile, dx, timeid=0):
    """
    Computes the one-dimensional zonal perturbation kinetic energy spectrum at a given time
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param pert_ncfile: netCDF Dataset: perturbed WRF restart file opened using netCDF4.Dataset()
    :param dx: float: grid spacing in x direction in meters
    :param time: int: time index (default 0)
    Returns a three-dimensional array
    """    
    u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=timeid))
    v_base = np.asarray(getvar(base_ncfile,'va',timeidx=timeid))
    u_pert = np.asarray(getvar(pert_ncfile,'ua',timeidx=timeid))
    v_pert = np.asarray(getvar(pert_ncfile,'va',timeidx=timeid))
    nz = u_base.shape[0]
    ny = u_base.shape[1]
    nx = u_base.shape[2]
    du = u_base - u_pert
    dv = v_base - v_pert
    ke = np.zeros((nz,ny,nx))
    for k in range(nz):
        for j in range(ny):
            du_fft = np.abs(fft(du[k,j,:]))
            dv_fft = np.abs(fft(dv[k,j,:]))
            for i in range(nx):
                ke[k,j,i] = (dx/(2*np.pi*nx))*(np.abs(du_fft[i]*du_fft[i]) + \
                                np.abs(dv_fft[i]*dv_fft[i]))
                
    return ke

@njit
def ke_bin(u_fft, v_fft, rhoz, nx, ny, nmax, dkx, dky, dkh, dx, kx, ky):
    
    e_p = np.zeros(nmax)
    for p in range(nmax):
        kp = (p+1)*dkh
        c = 1
        for k in range(nx):
            for l in range(ny):
                mag = (kx[k]**2 + ky[l]**2)**(1/2)
                if mag < kp + dkh/2 and mag >= kp - dkh/2:
                    val = np.real(u_fft[l,k]*np.conj(u_fft[l,k]) + v_fft[l,k]*np.conj(v_fft[l,k]))
                    e_p[p] = e_p[p] + val
                    c = c + 1

        e_p[p] = e_p[p]*rhoz*dx*dx*np.minimum(dkx,dky)*2*\
                    np.pi*kp/(np.minimum(dkx,dky)*c*8*nx*ny*np.pi**2)
            
    return e_p

def ke_back_2d(ncfile_ctl, dx, zbot=0.2, ztop=18, nz=90, timeid=0):

    ### Read in the data
    u_ctl = np.asarray(getvar(ncfile_ctl,'ua',timeidx=timeid))
    v_ctl = np.asarray(getvar(ncfile_ctl,'va',timeidx=timeid))
    p = np.asarray(getvar(ncfile_ctl,'p',timeidx=timeid))
    tv = np.asarray(getvar(ncfile_ctl,'tv',timeidx=timeid))
    rho = p/(RD*tv)
    ### Interpolate onto height surface
    z = np.linspace(zbot,ztop,nz)
    u_ctl = vinterp(ncfile_ctl,u_ctl,'ght_msl',z,timeidx=timeid)
    v_ctl = vinterp(ncfile_ctl,v_ctl,'ght_msl',z,timeidx=timeid)
    rho = vinterp(ncfile_ctl,rho,'ght_msl',z,timeidx=timeid)
    rho_mean = np.mean(rho,axis=tuple(range(1,3)))
    rho_mean = np.asarray(rho_mean)
    ### Extend periodically
    u_ctl = np.append(u_ctl,np.flip(u_ctl,axis=1),axis=1)
    v_ctl = np.append(v_ctl,-np.flip(v_ctl,axis=1),axis=1)
    ### Compute spectra
    ny = u_ctl.shape[1]
    nx = u_ctl.shape[2]
    lx = nx*dx
    ly = ny*dx
    dkx = 2*np.pi/lx
    dky = 2*np.pi/ly
    dkh = np.maximum(dkx,dky)
    kx = 2*np.pi*np.fft.fftfreq(nx,dx)
    ky = 2*np.pi*np.fft.fftfreq(ny,dx)
    nmax = int(np.ceil(np.sqrt(2)*np.maximum(nx/2,ny/2)))
    ke = np.zeros((nz,nmax))
    for m in range(nz):
        u_fft = np.fft.fft2(u_ctl[m,:,:])
        v_fft = np.fft.fft2(v_ctl[m,:,:])
        rhoz = rho_mean[m]
        ke[m,:] = ke_bin(u_fft, v_fft, rhoz, nx, ny, nmax, dkx, dky, dkh, dx, kx, ky)
    
    return ke

def ke_pert_2d(ncfile_ctl, ncfile_pert, dx, zbot=0.2, ztop=18, nz=90, timeid=0):

    ### Read in the data
    u_ctl = np.asarray(getvar(ncfile_ctl,'ua',timeidx=timeid))
    v_ctl = np.asarray(getvar(ncfile_ctl,'va',timeidx=timeid))
    u_pert = np.asarray(getvar(ncfile_pert,'ua',timeidx=timeid))
    v_pert = np.asarray(getvar(ncfile_pert,'va',timeidx=timeid))
    p = np.asarray(getvar(ncfile_ctl,'p',timeidx=timeid))
    tv = np.asarray(getvar(ncfile_ctl,'tv',timeidx=timeid))
    rho = p/(RD*tv)
    ### Interpolate onto height surface
    z = np.linspace(zbot,ztop,nz)
    u_ctl = vinterp(ncfile_ctl,u_ctl,'ght_msl',z,timeidx=timeid)
    v_ctl = vinterp(ncfile_ctl,v_ctl,'ght_msl',z,timeidx=timeid)
    u_pert = vinterp(ncfile_pert,u_pert,'ght_msl',z,timeidx=timeid)
    v_pert = vinterp(ncfile_pert,v_pert,'ght_msl',z,timeidx=timeid)
    rho = vinterp(ncfile_ctl,rho,'ght_msl',z,timeidx=timeid)
    rho_mean = np.mean(rho,axis=tuple(range(1,3)))
    rho_mean = np.asarray(rho_mean)
    ### Extend periodically
    u_ctl = np.append(u_ctl,np.flip(u_ctl,axis=1),axis=1)
    v_ctl = np.append(v_ctl,-np.flip(v_ctl,axis=1),axis=1)
    u_pert = np.append(u_pert,np.flip(u_pert,axis=1),axis=1)
    v_pert = np.append(v_pert,-np.flip(v_pert,axis=1),axis=1)
    ### Compute spectra
    ny = u_ctl.shape[1]
    nx = u_ctl.shape[2]
    lx = nx*dx
    ly = ny*dx
    dkx = 2*np.pi/lx
    dky = 2*np.pi/ly
    dkh = np.maximum(dkx,dky)
    kx = 2*np.pi*np.fft.fftfreq(nx,dx)
    ky = 2*np.pi*np.fft.fftfreq(ny,dx)
    nmax = int(np.ceil(np.sqrt(2)*np.maximum(nx/2,ny/2)))
    du = u_ctl - u_pert
    dv = v_ctl - v_pert
    ke = np.zeros((nz,nmax))
    for m in range(nz):
        du_fft = np.fft.fft2(du[m,:,:])
        dv_fft = np.fft.fft2(dv[m,:,:])
        rhoz = rho_mean[m]
        ke[m,:] = ke_bin(du_fft, dv_fft, rhoz, nx, ny, nmax, dkx, dky, dkh, dx, kx, ky)
    
    return ke

#-------------------- FSS --------------------#
@njit    
def _fss(io, im, kernel, norm=1.):
    """
    Nested function for computing FSS
    """
    ny = io.shape[0]
    nx = io.shape[1]
    nl = kernel.shape[0]
    nk = kernel.shape[1]
    ks = np.sum(kernel)
    nl2 = (nl-1) // 2
    nk2 = (nk-1) // 2
    o_array = np.zeros((ny-nl+1, nx-nk+1))
    m_array = np.zeros_like(o_array)
    for j in range(nl2, ny-nl2):
        for i in range(nk2, nx-nk2):
            test_io = np.sum(io[j-nl2:j+nl-nl2, i-nk2:i+nk-nk2])
            test_im = np.sum(im[j-nl2:j+nl-nl2, i-nk2:i+nk-nk2])
            if test_io > 1.e-10:
                for l in range(nl):
                    for k in range(nk):
                        o_array[j-nl2, i-nk2] += io[j+l-nl2, i+k-nk2] * kernel[l, k] / ks
            if test_im > 1.e-10:
                for l in range(nl):
                    for k in range(nk):
                        m_array[j-nl2, i-nk2] += im[j+l-nl2, i+k-nk2] * kernel[l, k] / ks
    mse_ = np.sum((o_array - m_array)**2)
    mse_ /= 1. * norm  # 1. * (ny-nl+1) * (nx-nk+1)
    ref_ = np.sum(o_array**2 + m_array**2)
    ref_ /= 1. * norm
    if ref_ <= 1.e-10:
        fss_ = 1.
    else:
        fss_ = 1. - mse_ / ref_
    return fss_, mse_

def fss(modeled, observed, threshold, neighborhood=1, kernel='square', inverse_threshold=False, return_mse=False,
        verbose=False):
    """
    Calculates the Fractions Skill Score of a modeled field given the observed field. The threshold parameter sets the
    threshold value for the FSS calculation, while the neighborhood is the number of points away from the center point
    to consider in the calculation (if it is zero, only the center point is used). The kernel can either be 'square',
    in which case all values within a square around each grid point are considered, or 'circle', where only points
    within a neighborhood radius away from the center are considered. If inverse_threshold is True, then we look for
    values LOWER than the threshold value.
    :param modeled: ndarray: modeled values. Acts on the last two dimensions.
    :param observed: ndarray: observed values. Must match dimensions of modeled.
    :param threshold: float: threshold value
    :param neighborhood: int: grid-point neighborhood radius (default 1)
    :param kernel: str: 'square' or 'circle' (default 'square')
    :param inverse_threshold: set to True if values BELOW threshold are desired (default False)
    :param return_mse: set to True to also return MSE values.
    :param verbose: set to True to get progress output (default False)
    :return: float or ndarray: FSS scores, and MSE scores if desired. Returns a single value if modeled is
    2-dimensional, otherwise returns an ndarray of the size modeled.shape[0].
    """
    # Check dimensions
    dims = modeled.shape
    dims_observed = observed.shape
    if dims != dims_observed:
        raise ValueError("Dimensions of 'modeled' must match those of 'observed'; got %s and %s" %
                         (dims, dims_observed))
    if len(dims) > 2:
        first_dims = dims[:-2]
        nz = int(np.prod(first_dims))
        multi_dims = True
    else:
        multi_dims = False
    ny, nx = dims[-2:]

    # Create the kernel array
    if verbose:
        print('fss: initializing kernel and binary arrays')
    kernel_dim = 2 * neighborhood + 1
    if kernel_dim > dims[-1] or kernel_dim > dims[-2]:
        raise ValueError('neighborhood size (%d) must be smaller than 1/2 the smallest modeled array dimension (%d)' %
                         (neighborhood, min(ny, nx)))
    if kernel == 'square':
        kernel_array = np.ones((kernel_dim, kernel_dim))
    elif kernel == 'circle':
        kernel_array = np.zeros((kernel_dim, kernel_dim))
        x, y = np.meshgrid(np.arange(kernel_dim), np.arange(kernel_dim))
        kernel_array[np.sqrt((x - neighborhood) ** 2 + (y - neighborhood) ** 2) <= neighborhood] = 1.
    else:
        raise ValueError("kernel must be 'square' or 'circle'")

    # Create the I_O  and I_M arrays
    I_M = np.zeros_like(modeled)
    I_O = np.zeros_like(observed)
    if inverse_threshold:
        I_M[modeled <= threshold] = 1.
        I_O[observed <= threshold] = 1.
    else:
        I_M[modeled >= threshold] = 1.
        I_O[observed >= threshold] = 1.

    # Calculate FSS
    if multi_dims:
        I_M = np.reshape(I_M, (nz, ny, nx))
        I_O = np.reshape(I_O, (nz, ny, nx))
        fss_ = np.zeros(nz, dtype=modeled.dtype)
        mse_ = np.zeros_like(fss_)
        normalization = np.maximum(np.sum(I_O, axis=(-2, -1)), np.ones(nz))
        for z in range(nz):
            if verbose:
                print('fss: calculating FSS for index %d of %d' % (z+1, nz))
            fss_[z], mse_[z] = _fss(I_O[z, :, :], I_M[z, :, :], kernel_array, normalization[z])
        fss_ = np.reshape(fss_, first_dims)
        mse_ = np.reshape(mse_, first_dims)
    else:
        normalization = np.max([1., np.sum(I_O)])
        if verbose:
            print('fss: calculating FSS')
        fss_, mse_ = _fss(I_O, I_M, kernel_array, normalization)

    if return_mse:
        return fss_, mse_
    else:
        return fss_
    
#-------------------- RESTART FILES --------------------#
def dte_rst(base_ncfile, pert_ncfile):
    """
    Computes the domain-integrated difference total energy from a WRF restart file
    :param base_ncfile: netCDF Dataset: control WRF restart file opened using netCDF4.Dataset()
    :param pert_ncfile: netCDF Dataset: perturbed WRF restart file opened using netCDF4.Dataset()
    Returns a float
    """
    u_base = np.asarray(getvar(base_ncfile,'U_2',timeidx=0))
    v_base = np.asarray(getvar(base_ncfile,'V_2',timeidx=0))
    theta_base = np.asarray(getvar(base_ncfile,'T_2',timeidx=0)) + TR
    p_base = np.asarray(getvar(base_ncfile,'P',timeidx=0))/100.
    pb_base = np.asarray(getvar(base_ncfile,'PB',timeidx=0))/100.
    tk_base = theta_base*((p_base + pb_base)/P0)**(RD/CP)
    u_pert = np.asarray(getvar(pert_ncfile,'U_2',timeidx=0))
    v_pert = np.asarray(getvar(pert_ncfile,'V_2',timeidx=0))
    theta_pert = np.asarray(getvar(pert_ncfile,'T_2',timeidx=0)) + TR
    p_pert = np.asarray(getvar(pert_ncfile,'P',timeidx=0))/100.
    pb_pert = np.asarray(getvar(pert_ncfile,'PB',timeidx=0))/100.
    tk_pert = theta_pert*((p_pert + pb_pert)/P0)**(RD/CP)
    du = u_pert - u_base
    dv = v_pert - v_base
    dtk = tk_pert - tk_base
    dte1 = np.sum(np.sum(np.sum(np.squeeze(du[:,:,:]**2.0),1),1))
    dte2 = np.sum(np.sum(np.sum(np.squeeze(dv[:,:,:]**2.0),1),1))
    dte3 = np.sum(np.sum(np.sum(np.squeeze(CP/TR*(dtk[:,:,:]**2.0)),1),1))
    dte = 0.5*(dte1 + dte2 + dte3)
    
    return dte

def dte_1d_rst(base_ncfile, pert_ncfile, dx):
    """
    Computes the one-dimensional zonal difference total energy spectrum from a WRF restart file
    :param base_ncfile: netCDF Dataset: control WRF restart file opened using netCDF4.Dataset()
    :param pert_ncfile: netCDF Dataset: perturbed WRF restart file opened using netCDF4.Dataset()
    :param dx: float: grid spacing in x direction in meters
    Returns a three-dimensional array
    """ 
    u_base = np.asarray(getvar(base_ncfile,'U_2',timeidx=0))
    v_base = np.asarray(getvar(base_ncfile,'V_2',timeidx=0))
    theta_base = np.asarray(getvar(base_ncfile,'T_2',timeidx=0)) + TR
    p_base = np.asarray(getvar(base_ncfile,'P',timeidx=0))/100.
    pb_base = np.asarray(getvar(base_ncfile,'PB',timeidx=0))/100.
    tk_base = theta_base*((p_base + pb_base)/P0)**(R/CP)
    u_pert = np.asarray(getvar(pert_ncfile,'U_2',timeidx=0))
    v_pert = np.asarray(getvar(pert_ncfile,'V_2',timeidx=0))
    theta_pert = np.asarray(getvar(pert_ncfile,'T_2',timeidx=0)) + TR
    p_pert = np.asarray(getvar(pert_ncfile,'P',timeidx=0))/100.
    pb_pert = np.asarray(getvar(pert_ncfile,'PB',timeidx=0))/100.
    tk_pert = theta_pert*((p_pert + pb_pert)/P0)**(R/CP)   
    nz = u_base.shape[0]
    ny = u_base.shape[1]
    nx = u_base.shape[2]
    du = u_base - u_pert
    dv = v_base - v_pert
    dt = tk_base - tk_pert
    dte_spec = np.zeros((nz,ny,nx))
    for k in range(nz):
        for j in range(ny):
            du_fft = np.abs(fft(du[k,j,:]))
            dv_fft = np.abs(fft(dv[k,j,:]))
            dt_fft = np.abs(fft(dt[k,j,:]))
            for i in range(nx):
                dte_spec[k,j,i] = (dx/(2*np.pi*nx))*(np.abs(du_fft[i]*du_fft[i]) + \
                                np.abs(dv_fft[i]*dv_fft[i]) + \
                                np.abs(CP/TR*dt_fft[i]*dt_fft[i]))
                
    return dte_spec

def ke_pert_1d_rst(base_ncfile, pert_ncfile, dx):
    """
    Computes the one-dimensional zonal perturbation kinetic energy spectrum at the initial time
    :param base_ncfile: netCDF Dataset: control WRF output file opened using netCDF4.Dataset()
    :param pert_ncfile: netCDF Dataset: perturbed WRF restart file opened using netCDF4.Dataset()
    :param dx: float: grid spacing in x direction in meters
    Returns a three-dimensional array
    """ 
    u_base = np.asarray(getvar(base_ncfile,'U_1',timeidx=0))
    v_base = np.asarray(getvar(base_ncfile,'V_1',timeidx=0))
    u_pert = np.asarray(getvar(pert_ncfile,'U_1',timeidx=0))
    v_pert = np.asarray(getvar(pert_ncfile,'V_1',timeidx=0))
    nz = u_base.shape[0]
    ny = u_base.shape[1]
    nx = u_base.shape[2]
    du = u_base - u_pert
    dv = v_base - v_pert
    ke = np.zeros((nz,ny,nx))
    for k in range(nz):
        for j in range(ny):
            du_fft = np.abs(fft(du[k,j,:]))
            dv_fft = np.abs(fft(dv[k,j,:]))
            for i in range(nx):
                ke[k,j,i] = (dx/(2*np.pi*nx))*(np.abs(du_fft[i]*du_fft[i]) + \
                                np.abs(dv_fft[i]*dv_fft[i]))
                
    return ke
