#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:31:29 2020

@author: Daniel Lloveras

Defines functions for computing moist soundings for fields on the WRF grid
"""
import numpy as np
from numba import njit

#%% Constants

G = 9.81  # gravitational acceleration, in m/s^2
T0 = 300.  # reference potential temperature, in K
P0 = 1.0e5  # reference pressure, in Pa
CP = 1004.  # specific heat at constant pressure, in J/(K*kg)       
CV = 717.  # specific heat at constant volume, in J/(K*kg)
RD = 287.  # ideal gas constant for dry air, in J/(K*kg)
RV = 461.6 # ideal gas constant for water vapor, in J/(K*kg)
F0 = 1.0e-4  # Coriolis parameter, in s^-1
SVPT0 = 273.15
GAMMA = CP/CV
KAPPA = RD/CP

#%% WRF grid

@njit
def wrf_grid(nx, ny, nz, hres, znw):
    
    dnw = np.zeros(nz)
    rdnw = np.zeros(nz)
    for k in range(nz):
        dnw[k] = znw[k+1] - znw[k]
        rdnw[k] = 1./dnw[k]
    
    dn = np.zeros(nz)
    rdn = np.zeros(nz)
    fnp = np.zeros(nz)
    fnm = np.zeros(nz)
    for k in range(1,nz):
        dn[k] = 0.5*(dnw[k] + dnw[k-1])
        rdn[k] = 1./dn[k]
        fnp[k] = 0.5*dnw[k]/dn[k]
        fnm[k] = 0.5*dnw[k-1]/dn[k]
        
    cof1 = (2.*dn[1] + dn[2])/(dn[1] + dn[2])*dnw[0]/dn[1]
    cof2 = dn[1]/(dn[1] + dn[2])*dnw[0]/dn[2]
    cf1 = fnp[1] + cof1
    cf2 = fnm[1] - cof1 - cof2
    cf3 = cof2
    cfn = (0.5*dnw[nz-1] + dn[nz-1])/dn[nz-1]
    cfn1 = -0.5*dnw[nz-1]/dn[nz-1]
    dx = hres*1000.
    dy = hres*1000.
    rdx = 1./dx
    rdy = 1./dy
    
    return dnw, rdnw, dn, rdn, fnp, fnm, cof1, cof2, cf1, cf2, cf3, cfn, \
            cfn1, rdx, rdy
            
#%% Interpolation function

@njit
def interp_0(v_in, z_in, z_val, nz_in):
    
    v_val = 0.
    height = z_val
    if (z_in[nz_in-1] > z_in[0]):
        if (height > z_in[nz_in-1]):
            w2 = (z_in[nz_in-1] - height)/(z_in[nz_in-1] - z_in[nz_in-2])
            w1 = 1. - w2
            v_val = w1*v_in[nz_in-1] + w2*v_in[nz_in-2]
        
        elif (height < z_in[0]):
            w2 = (z_in[1] - height)/(z_in[1] - z_in[0])
            w1 = 1. - w2
            v_val = w1*v_in[1] + w2*v_in[0]
            
        else:
            interp = False
            kp = nz_in-1
            while (interp == False and kp >= 1):
                if (z_in[kp] >= height and z_in[kp-1] <= height):
                    w2 = (height - z_in[kp])/(z_in[kp-1] - z_in[kp])
                    w1 = 1. - w2
                    v_val = w1*v_in[kp] + w2*v_in[kp-1]
                    interp = True

                kp = kp - 1                    
    
    else:
        if (height < z_in[nz_in-1]):
            w2 = (z_in[nz_in-1] - height)/(z_in[nz_in-1] - z_in[nz_in-2])
            w1 = 1. - w2
            v_val = w1*v_in[nz_in-1] + w2*v_in[nz_in-2]
        
        elif (height > z_in[0]):
            w2 = (z_in[1] - height)/(z_in[1] - z_in[0])
            w1 = 1. - w2
            v_val = w1*v_in[1] + w2*v_in[0]
            
        else:
            interp = False
            kp = nz_in-1
            while (interp == False and kp >= 1):
                if (z_in[kp] <= height and z_in[kp-1] >= height):
                    w2 = (height - z_in[kp])/(z_in[kp-1] - z_in[kp])
                    w1 = 1. - w2
                    v_val = w1*v_in[kp] + w2*v_in[kp-1]
                    interp = True

                kp = kp - 1          
            
    return v_val

#%% Add QGPV perturbation

@njit    
def add_pert(u_jet, theta_jet, rho_jet, z_jet, p_jet, z_pert, u_pert, v_pert, \
             rho_pert, theta_pert, p_pert, nx, ny, nz):
    
    ### Extend jet variables to 3D
    u = np.zeros((nz,ny,nx))
    v = np.zeros((nz,ny,nx))
    theta = np.zeros((nz,ny,nx))
    rho = np.zeros((nz,ny,nx))
    z = np.zeros((nz,ny,nx))
    pres = np.zeros((nz,ny,nx))

    for i in range(nx):
        u[:,:,i] = u_jet
        v[:,:,i] = np.zeros((nz,ny))
        theta[:,:,i] = theta_jet
        rho[:,:,i] = rho_jet
        z[:,:,i] = z_jet
        pres[:,:,i] = p_jet
        
    up_eta = np.zeros((nz,ny,nx))
    vp_eta = np.zeros((nz,ny,nx))
    thetap_eta = np.zeros((nz,ny,nx))
    rhop_eta = np.zeros((nz,ny,nx))
    presp_eta = np.zeros((nz,ny,nx))
    ### Add perturbation
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                thetap_eta[k,j,i] = interp_0(theta_pert[:,j,i], z_pert, z[k,j,i], nz)
                rhop_eta[k,j,i] = interp_0(rho_pert[:,j,i], z_pert, z[k,j,i], nz)            
                up_eta[k,j,i] = interp_0(u_pert[:,j,i], z_pert, z[k,j,i], nz)
                vp_eta[k,j,i] = interp_0(v_pert[:,j,i], z_pert, z[k,j,i], nz)
                presp_eta[k,j,i] = interp_0(p_pert[:,j,i], z_pert, z[k,j,i], nz)            
    
    u_eta = u + up_eta
    v_eta = v + vp_eta
    rho_eta = rho + rhop_eta
    theta_eta = theta + thetap_eta
    pres_eta = pres + presp_eta
    z_eta = z

    return u_eta, v_eta, theta_eta, rho_eta, pres_eta, z_eta, up_eta, vp_eta, thetap_eta, rhop_eta, presp_eta

#%% Moist sounding
    
@njit
def calc_sounding(i, j, nz, z_eta, theta_eta, rho_eta, u_eta, v_eta, pres_eta, p_bot, avg_rh):
    
    hprof = np.zeros(nz)
    thprof = np.zeros(nz)
    qvprof = np.zeros(nz)
    rhoprof = np.zeros(nz)
    uprof = np.zeros(nz)
    vprof = np.zeros(nz)
    pprof = np.zeros(nz)
    pmprof = np.zeros(nz)
    pdprof = np.zeros(nz)
    rel_hum = np.ones(nz)*avg_rh
    for k in range(nz):
        hprof[k] = z_eta[k,j,i]
        thprof[k] = theta_eta[k,j,i]
        rhoprof[k] = rho_eta[k,j,i]
        uprof[k] = u_eta[k,j,i]
        vprof[k] = v_eta[k,j,i]
        pprof[k] = pres_eta[k,j,i]
        tmppi = (pprof[k]/P0)**KAPPA
        temp = tmppi*thprof[k]
        if (temp > SVPT0):
            es = 1000.*0.6112*np.exp(17.67*(temp - SVPT0)/(temp - 29.65))
            qvs = (RD/RV)*es/(pprof[k] - es)
            
        else:
            es = 1000.*0.6112*np.exp(21.87*(temp - SVPT0)/(temp - 7.66))
            qvs = (RD/RV)*es/(pprof[k] - es)
            
        qvprof[k] = rel_hum[k]*qvs
        thprof[k] = thprof[k]/(1. + 0.61*qvprof[k])
        
    pprof_surf = p_bot[j,i]
    thprof_surf = 1.5*thprof[0] - 0.5*thprof[1]
    
    qvf = 1. + (RV/RD)*qvprof[0]
    qvf1 = 1. + qvprof[0]
    rhoprof_surf = 1./((RD/P0)*thprof_surf*qvf*((pprof_surf/P0)**(-1/GAMMA)))
    dz1 = hprof[0]
    for it in range(10):
        pmprof[0] = pprof_surf - 0.5*dz1*(rhoprof_surf + rhoprof[0])*G*qvf1
        rhoprof[0] = 1./((RD/P0)*thprof[0]*qvf*((pmprof[0]/P0)**(-1/GAMMA)))
        
    for k in range(1,nz):
        dz = hprof[k] - hprof[k-1]
        qvf1 = 1. + 0.5*(qvprof[k-1] + qvprof[k])
        qvf = 1. + (RV/RD)*qvprof[k]
        for it in range(10):
            pmprof[k] = pmprof[k-1] - 0.5*dz*(rhoprof[k] + rhoprof[k-1])*G*qvf1
            rhoprof[k] = 1./((RD/P0)*thprof[k]*qvf*((pmprof[k]/P0)**(-1/GAMMA)))
    
    pdprof[nz-1] = pmprof[nz-1]
    for k in range(nz-2,-1,-1):
        dz = hprof[k+1] - hprof[k]
        pdprof[k] = pdprof[k+1] + 0.5*dz*(rhoprof[k] + rhoprof[k+1])*G
        
    return hprof, pmprof, pdprof, thprof, rhoprof, uprof, vprof, qvprof
            
#%% Middle of domain
    
@njit
def middle_sounding(nx, ny, nz, p_bot, p_top, znu, dnw, z_eta, theta_eta, rho_eta, u_eta, v_eta, pres_eta, avg_rh):
        
    ip = int(nx/2) - 1
    jp = int(ny/2) - 1
    zmid, pmmid, pdmid, thmid, rhomid, umid, vmid, qvmid = calc_sounding(ip, \
        jp, nz, z_eta, theta_eta, rho_eta, u_eta, v_eta, pres_eta, p_bot, avg_rh)

    p_ex = np.append(p_bot[jp,ip],pmmid[:])
    z_ex = np.append(0.,zmid[:])
    mub = np.zeros((ny,nx))
    pb = np.zeros((nz,ny,nx))
    t_init = np.zeros((nz,ny,nx))
    alb = np.zeros((nz,ny,nx))
    phb = np.zeros((nz+1,ny,nx))
    for j in range(ny):
        for i in range(nx):
            ps = interp_0(p_ex, z_ex, 0., nz+1)
            mub[j,i] = ps - p_top[jp,ip]
            for k in range(nz):
                pb[k,j,i] = znu[k]*(ps - p_top[jp,ip]) + p_top[jp,ip]
                t_init[k,j,i] = interp_0(thmid, pmmid, pb[k,j,i], nz) - T0
                alb[k,j,i] = (RD/P0)*(t_init[k,j,i] + T0)*((pb[k,j,i]/P0)**(-1/GAMMA))
                
            for k in range(1,nz+1):    
                phb[k,j,i] = phb[k-1,j,i] - dnw[k-1]*mub[j,i]*alb[k-1,j,i]
    
    return mub, pb, t_init, alb, phb

#%% Full domain
    
@njit
def full_domain(nx, ny, nz, p_bot, p_top, znu, dnw, mub, pb, alb, dn, \
                z_eta, theta_eta, rho_eta, u_eta, v_eta, pres_eta, avg_rh):

    mu = np.zeros((ny,nx))
    moist = np.zeros((nz,ny,nx))
    t = np.zeros((nz,ny,nx))
    p = np.zeros((nz,ny,nx))
    alt = np.zeros((nz,ny,nx))
    al = np.zeros((nz,ny,nx))
    ph = np.zeros((nz+1,ny,nx))
    u = np.zeros((nz,ny,nx+1))
    v = np.zeros((nz,ny+1,nx))
    for j in range(ny):
        for i in range(nx):
            z1, pm1, pd1, th1, rho1, u1, v1, qv1 = calc_sounding(i, j, nz, z_eta,\
                theta_eta, rho_eta, u_eta, v_eta, pres_eta, p_bot, avg_rh)
             
            ### Dry sounding
            p_ex = np.append(p_bot[j,i],pd1[:])
            z_ex = np.append(0.,z1[:])
            pd_surf = interp_0(p_ex, z_ex, 0., nz+1)
        
            ### Add moisture to sounding and interpolate
            for k in range(nz):
                p_level = znu[k]*(pd_surf - p_top[j,i]) + p_top[j,i]
                moist[k,j,i] = interp_0(qv1, pd1, p_level, nz)
                t[k,j,i] = interp_0(th1, pd1, p_level, nz) - T0
            
            ### Compute moist fields at top of atmosphere    
            qv_avg = moist[nz-1,j,i]
            mu[j,i] = (1 + qv_avg)*(pd_surf - p_top[j,i]) - mub[j,i]
            p[nz-1,j,i] = -0.5*((1. + qv_avg)*mu[j,i] + qv_avg*mub[j,i])*dnw[nz-1]
            qvf = 1. + (RV/RD)*moist[nz-1,j,i]
            alt[nz-1,j,i] = (RD/P0)*(t[nz-1,j,i] + T0)*qvf*(((p[nz-1,j,i] + pb[nz-1,j,i])/P0)**(-1/GAMMA))
            al[nz-1,j,i] = alt[nz-1,j,i] - alb[nz-1,j,i]
        
            ### Compute moist fields down the column
            for k in range(nz-2,-1,-1):
                qv_avg = 0.5*(moist[k,j,i] + moist[k+1,j,i])
                p[k,j,i] = p[k+1,j,i] - ((1. + qv_avg)*mu[j,i] + qv_avg*mub[j,i])*dn[k+1]
                qvf = 1. + (RV/RD)*moist[k,j,i]
                alt[k,j,i] = (RD/P0)*(t[k,j,i] + T0)*qvf*(((p[k,j,i] + pb[k,j,i])/P0)**(-1/GAMMA))
                al[k,j,i] = alt[k,j,i] - alb[k,j,i]
            
            ### Compute geopotential
            for k in range(1,nz+1):
                ph[k,j,i] = ph[k-1,j,i] - dnw[k-1]*((mub[j,i] + mu[j,i])*al[k-1,j,i] + \
                        mu[j,i]*alb[k-1,j,i])
        
            ### Interpolate u and v    
            p_surf_val = interp_0(p_ex, z_ex, 0, nz+1)
            for k in range(nz):
                p_level = znu[k]*(p_surf_val - p_top[j,i]) + p_top[j,i]
                u[k,j,i] = interp_0(u1, pm1, p_level, nz)
                v[k,j,i] = interp_0(v1, pm1, p_level, nz)
      
    ### Zonal periodicity      
    for j in range(ny+1):
        for k in range(0,nz):
            v[k,j,nx-1] = v[k,j,0]
               
    for j in range(ny):
        for k in range(nz):
            u[k,j,nx] = u[k,j,0]         
            
    tsk = np.zeros((ny,nx))
    for j in range(0,ny):
        for i in range(0,nx):
            thtmp = t[0,j,i] + T0
            ptmp = p[0,j,i] + pb[0,j,i]
            tsk[j,i] = thtmp*((ptmp/P0)**(RD/CP))        
    
    return u, v, t, ph, mu, p, moist, tsk
