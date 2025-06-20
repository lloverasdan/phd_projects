#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 09:38:33 2020

@author: Daniel Lloveras

Defines functions for inverting a 2D EPV distribution to create a zonal jet and
a meridional temperature gradient
"""

import numpy as np
from numba import njit

#%% Constants

F0 = 1.0e-4  # Coriolis parameter, in s^-1
G = 9.81  # gravitational acceleration, in m/s^2
T0 = 300.  # reference potential temperature, in K
P0 = 1.0e5  # reference pressure, in Pa
CP = 1004.  # specific heat at constant pressure, in J/(K*kg)       
CV = 717.  # specific heat at constant volume, in J/(K*kg)
RD = 287.  # ideal gas constant for dry air, in J/(K*kg)
GAMMA = CP/CV
KAPPA = RD/CP

#%% Grids

@njit
def y_grid(ly, ny_l, ny_h):
    
    dy_l = ly/(ny_l - 1)  # grid res in y for low-res run
    ylvls_l = np.zeros(ny_l)  # low-res y grid values
    for j in range(ny_l):
        ylvls_l[j] = j*dy_l
        
    dy_h = ly/(ny_h - 1)  # grid res in y for high-res run
    ylvls_h = np.zeros(ny_h)  # high-res y grid values
    for j in range(ny_h):
        ylvls_h[j] = j*dy_h
        
    return dy_l, ylvls_l, dy_h, ylvls_h

@njit
def pi_grid(pbot, pibot, pitop, npi_l, npi_h):
    
    dpi_l = (pibot - pitop)/(npi_l - 1)  # grid res in pi for low-res run
    pilvls_l = np.zeros(npi_l)  # low-res pi grid values
    p_pi_l = np.zeros(npi_l)
    for k in range(npi_l):
        pilvls_l[k] = pibot - k*dpi_l
        #p_pi_l[k] = P0*(pilvls_l[k]/CP)**(1/KAPPA)
        p_pi_l[k] = pbot*(pilvls_l[k]/CP)**(1/KAPPA)
        
    dpi_h = (pibot - pitop)/(npi_h - 1)  # grid res in pi for high-res run  
    pilvls_h = np.zeros(npi_h)  # high-res pi grid values
    p_pi_h = np.zeros(npi_h)
    for k in range(npi_h):
        pilvls_h[k] = pibot - k*dpi_h
        #p_pi_h[k] = P0*(pilvls_h[k]/CP)**(1/KAPPA)
        p_pi_h[k] = pbot*(pilvls_h[k]/CP)**(1/KAPPA)
        
    return dpi_l, pilvls_l, dpi_h, pilvls_h, p_pi_l, p_pi_h

@njit
def eta_grid(nz):
    
    znw = np.zeros(nz+1)  # staggered eta levels
    for k in range(nz+1):
        znw[k] = (np.exp(-2*k/float(nz)) - np.exp(-2))/(1 - np.exp(-2))

    #znw = np.linspace(1,0,nz+1)
    
    znu = np.zeros(nz)  # unstaggered eta levels
    for k in range(nz):
        znu[k] = 0.5*(znw[k+1] + znw[k])
        
    return znw, znu

#%% Interpolation functions

### Increasing grid values
@njit    
def interp_increasing(x_in, x_out, y_in, n_in, n_out):
    
    y_out = np.zeros(n_out)
    i_in = 0
    for i in range(n_out):
        while True:
            if (x_out[i] >= x_in[i_in] and x_out[i] <= x_in[i_in+1]):
                dx = x_out[i] - x_in[i_in]
                dx_in = x_in[i_in+1] - x_in[i_in]
                y_out[i] = y_in[i_in] + (y_in[i_in+1] - y_in[i_in])*dx/dx_in
                
                break
            
            else:
                i_in = i_in+1
                if (i_in > n_in-1):
                    break
                
        else:
            continue
        
        if (i_in > n_in-1):
            break

    return y_out

### Decreasing grid values
@njit    
def interp_decreasing(x_in, x_out, y_in, n_in, n_out):
    
    y_out = np.zeros(n_out)
    i_in = 0
    for i in range(n_out):
        while True:
            if (x_out[i] <= x_in[i_in] and x_out[i] >= x_in[i_in+1]):
                dx = x_out[i] - x_in[i_in]
                dx_in = x_in[i_in+1] - x_in[i_in]
                y_out[i] = y_in[i_in] + (y_in[i_in+1] - y_in[i_in])*dx/dx_in
                
                break
            
            else:
                i_in = i_in+1
                if (i_in > n_in-1):
                    break
                
        else:
            continue
        
        if (i_in > n_in-1):
            break

    return y_out

#%% Theta at top of atmosphere
    
@njit
def theta_top(ly, ny, dy, dyth, thtop, dth):

    thtopline = np.zeros(ny)
    for j in range(ny):
        a2 = 1.5*(j*dy - ly/2.)/dyth
        if (np.abs(a2) < np.pi/2.):
            thtopline[j] = thtop + dth*np.sin(a2)
        elif (a2 > np.pi/2.):
            thtopline[j] = thtop + dth
        else:
            thtopline[j] = thtop - dth

    return thtopline

#%% Tropopause shape
    
@njit
def trop_shape(ly, ny, dy, dytr, pim, dpitr):
    
    pitp = np.zeros(ny)
    for j in range(ny):
        a1 = 2.*(j*dy - ly/2.)/dytr
        if (np.abs(a1) <= np.pi/2.):
            pitp[j] = pim + dpitr*np.sin(a1)
        elif (a1 > np.pi/2.):
            pitp[j] = pim + dpitr
        else:
            pitp[j] = pim - dpitr
    
    return pitp
    
#%% PV distribution

@njit
def pv_dist(pvt, pvs, dpipv, ny, npi, pibot, pitop, dpi, pitp):
    
    pv = np.zeros((npi,ny))
    for j in range(ny):
        for k in range(npi):
            pilev = pibot - k*dpi
            gamma_val = (pibot - pilev)/(pibot - pitp[j])
            beta_val = (pitp[j] - pilev)/(pitp[j] - pitop)
            if (pitp[j] > pilev):
                gamma_val = 0.
            else:
                beta_val = 0.
                
            a = 1. + 3.0*(gamma_val**2.)
            b = 1. + 5.0*(beta_val**2.)                   
            pv[k,j] = (pvt*a + pvs*b)/2. + \
                        (pvt*a - pvs*b)/2.*(np.tanh(2.*(pilev - pitp[j])/dpipv))
      
    return pv

#%% PV inversion
    
@njit
def pv_inv(pv, f, ny, pibot, dy, npi, dpi, om):    

    fct = ((dy*dpi)**2.)*F0*P0/(G*KAPPA*(CP**(1/KAPPA)))
    for k in range(1,npi-1):
        pilev = pibot - k*dpi
        for j in range(1,ny-1):
            b = -0.5*(f[k,j+1] + f[k,j-1] + f[k+1,j] + f[k-1,j] + (F0*dy)**2.)
            c = 0.25*(-fct*pilev**(-(1 - (1/KAPPA)))*pv[k,j] \
                - (1/16.)*((f[k+1,j+1] - f[k+1,j-1] - f[k-1,j+1] + f[k-1,j-1])**2.) \
                + (f[k,j+1] + f[k,j-1])*(f[k+1,j] + f[k-1,j]) \
                + ((F0*dy)**2.)*(f[k+1,j] + f[k-1,j]))
            d = np.abs(b**2. - 4.*c)
            r = 0.5*(-b - np.sqrt(d)) - f[k,j]
            f[k,j] = f[k,j] + om*r
    
    return f

#%% U, Theta, and PV computations
    
@njit
def u_calc(f, ny, npi, dy):    

    u = np.zeros((npi,ny))
    for j in range(1,ny-1):
        u[:,j] = -(1/F0)*(f[:,j+1] - f[:,j-1])/(2.*dy)   
    
    u[:,0] = -(1/F0)*(f[:,1] - f[:,0])/dy
    u[:,ny-1] = -(1/F0)*(f[:,ny-1] - f[:,ny-2])/dy     

    return u

@njit
def theta_calc(f, ny, npi, dpi):

    theta = np.zeros((npi,ny))
    for k in range(1,npi-1):
        theta[k,:] = (f[k+1,:] - f[k-1,:])/(2.*dpi)

    theta[0,:] = (f[1,:] - f[0,:])/dpi 
    theta[npi-1,:] = (f[npi-1,:] - f[npi-2,:])/dpi        

    return theta
    
@njit
def pv_calc(f, ny, dy, npi, dpi, pibot):

    pv_out = np.zeros((npi,ny))
    for k in range(1,npi-1):
        pilev = pibot - k*dpi
        for j in range(1,ny-1):
            pv_out[k,j] = (G*KAPPA*(CP**(1/KAPPA))/P0)*pilev**(1 - (1/KAPPA)) \
                        *((F0*(f[k+1,j] - 2*f[k,j] + f[k-1,j])/(dpi**2.)) \
                        - ((1/F0)*(1/(4*dy*dpi)**2.)*(f[k+1,j+1] - f[k+1,j-1] \
                        - f[k-1,j+1] + f[k-1,j-1])**2) \
                        + ((f[k,j+1] - 2*f[k,j] + f[k,j-1])/(F0*dy*dy)) \
                        *((f[k+1,j] - 2*f[k,j] + f[k-1,j])/(dpi**2.)))
                        
    pv_out[0,:] = pv_out[1,:]
    pv_out[npi-1,:] = pv_out[npi-2,:]
    pv_out[:,0] = pv_out[:,1]
    pv_out[:,ny-1] = pv_out[:,ny-2]
    
    return pv_out
    
#%% Iteration
 
@njit    
def solve_PV_inversion(ly, ny, dy, dytr, pim, dpitr, pvt, pvs, dpipv, npi, \
                       pibot, pitop, dpi, dyth, thtop, dth, f, om, nit):
    
    ### Compute tropopause shape, PV distribution, and top theta
    pitp = trop_shape(ly, ny, dy, dytr, pim, dpitr)
    pv = pv_dist(pvt, pvs, dpipv, ny, npi, pibot, pitop, dpi, pitp)
    thtopline = theta_top(ly, ny, dy, dyth, thtop, dth)
    
    ### Run PV inversion iterations and apply boundary conditions
    for i in range(nit):
        f = pv_inv(pv, f, ny, pibot, dy, npi, dpi, om)
        f[npi-1,:] = f[npi-2,:] + thtopline[:]*dpi
        f[:,0] = f[:,1]
        f[:,ny-1] = f[:,ny-2]

    ### Compute u, theta, and PV given solution for phi
    u = u_calc(f, ny, npi, dy)
    theta = theta_calc(f, ny, npi, dpi)
    pv_out = pv_calc(f, ny, dy, npi, dpi, pibot)
    
    return pv, f, u, theta, pv_out

#%% Computations on eta levels
 
@njit    
def eta_calc(ny, nz, pbot, ptop, znu, npi, p_pi, u_pi, theta_pi, phi_pi):
    
    p = np.zeros((nz,ny))
    theta_eta = np.zeros((nz,ny))
    u_eta = np.zeros((nz,ny))

    p_lev = znu*(pbot - ptop) + ptop
    for j in range(ny):
        p[:,j] = p_lev
        theta_eta[:,j] = interp_decreasing(p_pi,p_lev,theta_pi[:,j],npi,nz)
        u_eta[:,j] = interp_decreasing(p_pi,p_lev,u_pi[:,j],npi,nz)
    
    rho = P0*(p/P0)**(1/GAMMA)/(RD*theta_eta)
    
    theta_surf = np.zeros(ny)
    theta_surf[:] = 1.5*theta_eta[0,:] - 0.5*theta_eta[1,:]
    rho_surf = ((pbot/P0)**(1/GAMMA))*P0/(RD*theta_surf)

    z = np.zeros((nz,ny))
    z[0,:] = -(p[0,:] - pbot)/(0.5*(rho_surf + rho[0,:])*G)
    for k in range(1,nz):
        z[k,:] = z[k-1,:] - (p[k,:] - p[k-1,:])/(0.5*(rho[k-1,:] + rho[k,:])*G)
    
    N = np.zeros((nz-1,ny))
    dtdz = np.zeros((nz-1,ny))
    for k in range(1,nz-1):
        dtdz[k,:] = (theta_eta[k+1,:] - theta_eta[k-1,:])/(z[k+1,:] - z[k-1,:])
        N[k,:] = np.sqrt(G/theta_eta[k]*dtdz[k,:])
    
    dtdz[0,:] = dtdz[1,:]
    N[0,:] = N[1,:]
    
    return p, theta_eta, u_eta, rho, z, dtdz, N

#%% PV check
 
@njit    
def pv_check(nz, ny, u, theta, z, dy, rho):

    pv_check = np.zeros((nz,ny))
    for k in range(1,nz-1):
        for j in range(1,ny-1):
            dzpv = z[k+1,j] - z[k-1,j]
            val1 = F0 - (u[k,j+1] - u[k,j-1])/(2*dy)
            val2 = (theta[k+1,j] - theta[k-1,j])/(2*dzpv)
            val3 = ((u[k+1,j] - u[k-1,j])/(2*dzpv))*((theta[k,j+1] - theta[k,j-1])/(2*dy))
            pv_check[k,j] = (val1*val2 + val3)/rho[k,j]
            
    pv_check[:,0] = pv_check[:,1]
    pv_check[:,ny-1] = pv_check[:,ny-2]
    pv_check[0,:] = pv_check[1,:]
    pv_check[nz-1,:] = pv_check[nz-2,:]

    return pv_check

#%% Faster jet with surface pressure gradient
@njit
def faster_jet(nz, ny, u, theta, p, p_mid, c, hres, p_top, z, znu):
    
    u_out = u + c
    p_surf_out = np.zeros(ny)
    py = p_mid
    
    for i in range(int(ny/2. + 1)):
        py = +F0*u_out[0,int(ny/2. - i)]*((RD*theta[0,int(ny/2. - i)]/P0)\
                      *((py/P0)**(-1./GAMMA)))**(-1)*hres*1000. + py
        p_surf_out[int(ny/2.) - i] = py

    py = p_mid
    for i in range(int(ny/2.),ny):
        py = -F0*u_out[0,i]*((RD*theta[0,i]/P0)*((py/P0)**(-1./GAMMA)))**(-1)*hres*1000. + py
        p_surf_out[i] = py

    p_out = np.zeros((nz,ny))
    for j in range(ny):
        for k in range(nz):
            p_out[k,j] = znu[k]*(p_surf_out[j] - p_top) + p_top
            
    rho_out = P0*(p_out/P0)**(1/GAMMA)/(RD*theta)

    z_out = np.zeros((nz,ny))
    z_out[0,:] = -(p_out[0,:] - p_surf_out[:])/(rho_out[0,:]*G)
    for k in range(1,nz):
        z_out[k,:] = z_out[k-1,:] - (p_out[k,:] - p_out[k-1,:])/(0.5*(rho_out[k-1,:] + rho_out[k,:])*G)
            
    return u_out, rho_out, p_out, p_surf_out, z_out

#%% Uniform barotropic shear using a linear function
@njit
def baro_shear_linear(nz, ny, u, theta, p, p_mid, c, hres, p_top, z, znu):
    
    u_out = np.zeros(np.shape(u))
    u_shear = np.zeros(ny)
    p_surf_out = np.zeros(ny)
    py = p_mid
    
    for i in range(int(ny/2. + 1)):
        u_shear[int(ny/2. - i)] = u_shear[int(ny/2. - i + 1)] - c*hres*1000.
        u_out[:,int(ny/2. - i)] = u[:,int(ny/2. - i)] + u_shear[int(ny/2. - i)]
        py = +F0*u_out[0,int(ny/2. - i)]*((RD*theta[0,int(ny/2. - i)]/P0)\
                      *((py/P0)**(-1./GAMMA)))**(-1)*hres*1000. + py
        p_surf_out[int(ny/2. - i)] = py

    py = p_mid
    for i in range(int(ny/2.),ny):
        u_shear[i] = u_shear[i-1] + c*hres*1000.
        u_out[:,i] = u[:,i] + u_shear[i]
        py = -F0*u_out[0,i]*((RD*theta[0,i]/P0)*((py/P0)**(-1./GAMMA)))**(-1)*hres*1000. + py
        p_surf_out[i] = py

    p_out = np.zeros((nz,ny))
    for j in range(ny):
        for k in range(nz):
            p_out[k,j] = znu[k]*(p_surf_out[j] - p_top) + p_top
            
    rho_out = P0*(p_out/P0)**(1/GAMMA)/(RD*theta)

    z_out = np.zeros((nz,ny))
    z_out[0,:] = -(p_out[0,:] - p_surf_out[:])/(rho_out[0,:]*G)
    for k in range(1,nz):
        z_out[k,:] = z_out[k-1,:] - (p_out[k,:] - p_out[k-1,:])/(0.5*(rho_out[k-1,:] + rho_out[k,:])*G)
            
    return u_out, rho_out, p_out, p_surf_out, z_out

#%% Uniform barotropic shear using a piecewise sine function
@njit
def baro_shear_sin(nz, ny, u, theta, p, p_mid, c, hres, p_top, z, znu, decay):
    
    u_out = np.zeros(np.shape(u))
    u_shear = np.zeros(ny)
    p_surf_out = np.zeros(ny)
    py = p_mid
    
    for i in range(ny):
        a = 1.5*(i*hres*1000. - (ny*hres*1000.)/2.)/decay
        if (np.abs(a) < np.pi/2.):
            u_shear[i] = c*np.sin(a)
        elif (a > np.pi/2.):
            u_shear[i] = c
        else:
            u_shear[i] = -c
    
    for i in range(int(ny/2. + 1)):
        u_out[:,int(ny/2. - i)] = u[:,int(ny/2. - i)] + u_shear[int(ny/2. - i)]
        py = +F0*u_out[0,int(ny/2. - i)]*((RD*theta[0,int(ny/2. - i)]/P0)\
                      *((py/P0)**(-1./GAMMA)))**(-1)*hres*1000. + py
        p_surf_out[int(ny/2. - i)] = py

    py = p_mid
    for i in range(int(ny/2.),ny):
        u_out[:,i] = u[:,i] + u_shear[i]
        py = -F0*u_out[0,i]*((RD*theta[0,i]/P0)*((py/P0)**(-1./GAMMA)))**(-1)*hres*1000. + py
        p_surf_out[i] = py

    p_out = np.zeros((nz,ny))
    for j in range(ny):
        for k in range(nz):
            p_out[k,j] = znu[k]*(p_surf_out[j] - p_top) + p_top
            
    rho_out = P0*(p_out/P0)**(1/GAMMA)/(RD*theta)

    z_out = np.zeros((nz,ny))
    z_out[0,:] = -(p_out[0,:] - p_surf_out[:])/(rho_out[0,:]*G)
    for k in range(1,nz):
        z_out[k,:] = z_out[k-1,:] - (p_out[k,:] - p_out[k-1,:])/(0.5*(rho_out[k-1,:] + rho_out[k,:])*G)
            
    return u_out, rho_out, p_out, p_surf_out, z_out

#%% Baroclinic shear using a piecewise sine function
@njit
def baroclinic_shear(nz, ny, u, theta, p, p_mid, c, hres, p_top, z, znu):
    
    u_out = np.zeros(np.shape(u))
    u_shear = np.zeros(np.shape(u))
    p_surf_out = np.zeros(ny)
    py = p_mid
    z_mid = int(nz/2)+1
    y_mid = int(ny/2)+1
    az = 2000.
    const = np.zeros(nz)
    
    for i in range(nz):
        ## Gaussian
        #const[i] = c*np.exp(-((z[i,y_mid] - z[z_mid,y_mid])/az)**2)
        # Linear
        if i < z_mid:
            const[i] = c*z[i,y_mid]/z[z_mid,y_mid]
        else:
            const[i] = c*np.exp(-((z[i,y_mid] - z[z_mid,y_mid])/az)**2)
        for j in range(ny):
            a = 1.5*(j*hres*1000. - (ny*hres*1000.)/2.)/1.5e6
            if (np.abs(a) < np.pi/2.):
                u_shear[i,j] = const[i]*np.sin(a)
            elif (a > np.pi/2.):
                u_shear[i,j] = const[i]
            else:
                u_shear[i,j] = -const[i]
    
    for i in range(int(ny/2. + 1)):
        u_out[:,int(ny/2. - i)] = u[:,int(ny/2. - i)] + u_shear[:,int(ny/2. - i)]
        py = +F0*u_out[0,int(ny/2. - i)]*((RD*theta[0,int(ny/2. - i)]/P0)\
                      *((py/P0)**(-1./GAMMA)))**(-1)*hres*1000. + py
        p_surf_out[int(ny/2. - i)] = py

    py = p_mid
    for i in range(int(ny/2.),ny):
        u_out[:,i] = u[:,i] + u_shear[:,i]
        py = -F0*u_out[0,i]*((RD*theta[0,i]/P0)*((py/P0)**(-1./GAMMA)))**(-1)*hres*1000. + py
        p_surf_out[i] = py

    p_out = np.zeros((nz,ny))
    for j in range(ny):
        for k in range(nz):
            p_out[k,j] = znu[k]*(p_surf_out[j] - p_top) + p_top
            
    rho_out = P0*(p_out/P0)**(1/GAMMA)/(RD*theta)

    z_out = np.zeros((nz,ny))
    z_out[0,:] = -(p_out[0,:] - p_surf_out[:])/(rho_out[0,:]*G)
    for k in range(1,nz):
        z_out[k,:] = z_out[k-1,:] - (p_out[k,:] - p_out[k-1,:])/(0.5*(rho_out[k-1,:] + rho_out[k,:])*G)
            
    return u_out, rho_out, p_out, p_surf_out, z_out
