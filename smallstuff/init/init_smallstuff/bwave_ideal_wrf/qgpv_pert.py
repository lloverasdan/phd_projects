#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:40:27 2020

@author: Daniel Lloveras

Defines functions for inverting a 3D QGPV anomaly to trigger baroclinic growth
"""
import numpy as np
from numba import njit

#%% Constants

F0 = 1.0e-4  # Coriolis parameter, in s^-1
G = 9.81  # gravitational acceleration, in m/s^2

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

#%% Grids

def cartesian_mesh(nx, ny, nz, hres, zl):
    
    xl = nx*hres*1000.  # x domain length, in m
    yl = ny*hres*1000. # y domain length, in m
    ddx = xl/nx  # grid spacing in x direction
    ddy = yl/ny  # grid spacing in y direction
    xx = np.arange(0,xl+ddx,ddx)  # unstaggered grid in x direction
    x = xx[1:nx+1] - xl/2.  # staggered grid in x direction 
    yy = np.arange(0,yl+ddy,ddy)  # unstaggered grid in y direction
    y = yy[1:ny+1] - yl/2.  # staggered grid in y direction
    xg,yg = np.meshgrid(x,y)  # staggered grid mesh
    dz = zl*1000./nz  # grid spacing in z direction
    z = np.arange(1,nz+1)*dz - (dz/2) # staggered grid in z direction
    facz = np.ones(nz)

    return xl, yl, x, y, xg, yg, dz, z, facz
    
def spectral_mesh(nx, ny, xl, yl):
    
    kmax = nx/2.  # number of x waves
    lmax = ny/2.  # number of y waves
    facx = 2.*np.pi/xl  # Fourier factor in x direction
    facy = 2.*np.pi/yl  # Fourier factor in y direction
    dx = np.arange(-kmax,kmax)*facx  # spectral grid in x direction
    dy = np.arange(-lmax,lmax)*facy  # spectral grid in y direction
    dxg,dyg = np.meshgrid(dx,dy)  # spectral grid mesh
    dxg = np.fft.fftshift(dxg)  # shift zero-frequency component to center
    dyg = np.fft.fftshift(dyg)  # shift zero-frequency component to center
    
    return kmax, lmax, facx, facy, dxg, dyg

#%% Original QGPV anomaly
    
def qgpv_anom_gauss(qgpv_mag, theta_mag, x_pert, y_pert, z_pert, az, ax, ay, x, y, z, \
              xg, yg, nx, ny, nz, n_pert, facz, dz, dtdz_pert):
    
    ax = ax*1000.  # horizontal decay scale, in m
    ay = ay*1000.
    az = az*1000.  # vertical decay scale, m
    
    xnot = x[x_pert]
    ynot = y[y_pert]
    
    rr = np.sqrt(((xg - xnot)/ax)**2 + ((yg - ynot)/ay)**2)
    ubcxy = theta_mag*np.exp(-(rr)**2)*np.exp(-((z[nz-1] - z[z_pert-1])/az)**2)
    lbcxy = theta_mag*np.exp(-(rr)**2)*np.exp(-((z[0] - z[z_pert-1])/az)**2)
    
    pvxy = np.zeros((nz,ny,nx))
    for k in range(nz):
        pvxy[k,:,:] = qgpv_mag*np.exp(-(rr)**2)*np.exp(-((z[k-1] - z[z_pert-1])/az)**2)

    ### Volume-integrated PV must be zero
    pvxy = pvxy - np.mean(pvxy)
    
    ### Coefficents in vertical direction
    bu = n_pert/F0
    bu_fac = bu**2
    
    ### Transform PV anomaly to spectral space in horizontal
    lbcsp = np.fft.fft2(lbcxy)
    ubcsp = np.fft.fft2(ubcxy)
    pvsp = np.zeros((nz,ny,nx),dtype=np.dtype(np.complex64))
    for k in range(0,nz):
        pvsp[k,:,:] = np.fft.fft2(np.squeeze(pvxy[k,:,:]))
    
    ### Apply Neumann boundary conditions    
    pvsp[0,:,:] = np.squeeze(pvsp[0,:,:]) + F0*(lbcsp/(facz[0]*dz*dtdz_pert[0]))
    pvsp[nz-1,:,:] = np.squeeze(pvsp[nz-1,:,:]) - F0*(ubcsp/(facz[nz-1]*dz*dtdz_pert[nz-2]))
    
    return pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac

#%% Cosine QGPV anomaly
    
def qgpv_anom_cos(qgpv_mag, theta_mag, x_pert, y_pert, z_pert, az, ax, ay, x, y, z, \
              xg, yg, nx, ny, nz, n_pert, facz, dz, dtdz_pert):
    
    ax = ax*1000.  # horizontal decay scale, in m
    ay = ay*1000.
    az = az*1000.  # vertical decay scale, m
    
    xnot = x[x_pert]
    ynot = y[y_pert]
    
    rr_xy = np.sqrt(((xg - xnot)/ax)**2 + ((yg - ynot)/ay)**2)   
    rr_xy[rr_xy > np.pi/2.] = np.pi/2.
    
    rr_zt = np.abs((z[nz-1] - z[z_pert-1])/az)
    if rr_zt > np.pi/2.: rr_zt = np.pi/2.
    rr_zb = np.abs((z[0] - z[z_pert-1])/az)
    if rr_zb > np.pi/2.: rr_zb = np.pi/2.
    
    ubcxy = theta_mag*np.cos(-rr_xy)*np.cos(-rr_zt)
    lbcxy = theta_mag*np.cos(-rr_xy)*np.cos(-rr_zb)
    
    pvxy = np.zeros((nz,ny,nx))
    for k in range(nz):
        rr_z = np.abs((z[k-1] - z[z_pert-1])/az)
        if rr_z > np.pi/2.: rr_z = np.pi/2.
        pvxy[k,:,:] = qgpv_mag*np.cos(-rr_xy)*np.cos(-rr_z)

    ### Volume-integrated PV must be zero
    pvxy = pvxy - np.mean(pvxy)
    
    ### Coefficents in vertical direction
    bu = n_pert/F0
    bu_fac = bu**2
    
    ### Transform PV anomaly to spectral space in horizontal
    lbcsp = np.fft.fft2(lbcxy)
    ubcsp = np.fft.fft2(ubcxy)
    pvsp = np.zeros((nz,ny,nx),dtype=np.dtype(np.complex64))
    for k in range(0,nz):
        pvsp[k,:,:] = np.fft.fft2(np.squeeze(pvxy[k,:,:]))
    
    ### Apply Neumann boundary conditions    
    pvsp[0,:,:] = np.squeeze(pvsp[0,:,:]) + F0*(lbcsp/(facz[0]*dz*dtdz_pert[0]))
    pvsp[nz-1,:,:] = np.squeeze(pvsp[nz-1,:,:]) - F0*(ubcsp/(facz[nz-1]*dz*dtdz_pert[nz-2]))
    
    return pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac

#%% Theta QGPV anomaly
    
def qgpv_theta(qgpv_mag, theta_mag, x_pert, y_pert, th_pert, ath, ax, ay, x, y, z, \
              xg, yg, nx, ny, nz, n_pert, facz, dz, theta_jet, z_jet, dtdz_pert):
    
    ### X-Y Plane
    ax = ax*1000.  # horizontal decay scale, in m
    ay = ay*1000.
    xnot = x[x_pert]
    ynot = y[y_pert]
    rr_xy = np.sqrt(((xg - xnot)/ax)**2 + ((yg - ynot)/ay)**2)
    rr_xy[rr_xy > np.pi/2.] = np.pi/2.
    
    ### Theta Plane and Boundaries
    th = np.linspace(np.min(theta_jet),np.max(theta_jet),nz)
    th_ind = min(range(len(th)), key=lambda i: abs(th[i]-th_pert))
    rr_tt = np.abs((th[nz-1] - th[th_ind])/ath)
    if rr_tt > np.pi/2.: rr_tt = np.pi/2.
    rr_tb = np.abs((th[0] - th[th_ind])/ath)
    if rr_tb > np.pi/2.: rr_tb = np.pi/2.
    
    ubcxy = theta_mag*np.cos(-rr_xy)*np.cos(-rr_tt)
    lbcxy = theta_mag*np.cos(-rr_xy)*np.cos(-rr_tb)
    
    ### Function
    pvxyt = np.zeros((nz,ny,nx))
    for k in range(nz):
        rr_th = np.abs((th[k] - th[th_ind])/ath)
        if rr_th > np.pi/2.: rr_th = np.pi/2.
        pvxyt[k,:,:] = qgpv_mag*np.cos(-rr_xy)*np.cos(-rr_th)

    ### Interpolate onto z coordinates
    theta_z = np.zeros((nz,ny))
    for j in range(ny):
        for k in range(nz):
            theta_z[k,j] = interp_0(theta_jet[:,j],z_jet[:,j],z[k],nz)
            
    pvxy = np.zeros((nz,ny,nx))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                pvxy[k,j,i] = interp_0(pvxyt[:,j,i],th,theta_z[k,j],nz)

    ### Volume-integrated PV must be zero
    pvxy = pvxy - np.mean(pvxy)
    
    ### Coefficents in vertical direction
    bu = n_pert/F0
    bu_fac = bu**2
    
    ### Transform PV anomaly to spectral space in horizontal
    lbcsp = np.fft.fft2(lbcxy)
    ubcsp = np.fft.fft2(ubcxy)
    pvsp = np.zeros((nz,ny,nx),dtype=np.dtype(np.complex64))
    for k in range(0,nz):
        pvsp[k,:,:] = np.fft.fft2(np.squeeze(pvxy[k,:,:]))
    
    ### Apply Neumann boundary conditions    
    pvsp[0,:,:] = np.squeeze(pvsp[0,:,:]) + F0*(lbcsp/(facz[0]*dz*dtdz_pert[0]))
    pvsp[nz-1,:,:] = np.squeeze(pvsp[nz-1,:,:]) - F0*(ubcsp/(facz[nz-1]*dz*dtdz_pert[nz-2]))
    
    return pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac

def qgpv_surf(qgpv_mag, theta_mag, x_pert, y_pert, az, ax, ayn, ays, x, y, z, \
              xg, yg, nx, ny, nz, n_pert, facz, dz, dtdz_pert):
    
    ax = ax*1000.  # horizontal decay scale, in m
    ayn = ayn*1000.
    ays = ays*1000.
    az = az*1000.  # vertical decay scale, m
    
    xnot = x[x_pert]
    ynot = y[y_pert]
    
    rr_xy = np.zeros((ny,nx))
    for i in range(nx):
        for j in range(ny):
            if j > y_pert: ay = ayn
            else: ay = ays
            rr_xy[j,i] = np.sqrt(((xg[j,i] - xnot)/ax)**2 + ((yg[j,i] - ynot)/ay)**2)
    
    rr_xy[rr_xy > np.pi/2.] = np.pi/2.
    
    ubcxy = np.zeros((ny,nx))
    lbcxy = theta_mag*np.cos(-rr_xy)
    
    pvxy = np.zeros((nz,ny,nx))
    for k in range(nz):
        pvxy[k,:,:] = qgpv_mag*np.cos(-rr_xy)*np.exp(-z[k]/az)

    ### Volume-integrated PV must be zero
    pvxy = pvxy - np.mean(pvxy)
    
    ### Coefficents in vertical direction
    bu = n_pert/F0
    bu_fac = bu**2
    
    ### Transform PV anomaly to spectral space in horizontal
    lbcsp = np.fft.fft2(lbcxy)
    ubcsp = np.fft.fft2(ubcxy)
    pvsp = np.zeros((nz,ny,nx),dtype=np.dtype(np.complex64))
    for k in range(0,nz):
        pvsp[k,:,:] = np.fft.fft2(np.squeeze(pvxy[k,:,:]))
    
    ### Apply Neumann boundary conditions    
    pvsp[0,:,:] = np.squeeze(pvsp[0,:,:]) + F0*(lbcsp/(facz[0]*dz*dtdz_pert[0]))
    pvsp[nz-1,:,:] = np.squeeze(pvsp[nz-1,:,:]) - F0*(ubcsp/(facz[nz-1]*dz*dtdz_pert[nz-2]))
    
    return pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac


#%% QGPV inversion

@njit
def qgpv_inversion(nx, ny, nz, bu_fac, facx, facy, facz, kmax, lmax, pvsp, \
                   ubcsp, lbcsp, dz):
    
    fsp = np.zeros((nz,ny,nx),dtype=np.dtype(np.complex64))  # spectral phi     
    fbsp = np.zeros((ny,nx),dtype=np.dtype(np.complex64))    # spectral phi on bottom boundary
    ftsp = np.zeros((ny,nx),dtype=np.dtype(np.complex64))    # spectral phi on top boundary 
    fzbsp = np.zeros((ny,nx),dtype=np.dtype(np.complex64))   # spectral d(phi)/dz on bottom boundary
    fztsp = np.zeros((ny,nx),dtype=np.dtype(np.complex64))   # spectral d(phi)/dz on top boundary
    for k in range(0,nx):
        for l in range(0,ny):
            ### Compute x and y wavenumbers
            ak = facx*k  # lower left quadrant
            bl = facy*l  # lower left quadrant
            if ((k+1) >= kmax and (l+1) <= lmax):  # lower right quadrant
                ak = -facx*(nx - k)
            elif ((k+1) <= kmax and (l+1) >= lmax):  # upper left quadrant
                bl = -facy*(ny - l)
            elif ((k+1) >= kmax and (l+1) >= lmax):  # upper right quadrant
                ak = -facx*(nx - k)
                bl = -facy*(ny - l)
                            
            ### Set phi in bottom left corner to 0
            if (k == 0 and l == 0):
                fbsp[l,k] = 0.
                ftsp[l,k] = 0.
                fzbsp[l,k] = 0.
                fztsp[l,k] = 0.
                fsp[:,l,k] = 0.            
            ### Construct tridiagonal matrix and invert it    
            else:
                d1 = np.zeros(nz)  # main diagonal factor 
                u1 = np.zeros(len(bu_fac))  # upper diagonal factor
                l1 = np.zeros(len(bu_fac))  # lower diagonal factor
                for al in range(0,nz-1):
                    d1[al] = 1./bu_fac[al]
                    u1[al] = d1[al]
            
                d1[nz-1] = d1[nz - 2]
            
                for al in range(0,nz-2):
                    l1[al] = 1./bu_fac[al]
                    d1[al+1] = d1[al+1] + 1./bu_fac[al]
            
                l1[nz-2] = 1./bu_fac[nz-2]
    
                dd = np.diag(-(ak**2 + bl**2 + d1*(1/(dz**2))))  # main diagonal
                ud = np.diag(u1*(1/(dz**2)),k=1)  # upper diagonal
                ld = np.diag(l1*(1/(dz**2)),k=-1)  # lower diagonal
                Atemp = dd + ud +ld
                
                ### Fix upper left, lower right corners to satisfy Neumann boundary conditions
                Atemp[0,0] = -(ak**2 + bl**2 + (1./bu_fac[0])*(1./(dz**2)))
                Atemp[-1,-1] = -(ak**2 + bl**2 + (1./bu_fac[nz-2])*(1./(dz**2)))
                
                ### Invert the matrix
                A = Atemp.astype(np.dtype(np.complex64))
                psi = np.linalg.solve(A,pvsp[:,l,k])
                fsp[:,l,k] = psi[:]
                
                ### Apply Neumann boundary conditions
                fbsp[l,k] = psi[0] - (0.5*dz*facz[0]*lbcsp[l,k])
                ftsp[l,k] = psi[nz-1] + (0.5*dz*facz[nz-1]*ubcsp[l,k])
                fzbsp[l,k] = lbcsp[l,k]
                fztsp[l,k] = ubcsp[l,k]
                
    return fbsp, ftsp, fzbsp, fztsp, fsp                
    
#%% Solve for perturbation fields
                
def qgpv_solver(fsp, nx, ny, nz, bu_fac, dxg, dyg, dz, lbcxy, ubcxy, dtdz_pert):               
    
    fxy = np.zeros((nz,ny,nx))
    for k in range(0,nz):
        fxy[k,:,:] = np.real(np.fft.ifft2(np.squeeze(fsp[k,:,:])))
        
    u_pert = np.zeros(np.shape(fxy))
    v_pert = np.zeros(np.shape(fxy))
    theta_pert = np.zeros(np.shape(fxy))
    rho_pert = np.zeros(np.shape(fxy))
    
    for k in range(0,nz):
        u_pert[k,:,:] = np.real(np.fft.ifft2((-1j*dyg)*np.squeeze(fsp[k,:,:])))
        v_pert[k,:,:] = np.real(np.fft.ifft2((1j*dxg)*np.squeeze(fsp[k,:,:])))
        
    for k in range(1,nz-1):
        rho_pert[k,:,:] = -F0*(fxy[k+1,:,:] - fxy[k-1,:,:])/(2*G*dz)            
        
    rho_pert[nz-1,:,:] = -(F0*(fxy[nz-1,:,:] - fxy[nz-2,:,:])/(G*dz))
    rho_pert[0,:,:] = -(F0*(fxy[1,:,:] - fxy[0,:,:])/(G*dz))
    
    for k in range(1,nz-1):
        theta_pert[k,:,:] = dtdz_pert[k]*(1/bu_fac[k])*\
                            (fxy[k+1,:,:] - fxy[k-1,:,:])/(2*dz*F0)  
    
    theta_pert[0,:,:] = lbcxy
    theta_pert[nz-1,:,:] = ubcxy
    
    return u_pert, v_pert, theta_pert, rho_pert, fxy


#%% QGPV Check
    
def qgpv_check(fbsp, ftsp, fzbsp, fztsp, fsp, nx, ny, nz, dz, dxg, dyg, bu_fac):
    
    fxy = np.zeros((nz,ny,nx))
    for k in range(0,nz):
        fxy[k,:,:] = np.real(np.fft.ifft2(np.squeeze(fsp[k,:,:])))
    
    ### Compute d(phi)/dz in spectral space
    fz = np.zeros((nz-1,ny,nx),dtype=np.dtype(np.complex64))    
    for k in range(0,nz-1):
        fz[k,:,:] = (fsp[k+1,:,:] - fsp[k,:,:])/dz
    
    ### Compute PV in spectral space given the phi field
    pvchecksp = np.zeros((nz,ny,nx),dtype=np.dtype(np.complex64))
    for k in range(0,nz):
        fxx = ((1j*dxg)**2)*np.squeeze(fsp[k,:,:]) 
        fyy = ((1j*dyg)**2)*np.squeeze(fsp[k,:,:])
        if (k==0):
            fzz = ((1./bu_fac[0])*(np.squeeze(fsp[k+1,:,:]) - np.squeeze(fsp[k,:,:]))/dz \
                    - (1./bu_fac[0])*fzbsp)/dz
        elif (k==nz-1):
            fzz = ((1/bu_fac[k-1])*fztsp - (1/bu_fac[k-1])*(np.squeeze(fsp[k,:,:]) \
                    - np.squeeze(fsp[k-1,:,:]))/dz)/dz
        else:
             fzz = np.squeeze(((1/bu_fac[k])*fsp[k+1,:,:] - (1/bu_fac[k] + 1/bu_fac[k-1])*fsp[k,:,:] \
                               + (1/bu_fac[k-1])*fsp[k-1,:,:])/dz)/dz     
        
        fzz[0,0] = 0.
        pvchecksp[k,:,:] = fxx + fyy + fzz
    
    ### Transfrom the PV check into cartesian space
    pvcheckxy = np.zeros((nz,ny,nx))
    for k in range(0,nz):
        pvcheckxy[k,:,:] = np.real(np.fft.ifft2(np.squeeze(pvchecksp[k,:,:])))
    
    return pvcheckxy
