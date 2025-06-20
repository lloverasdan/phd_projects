#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  18 09:18:04 2022

@author: lloverasdan
"""

import numpy as np
from netCDF4 import Dataset
from bwave_ideal_wrf import epv_jet, qgpv_pert
from wrf import destagger, getvar

### Constants
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

### Grid parameters
nx = 250 # number of grid points in x direction
ny = 225 # number of grid points in y direction
nz = 100 # number of grid points in z direction
hres = 32. # horizontal grid resolution in km
zl = 20. # model top in km
dz = 200. # approximate z grid spacing in m
ly = hres*ny*1000. # length in y direction in m

### Tropopause anomaly parameters
up_pert = False # tropopause anomaly option
up_pert_mag = 1.0e-4 # magnitude of tropopause anomaly in s^-1
x_up_pert = 94 #150 # x center gridpoint for tropopause anomaly
y_up_pert = 112 #180 # y center gridpoint for tropopause anomaly
z_up_pert = 50  # z center gridpoint for tropopause anomaly
ax_up_pert = 200. #700. #200. # x decay scale of tropopause anomaly in km
ay_up_pert = 600. #700. #600. # y decay scale of tropopause anomaly in km
az_up_pert = 1.5 # z decay scale of tropopause anomaly in km

### Surface anomaly parameters
surf_pert = True # surface anomaly option
surf_pert_mag = 5.0 # magnitude of surface theta anomaly in K
x_surf_pert = 130 # x center gridpoint for surface anomaly
y_surf_pert = 85 # y center gridpoint for surface anomaly
ax_surf_pert = 600. #700. #600. # x decay scale for surface anomaly in km
ay_surf_pert = 200. #700. #200. # y decay scale for surface anomaly in km

### Moisture parameters
moist = True # moisture option, Yes=True, No=False
rh_0 = 0.85 # reference relative humidity
zrh = 8000. # height decay scale in m
delta = 1.25 # height decay parameter
nit = 10 # number of iterations for regula falsi
# dy = 1e6 # length scale for meridional relative humidity transition
# rm = 0.9 # center relative humidity
# dr = 0.1 # difference in relative humidity

### Read in the file
nc_in = Dataset('/p/work/lloveras/bwave/32km_files/anticyc_new/input/wrfin_moist_surfonly','r+')

### Base variables
znu = np.asarray(nc_in.variables['ZNU'][0,:])
znw = np.asarray(nc_in.variables['ZNW'][0,:])
dnw = np.asarray(nc_in.variables['DNW'][0,:])
dn = np.asarray(nc_in.variables['DN'][0,:])
phb = np.asarray(nc_in.variables['PHB'][0,:,:,:])
pb = np.asarray(nc_in.variables['PB'][0,:,:,:])
mub = np.asarray(nc_in.variables['MUB'][0,:,:])
t_init = np.asarray(nc_in.variables['T_INIT'][0,:,:])
alb = (RD/P0)*(t_init + T0)*((pb/P0)**(-1/GAMMA))

### Full variables
u_jet = destagger(np.asarray(nc_in.variables['U'][0,:,:,:]),-1)
v_jet = destagger(np.asarray(nc_in.variables['V'][0,:,:,:]),-2)
th_jet = np.asarray(nc_in.variables['T'][0,:,:,:]) + T0
ph_jet = np.asarray(nc_in.variables['PH'][0,:,:,:]) + phb
p_jet = np.asarray(nc_in.variables['P'][0,:,:,:]) + pb
mu_jet = np.asarray(nc_in.variables['MU'][0,:,:]) + mub
z_jet = destagger(ph_jet/G,0)
rho_jet = P0/(RD*th_jet)*(p_jet/P0)**(CV/CP)

### Set up grids
xl, yl, x, y, xg, yg, dz, zp, facz = qgpv_pert.cartesian_mesh(nx, ny, nz, hres, zl)
kmax, lmax, facx, facy, dxg, dyg = qgpv_pert.spectral_mesh(nx, ny, xl, yl)

### Compute reference static stability
n_jet = np.zeros((nz-1,ny,nx))
dtdz_jet = np.zeros((nz-1,ny,nx))
for k in range(1,nz-1):
    dtdz_jet[k,:,:] = (th_jet[k+1,:,:] - th_jet[k-1,:,:])/(z_jet[k+1,:,:] - z_jet[k-1,:,:])
    n_jet[k,:,:] = np.sqrt(G/th_jet[k,:,:]*dtdz_jet[k,:,:])
    
dtdz_jet[0,:,:] = dtdz_jet[1,:,:]
n_jet[0,:,:] = n_jet[1,:,:]

dtdz_interp = np.zeros((nz-1,ny,nx))
n_interp = np.zeros((nz-1,ny,nx))
for i in range(nx):
    for j in range(ny):
        for k in range(nz-1):
            dtdz_interp[k,j,i] = qgpv_pert.interp_0(dtdz_jet[:,j,i], z_jet[:-1,j,i], zp[k], nz)
            n_interp[k,j,i] = qgpv_pert.interp_0(n_jet[:,j,i], z_jet[:-1,j,i], zp[k], nz)

dtdz_z = np.mean(dtdz_interp,axis=-1)
n_z = np.mean(n_interp,axis=-1)

# n_up_pert = n_z[:,y_up_pert]
# n_surf_pert = n_z[:,y_surf_pert]
# dtdz_up_pert = dtdz_z[:,y_up_pert]
# dtdz_surf_pert = dtdz_z[:,y_surf_pert]

n_up_pert = n_z[:,0]
n_surf_pert = n_z[:,0]
dtdz_up_pert = dtdz_z[:,0]
dtdz_surf_pert = dtdz_z[:,0]

### Tropopause anomaly
if up_pert:
    
    ### Initialize
    pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac = qgpv_pert.trop_anom(up_pert_mag,\
        x_up_pert, y_up_pert, z_up_pert, az_up_pert, ax_up_pert, ay_up_pert, x, y, zp, xg, yg, nx, ny, \
        nz, n_up_pert, facz, dz, dtdz_up_pert)
    
    ### Invert
    fbsp, ftsp, fzbsp, fztsp, fsp = qgpv_pert.qgpv_inversion(nx, ny, nz, bu_fac, \
            facx, facy, facz, kmax, lmax, pvsp, ubcsp, lbcsp, dz)

    ### Compute perturbations
    u_up_pert, v_up_pert, theta_up_pert, rho_up_pert, fxy_up = qgpv_pert.qgpv_solver(fsp, nx, ny, nz, \
            bu_fac, dxg, dyg, dz, lbcxy, ubcxy, dtdz_up_pert)

    p_up_pert = fxy_up*F0
    
else:
    u_up_pert = np.zeros((nz,ny,nx))
    v_up_pert = np.zeros((nz,ny,nx))
    theta_up_pert = np.zeros((nz,ny,nx))
    rho_up_pert = np.zeros((nz,ny,nx))
    p_up_pert = np.zeros((nz,ny,nx))

### Surface anomaly
if surf_pert:
    
    ### Initialize
    pvs, ubcs, lbcs, pvsps, ubcsps, lbcsps, bu_facs = qgpv_pert.surf_anom(surf_pert_mag,\
        x_surf_pert, y_surf_pert, ax_surf_pert, ay_surf_pert, x, y, zp, xg, yg, nx, ny, nz,\
        n_surf_pert, facz, dz, dtdz_surf_pert)
    
    ### Invert
    fbsps, ftsps, fzbsps, fztsps, fsps = qgpv_pert.qgpv_inversion(nx, ny, nz, bu_facs, \
        facx, facy, facz, kmax, lmax, pvsps, ubcsps, lbcsps, dz)
    
    ### Compute perturbations
    u_surf_pert, v_surf_pert, theta_surf_pert, rho_surf_pert, fxy_surf = qgpv_pert.qgpv_solver(fsps, nx, ny, nz, \
        bu_facs, dxg, dyg, dz, lbcs, ubcs, dtdz_surf_pert)

    p_surf_pert = fxy_surf*F0

else:
    u_surf_pert = np.zeros((nz,ny,nx))
    v_surf_pert = np.zeros((nz,ny,nx))
    theta_surf_pert = np.zeros((nz,ny,nx))
    rho_surf_pert = np.zeros((nz,ny,nx))
    p_surf_pert = np.zeros((nz,ny,nx))

### Add perturbations together
u_pert = u_up_pert + u_surf_pert
v_pert = v_up_pert + v_surf_pert
theta_pert = theta_up_pert + theta_surf_pert
rho_pert = rho_up_pert + rho_surf_pert
p_pert = p_up_pert + p_surf_pert

### Interpolate onto eta levels
up_eta = np.zeros((nz,ny,nx))
vp_eta = np.zeros((nz,ny,nx))
thetap_eta = np.zeros((nz,ny,nx))
rhop_eta = np.zeros((nz,ny,nx))
presp_eta = np.zeros((nz,ny,nx))
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            thetap_eta[k,j,i] = qgpv_pert.interp_0(theta_pert[:,j,i], zp, z_jet[k,j,i], nz)
            rhop_eta[k,j,i] = qgpv_pert.interp_0(rho_pert[:,j,i], zp, z_jet[k,j,i], nz)            
            up_eta[k,j,i] = qgpv_pert.interp_0(u_pert[:,j,i], zp, z_jet[k,j,i], nz)
            vp_eta[k,j,i] = qgpv_pert.interp_0(v_pert[:,j,i], zp, z_jet[k,j,i], nz)
            presp_eta[k,j,i] = qgpv_pert.interp_0(p_pert[:,j,i], zp, z_jet[k,j,i], nz)            

u_eta = u_jet + up_eta
v_eta = v_jet + vp_eta
rho_eta = rho_jet + rhop_eta
theta_eta = th_jet + thetap_eta
p_eta = p_jet + presp_eta
ph_eta = z_jet*G

### Compute new perturbation variables
p_bot = 1.5*p_eta[0,:,:] - 0.5*p_eta[1,:,:]
p_top = 1.5*p_eta[-1,:,:] - 0.5*p_eta[-2,:,:]
u, v, ph, p, mu, t, tsk = epv_jet.pert_wrf(u_eta, v_eta, ph_eta, p_eta, theta_eta, pb, alb, mub,\
                                           p_bot, p_top, nx, ny, nz, znu, znw, dn, dnw)
        
### Add moisture and update theta
if moist:
    rh = np.zeros((nz,ny,nx))
#     r = np.zeros(ny)
    for i in range(nx):
        for j in range(ny):
#             a1 = 2.*(j*hres*1000. - ly/2.)/dy
#             if (np.abs(a1) <= np.pi/2.): r[j] = rm + dr*np.sin(a1)
#             elif (a1 > np.pi/2.): r[j] = rm + dr
#             else: r[j] = rm - dr
            for k in range(nz):
                z = k*dz
#                 if z < zrh: rh[k,j,i] = rh_0*r[j]*(1 - 0.9*(z/zrh)**delta)
#                 else: rh[k,j,i] = rh_0*r[j]*0.1
                if z < zrh: rh[k,j,i] = rh_0*(1 - 0.9*(z/zrh)**delta)
                else: rh[k,j,i] = rh_0*0.1

    theta_new = np.zeros((nit,nz,ny,nx))
    theta_new[0,:,:,:] = t + T0
    for it in range(nit-1):
        temp = theta_new[it,:,:,:]*((p + pb)/P0)**(RD/CP)
        es = 611.2*np.exp(17.67*(temp - SVPT0)/(temp - 29.65))
        qv = rh*(RD/RV)*(es/((p + pb) - es))
        theta_new[it+1,:,:,:] = (t + T0)/(1 + RV/RD*qv)

    t_moist = theta_new[nit-1,:,:,:] - T0
    tsk_moist = theta_new[nit-1,0,:,:]*(((p + pb)[0,:,:]/P0)**(RD/CP))
    
else:
    t_moist = t
    tsk_moist = tsk
    qv = np.zeros((nz,ny,nx))

### Write new variables to the input file
nc_in.variables['U'][0,:,:,:] = u
nc_in.variables['V'][0,:,:,:] = v
nc_in.variables['PH'][0,:,:,:] = ph
nc_in.variables['MU'][0,:,:] = mu
nc_in.variables['P'][0,:,:,:] = p
nc_in.variables['T'][0,:,:,:] = t_moist
nc_in.variables['TSK'][0,:,:] = tsk_moist
nc_in.variables['QVAPOR'][0,:,:,:] = qv
nc_in.close()
