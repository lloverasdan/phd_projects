#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  18 09:18:04 2022

@author: lloverasdan
"""
### Modules
import numpy as np
import netCDF4 as nc
from bwave_ideal_wrf import epv_jet
from wrf import cape_2d

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
nx = 400 # number of grid points in x direction
ny = 360 # number of grid points in y direction
nz = 100 # number of grid points in z direction
hres = 20. # horizontal grid spacing in km
zl = 20. # model top in km
p_bot = 103000. # bottom pressure for WRF in Pa
p_top = 5000. # top pressure for WRF in Pa
pibot = 1012.5 # bottom pi for EPV inversion
pitop = 424. # top pi for EPV inversion
npi_h = 180 # number of grid points in pi for high-resolution EPV inversion
npi_l = 40 # number of grid points in pi for low-resolution EPV inversion
ny_l = 30 # number of grid points in y for low-resolution EPV inversion
nit = 50000  # number of iterations for EPV inversion
om = 1.8  # successive over-relaxation coefficient for EPV inversion

### Background jet parameters
pvt = 0.2e-6 # tropospheric PV in (m^2*K)/(s*kg)
pvs = 7.0e-6 # stratospheric PV in (m^2*K)/(s*kg)
pim = 680. # mean height of tropopause in pi
dpitr = 30. # max displacement from mean tropopause height in pi
dpipv = 15. # tropopause depth in pi over which strongest PV change occurs
thtop = 535. # mean top potential temperature in K
dth = 10. # max displacement from mean top potential temperature in K
dyth = 1.0e6 # meridional scale for top potential temperature transition in m
dytr = 1.0e6  # meridional scale for tropopause height transition in m

### Horizontal shear option
shear_type = 'baro_anticyc' # choices are baro_cyc, baro_anticyc, surf_cyc, surf_anticyc

### Barotropic shear parameters
c = 10. # maximum perturbation wind speed in m/s
dyu = 1.0e6 # meridional scale for tropopause shear transition in m

### Surface-based shear parameters
dph = 750. # half displacement from bottom phi in m^2/s^2
dyph = 1.0e6 # meridional scale for bottom phi transition in m

### Moisture parameters
rh_0 = 0.75 # surface relative humidity
zrh = 8000. # height decay scale in m
delta = 1.25 # height decay parameter
nit_moist = 10 # number of iterations for virtual temperature calculation
dz = 200. # approximate z grid spacing in m

### Set up grids
ny_h = ny
ly = ny_h*hres*1000
dy_l, ylvls_l, dy_h, ylvls_h = epv_jet.y_grid(ly, ny_l, ny_h)
dpi_l, pilvls_l, dpi_h, pilvls_h, p_pi_l, p_pi_h = epv_jet.pi_grid(pibot, pitop, npi_l, npi_h)
znw, znu, dnw, rdnw, dn, rdn, fnp, fnm, cof1, cof2, cf1,\
    cf2, cf3, cfn, cfn1, dx, dy, rdx, rdy = epv_jet.eta_grid(nz, hres)

### Low-resolution run
f_l = np.zeros((npi_l,ny_l))
pv_l, f_l, u_l, theta_l, pv_out_l = epv_jet.solve_PV_inversion(ly, ny_l, \
    dy_l, dytr, pim, dpitr, pvt, pvs, dpipv, npi_l, pibot, pitop, dpi_l, \
    dyth, thtop, dth, f_l, om, nit, dyph, dph, shear_type)

### Interpolate onto high-resolution grid
f_temp = np.zeros((npi_h,ny_l))
for j in range(ny_l):
    f_temp[:,j] = epv_jet.interp_decreasing(pilvls_l,pilvls_h,f_l[:,j],npi_l,npi_h)

f_interp = np.zeros((npi_h,ny_h))
for k in range(npi_h):
    f_interp[k,:] = epv_jet.interp_increasing(ylvls_l,ylvls_h,f_temp[k,:],ny_l,ny_h)

### High-resolution run
pv_h, f_h, u_h, theta_h, pv_out_h = epv_jet.solve_PV_inversion(ly, ny_h, \
    dy_h, dytr, pim, dpitr, pvt, pvs, dpipv, npi_h, pibot, pitop, dpi_h, \
    dyth, thtop, dth, f_interp, om, nit, dyph, dph, shear_type)

### Barotropic shear
if shear_type == 'baro_cyc' or shear_type == 'baro_anticyc':
    u_h, f_h = epv_jet.barotropic_shear(u_h, f_h, ly, ny_h, npi_h, dy_h, dpi_h, c, dyu, shear_type)

### Unstaggered pressure for WRF levels
p_unstag = znu*(p_bot - p_top) + p_top

### Interpolate from pi to WRF levels
u_eta = np.zeros((nz,ny))
theta_eta = np.zeros((nz,ny))
ph_eta = np.zeros((nz,ny))
p_eta = np.zeros((nz,ny))
pv_eta = np.zeros((nz,ny))
for j in range(ny):
    u_eta[:,j] = epv_jet.interp_decreasing(p_pi_h,p_unstag,u_h[:,j],npi_h,nz)
    theta_eta[:,j] = epv_jet.interp_decreasing(p_pi_h,p_unstag,theta_h[:,j],npi_h,nz)
    ph_eta[:,j] = epv_jet.interp_decreasing(p_pi_h,p_unstag,f_h[:,j],npi_h,nz)
    pv_eta[:,j] = epv_jet.interp_decreasing(p_pi_h,p_unstag,pv_h[:,j],npi_h,nz)
    p_eta[:,j] = p_unstag

### Moisture
rh = np.zeros((nz,ny))
for k in range(nz):
    z = k*dz
    if z < zrh: rh[k,:] = np.ones(ny)*rh_0*(1 - 0.9*(z/zrh)**delta)
    else: rh[k,:] = np.ones(ny)*rh_0*0.1

theta_new = np.zeros((nit_moist,nz,ny))
qv = np.zeros((nit_moist-1,nz,ny))
theta_new[0,:,:] = theta_eta
for it in range(nit_moist-1):
    temp = theta_new[it,:,:]*(p_eta/P0)**(RD/CP)
    es = 611.2*np.exp(17.67*(temp - SVPT0)/(temp - 29.65))
    qv[it,:,:] = rh*(RD/RV)*(es/(p_eta - es))
    theta_new[it+1,:,:] = theta_eta/(1 + (RV/RD)*qv[it,:,:])
    
rh_it = np.zeros((nit_moist-1,nz,ny))
for it in range(nit_moist-1):
    temp_it = theta_new[it+1,:,:]*(p_eta/P0)**KAPPA
    es_it = 611.2*np.exp(17.67*(temp_it - SVPT0)/(temp_it - 29.65))
    qs_it = (RD/RV)*(es_it/(p_eta - es_it))
    rh_it[it,:,:] = qv[it,:,:]/qs_it*100.  

### CAPE
p_full = np.zeros((nz,ny,nx))
t_full = np.zeros((nz,ny,nx))
qv_full = np.zeros((nz,ny,nx))
z_full = np.zeros((nz,ny,nx))
for i in range(nx):
    p_full[:,:,i] = p_eta/100.
    t_full[:,:,i] = theta_new[-1,:,:]*(p_eta/P0)**KAPPA
    qv_full[:,:,i] = qv[-1,:,:]
    z_full[:,:,i] = ph_eta/G

terr = np.zeros((ny, nx))
ter_follow = False
cape = cape_2d(p_full,t_full,qv_full,z_full,terr,p_full[0,:,:],ter_follow)
cape1d = np.nan_to_num(cape[0,:,0])

### Save files
np.save('/p/work1/lloveras/bwave/processed/jets/ubot_barotropic_rh75',u_h[0])
np.save('/p/work1/lloveras/bwave/processed/jets/phbot_barotropic_rh75',f_h[0])
np.save('/p/work1/lloveras/bwave/processed/jets/u_barotropic_rh75',u_eta)
np.save('/p/work1/lloveras/bwave/processed/jets/th_barotropic_rh75',theta_new)
np.save('/p/work1/lloveras/bwave/processed/jets/thv_barotropic_rh75',theta_eta)
np.save('/p/work1/lloveras/bwave/processed/jets/rh_barotropic_rh75',rh_it)
np.save('/p/work1/lloveras/bwave/processed/jets/qv_barotropic_rh75',qv)
np.save('/p/work1/lloveras/bwave/processed/jets/pv_barotropic_rh75',pv_eta)
np.save('/p/work1/lloveras/bwave/processed/jets/z_barotropic_rh75',ph_eta/G)
np.save('/p/work1/lloveras/bwave/processed/jets/p_barotropic_rh75',p_eta)
np.save('/p/work1/lloveras/bwave/processed/jets/cape_barotropic_rh75',cape1d)
