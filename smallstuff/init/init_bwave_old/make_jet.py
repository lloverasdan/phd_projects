#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  18 09:18:04 2022

@author: lloverasdan
"""
### Modules
import numpy as np
import netCDF4 as nc
from bwave_ideal_wrf import epv_jet, write_wrfinputd01

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

### WRF file strings
file_name = '/p/work/lloveras/bwave/32km_files/scratch/setup/wrfin_start'
netcdf_type = 'NETCDF3_64BIT_OFFSET'
title_str = "OUTPUT FROM IDEAL V3.6.1 PREPROCESSOR"
time_str = "2021-01-01_00:00:00"

### Grid parameters
nx = 250 # number of grid points in x direction
ny = 225 # number of grid points in y direction
nz = 100 # number of grid points in z direction
hres = 32. # horizontal grid spacing in km
zl = 20. # model top in km
p_bot = 101200. #100800. or 101200. # bottom pressure in Pa
p_top = 5000. # top pressure in Pa
pibot = 1004. # bottom height for EPV inversion in pi
pitop = 424. # top height for EPV inversion in pi
npi_h = 180 # number of grid points in pi for high-resolution EPV inversion
npi_l = 40 # number of grid points in pi for low-resolution EPV inversion
ny_l = 30 # number of grid points in y for low-resolution EPV inversion
nit = 50000  # number of iterations for EPV inversion
om = 1.8  # successive over-relaxation coefficient for EPV inversion

### Background jet parameters
pv_dist = 'ctss' # PV distribution option, ctss or 2lpv
pvt = 0.2e-6 # tropospheric PV in (m^2*K)/(s*kg), CTSS=0.2e-6, 2LPV=0.4e-6, old = 0.3e-6
pvs = 7.0e-6 # stratospheric PV in (m^2*K)/(s*kg), CTSS=7.0e-6, 2LPV=4.0e-6, old = 5.0e-6
pim = 665. # mean height of tropopause in pi, CTSS=680, 2LPV=720, old = 660
dpitr = 30. # max displacement from mean tropopause height in pi, CTSS=30, 2LPV=50
dpipv = 15. # tropopause depth in pi over which strongest PV change occurs, CTSS=15, 2LPV=15
thtop = 517. # mean top potential temperature in K, CTSS=535, 2LPV=420, old = 526
dth = 10. # max displacement from mean top potential temperature in K, CTSS=10, 2LPV=10
dyth = 1.0e6 # meridional scale for top potential temperature transition in m, CTSS=1.0e6, 2LPV=1.0e6
dytr = 5.0e5  # meridional scale for tropopause height transition in m, CTSS=5.0e5, 2LPV=5.0e5
shear_type = 'none' # shear option, none or cyc or anticyc
dyph = 1.0e6 # meridional scale for bottom phi
phbot = 0. # mean phi at bottom (pressure level p_bot)
dph = 750. # max displacement from mean bottom phi

### Set up grids
ny_h = ny
ly = ny_h*hres*1000
dy_l, ylvls_l, dy_h, ylvls_h = epv_jet.y_grid(ly, ny_l, ny_h)
dpi_l, pilvls_l, dpi_h, pilvls_h, p_pi_l, p_pi_h = epv_jet.pi_grid(p_bot, pibot, pitop, npi_l, npi_h)
znw, znu, dnw, rdnw, dn, rdn, fnp, fnm, cof1, cof2, cf1,\
    cf2, cf3, cfn, cfn1, dx, dy, rdx, rdy = epv_jet.eta_grid(nz, hres)

### Low-resolution run
f_l = np.zeros((npi_l,ny_l))
pv_l, f_l, u_l, theta_l, pv_out_l = epv_jet.solve_PV_inversion(ly, ny_l, \
    dy_l, dytr, pim, dpitr, pvt, pvs, dpipv, npi_l, pibot, pitop, dpi_l, \
    dyth, thtop, dth, f_l, om, nit, pv_dist, dyph, phbot, dph, shear_type)

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
    dyth, thtop, dth, f_interp, om, nit, pv_dist, dyph, phbot, dph, shear_type)

### Unstaggered pressure
p_unstag = znu*(p_bot - p_top) + p_top

### Interpolate onto eta coordinates
u_eta = np.zeros((nz,ny))
theta_eta = np.zeros((nz,ny))
ph_eta = np.zeros((nz,ny))
p_eta = np.zeros((nz,ny))
for j in range(ny):
    u_eta[:,j] = epv_jet.interp_decreasing(p_pi_h,p_unstag,u_h[:,j],npi_h,nz)
    theta_eta[:,j] = epv_jet.interp_decreasing(p_pi_h,p_unstag,theta_h[:,j],npi_h,nz)
    ph_eta[:,j] = epv_jet.interp_decreasing(p_pi_h,p_unstag,f_h[:,j],npi_h,nz)
    p_eta[:,j] = p_unstag
    
### Bottom and top pressure for shear cases
p_bot_inc = np.zeros(ny)
p_top_inc = np.zeros(ny)
if shear_type == 'anticyc' or shear_type == 'cyc':
    rho_eta = P0/(RD*theta_eta)*(p_eta/P0)**(CV/CP)
    u_bot = 1.5*u_eta[0,:] - 0.5*u_eta[1,:]
    u_top = 1.5*u_eta[-1,:] - 0.5*u_eta[-2,:]
    rho_bot = 1.5*rho_eta[0,:] - 0.5*rho_eta[1,:]
    rho_top = 1.5*rho_eta[-1,:] - 0.5*rho_eta[-2,:]
    for j in range(1,ny-1):
        p_bot_inc[j] = p_bot_inc[j-1] - hres*1000.*u_bot[j]*rho_bot[j]*F0
        p_top_inc[j] = p_top_inc[j-1] - hres*1000.*u_top[j]*rho_top[j]*F0

    p_bot_inc[0] = p_bot_inc[1]
    p_bot_inc[-1] = p_bot_inc[-2]
    p_top_inc[0] = p_top_inc[1]
    p_top_inc[-1] = p_top_inc[-2]
    
### Extend into x
u_jet = np.zeros((nz,ny,nx))
v_jet = np.zeros((nz,ny,nx))
theta_jet = np.zeros((nz,ny,nx))
ph_jet = np.zeros((nz,ny,nx))
p_jet = np.zeros((nz,ny,nx))
pbot = np.zeros((ny,nx))
ptop = np.zeros((ny,nx))

for i in range(nx):
    u_jet[:,:,i] = u_eta
    theta_jet[:,:,i] = theta_eta
    ph_jet[:,:,i] = ph_eta
    p_jet[:,:,i] = p_eta
    pbot[:,i] = p_bot_inc + p_bot
    ptop[:,i] = p_top_inc + p_top

### Compute WRF base variables
qv = np.zeros((nz,ny,nx))
mub, pb, t_init, alb, phb = epv_jet.base_wrf(ph_jet, p_jet, theta_jet, pbot, ptop,\
                                             nx, ny, nz, znu, znw, dn, dnw)

### Compute WRF perturbation variables
u, v, ph, p, mu, t, tsk = epv_jet.pert_wrf(u_jet, v_jet, ph_jet, p_jet, theta_jet,\
                                           pb, alb, mub, pbot, ptop, nx, ny, nz, \
                                           znu, znw, dn, dnw)

hres_m = hres*1000.
ncfile = nc.Dataset(file_name,'w',format=netcdf_type)
ncfile = write_wrfinputd01.write(ncfile, nx, ny, nz, hres_m, title_str, \
        time_str, u, v, t, ph, phb, t_init, mu, mub, p, pb, \
        fnm, fnp, rdnw, rdn, dnw, dn, cfn, cfn1, rdx, rdy, cf1, \
        cf2, cf3, qv, znw, znu, p_top, tsk)
ncfile.close()
