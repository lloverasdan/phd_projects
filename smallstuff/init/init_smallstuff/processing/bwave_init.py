#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 09:18:04 2020

@author: lloverasdan
"""

#%%%%%%%%%% IMPORT MODULES

import numpy as np
import netCDF4 as nc
import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from bwave_ideal_wrf import epv_jet, qgpv_pert, wrf_fields, write_wrfinputd01

#%%%%%%%%%% USER OPTIONS AND PARAMETERS

### Background state options
psurf_gradient = False
shear = False
baroclinic = False

### PV perturbation options
nopert = False
pert_type = 3 # 1: Gaussian, 2: Cosine, 3: Theta
surf_pert = True
extra_surfpert = False

### Grid parameters
nx = 2000  # number of grid points in x direction
ny = 1800  # number of grid points in y direction
nz = 100  # number of grid points in z direction
hres = 4.  # horizontal grid resolution, in km
zl = 20.  # height of top of atmosphere, in km
pbot = 102000.  # bottom pressure, in Pa
ptop = 5000.  # top pressure, in Pa

### Anomaly grid parameters
x_pert = 600 # upper
y_pert = 850 # upper
th_pert = 310 # upper
z_pert = 40  # upper

xs_pert = 750 # lower
ys_pert = 650 # lower

xsp_pert = 625
ysp_pert = 700

### Jet parameters
pibot = 1004.  # bottom height, in pi
pitop = 424.  # top height, in pi
npi_h = 180  # number of grid points in pi
npi_l = 40  # number of grid points in pi for low-resolution run
ny_l = 30  # number of grid points in y for low-resolution EPV inversion
thtop = 525  # top potential temperature, in K
dth = 10.  # max displacement from thtop
dyth = 1e6  # meridional scale for theta transition, in m
pim = 660.  # mean height of tropopause, in pi coordinates
dpitr = 32.  # max displacement from pim
dytr = 5e5  # meridional scale for tropopause transition, in m
pvs = 5.0e-6  # PV at bottom of stratosphere, in (m^2*K)/(s*kg)
pvt = 0.3e-6  # PV at bottom of troposhere
dpipv = 15.  # tropopause depth in pi over which strongest PV change occurs
nit = 50000  # number of iterations
om = 1.8  # successive over-relaxation coefficient
avg_rh = 0.72 # average relative humidity in the domain for moist soundings
c = 1.0  # additional jet speed if pressure gradient option is enabled
c_shear_linear = 2.3e-7 # s^{-1}, positive = anticyclonic, negative = cyclonic
c_shear_sin = 2.0 # ms^{-1}, positive = anticyclonic, negative = cyclonic
p_mid = pbot

### Perturbation parameters
qgpv_mag = 3.0e-4  # magnitude of QGPV anomaly, in s-1
theta_mag = 0.  # magnitude of boundary theta, in K
az = 3.0  # vertical decay scale of anomaly, in km
ax = 200.  # x decay scale of anomaly, in km
ay = 600.  # y decay scale of anomaly, in km
ath = 2.0  # theta decay scale of anomaly, in K
dtheta_surf = 5.
dqgpv_surf = 0.
axs = 600.
ayn = 200.
ays = 200.
azs = 1.0
dtheta_pert = 1.
dqgpv_pert = 0.
axsp = 200.
aynp = 200.
aysp = 200.
azsp = 1.0

### Input file
file_name = '/p/work/lloveras/new/in_files/wrfin_ctl'
netcdf_type = 'NETCDF3_64BIT_OFFSET'
title_str = "OUTPUT FROM IDEAL V3.6.1 PREPROCESSOR"  # title of netCDF file
time_str = "2021-01-01_00:00:00"  # start time of simulation

### Namelist variables that appear in wrfinput_d01 file
dt = 20.
mp_physics = 17
ra_lw_physics = 0
ra_sw_physics = 0
sf_sfclay_physics = 1
sf_surface_physics = 0
bl_pbl_physics = 1
cu_physics = 0
num_soil_layers = 5
diff_opt = 1
km_opt = 4
damp_opt = 3
dampcoef = 0.4
khdif = 0.
kvdif = 0.
spec_bdy_width = 8
sf_lake_physics = 0
surface_input_source = 1
hypsometric_opt = 0
num_land_cat = 24
num_soil_cat = 16

#%%%%%%%%%% INITIAL JET

### Set up grids
ny_h = ny
ly = ny_h*hres*1000
dy_l, ylvls_l, dy_h, ylvls_h = epv_jet.y_grid(ly, ny_l, ny_h)
dpi_l, pilvls_l, dpi_h, pilvls_h, p_pi_l, p_pi_h = epv_jet.pi_grid(pbot, \
                                                pibot, pitop, npi_l, npi_h)
znw, znu = epv_jet.eta_grid(nz)

### Low-resolution run
f_l = np.zeros((npi_l,ny_l))
pv_l, f_l, u_l, theta_l, pv_out_l = epv_jet.solve_PV_inversion(ly, ny_l, \
    dy_l, dytr, pim, dpitr, pvt, pvs, dpipv, npi_l, pibot, pitop, dpi_l, \
    dyth, thtop, dth, f_l, om, nit)

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
    dyth, thtop, dth, f_interp, om, nit)

### Compute on eta grid
p_jet, theta_jet, u_jet, rho_jet, z_jet, dtdz_jet, n_jet = \
    epv_jet.eta_calc(ny_h, nz, pbot, ptop, znu, npi_h, p_pi_h, u_h, theta_h)

### Surface pressure
p_surf_jet = np.ones(ny)
p_surf_jet = p_surf_jet*pbot
if psurf_gradient:
    u_jet, rho_jet, p_jet, p_surf_jet, z_jet = epv_jet.faster_jet(nz, ny, u_jet,\
                            theta_jet, p_jet, p_mid, c, hres, ptop, z_jet,znu)

if shear:
    u_jet, rho_jet, p_jet, p_surf_jet, z_jet = epv_jet.baro_shear_sin(nz, ny, u_jet,\
                            theta_jet, p_jet, p_mid, c_shear_sin, hres, ptop, z_jet,znu)
    
if baroclinic:
    u_jet, rho_jet, p_jet, p_surf_jet, z_jet = epv_jet.baroclinic_shear(nz, ny, u_jet,\
                            theta_jet, p_jet, p_mid, c_shear_sin, hres, ptop, z_jet,znu)

#%%%%%%%%%% QGPV ANOMALY

### Set up grids
xl, yl, x, y, xg, yg, dz, z, facz = qgpv_pert.cartesian_mesh(nx, ny, nz, hres, zl)
kmax, lmax, facx, facy, dxg, dyg = qgpv_pert.spectral_mesh(nx, ny, xl, yl)

### Initialize anomaly
n_pert = n_jet[:,y_pert]
ns_pert = n_jet[:,ys_pert]
dtdz_pert = dtdz_jet[:,y_pert]
dtdzs_pert = dtdz_jet[:,ys_pert]
if pert_type==1:
    pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac = qgpv_pert.qgpv_anom_gauss(qgpv_mag,\
        theta_mag, x_pert, y_pert, z_pert, az, ax, ay, x, y, z, xg, yg, nx, ny, \
        nz, n_pert, facz, dz, dtdz_pert)
elif pert_type==2:
    pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac = qgpv_pert.qgpv_anom_cos(qgpv_mag,\
        theta_mag, x_pert, y_pert, z_pert, az, ax, ay, x, y, z, xg, yg, nx, ny, \
        nz, n_pert, facz, dz, dtdz_pert)    
elif pert_type==3:
    pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac = qgpv_pert.qgpv_theta(qgpv_mag,\
        theta_mag, x_pert, y_pert, th_pert, ath, ax, ay, x, y, z, xg, yg, nx, ny, \
        nz, n_pert, facz, dz, theta_jet, z_jet, dtdz_pert)
    
### Invert anomaly
fbsp, ftsp, fzbsp, fztsp, fsp = qgpv_pert.qgpv_inversion(nx, ny, nz, bu_fac, \
        facx, facy, facz, kmax, lmax, pvsp, ubcsp, lbcsp, dz)

### Compute u, v, rho, theta perturbation fields
u_pert, v_pert, theta_pert, rho_pert, fxy = qgpv_pert.qgpv_solver(fsp, nx, ny, nz, \
        bu_fac, dxg, dyg, dz, lbcxy, ubcxy, dtdz_pert)

if nopert:
    u_pert = np.zeros(np.shape(fxy))
    v_pert = np.zeros(np.shape(fxy))
    theta_pert = np.zeros(np.shape(fxy))
    rho_pert = np.zeros(np.shape(fxy))

if surf_pert:
    pvs, ubcs, lbcs, pvsps, ubcsps, lbcsps, bu_facs = qgpv_pert.qgpv_surf(dqgpv_surf, dtheta_surf,\
        xs_pert, ys_pert, azs, axs, ayn, ays, x, y, z, xg, yg, nx, ny, nz, ns_pert, facz, dz, dtdzs_pert)
    fbsps, ftsps, fzbsps, fztsps, fsps = qgpv_pert.qgpv_inversion(nx, ny, nz, bu_facs, \
        facx, facy, facz, kmax, lmax, pvsps, ubcsps, lbcsps, dz)
    us_pert, vs_pert, thetas_pert, rhos_pert, fxys = qgpv_pert.qgpv_solver(fsps, nx, ny, nz, \
        bu_facs, dxg, dyg, dz, lbcs, ubcs, dtdzs_pert)

    u_pert = u_pert + us_pert
    v_pert = v_pert + vs_pert
    theta_pert = theta_pert + thetas_pert
    rho_pert = rho_pert + rhos_pert
    
if extra_surfpert:
    pvsp, ubcsp, lbcsp, pvspsp, ubcspsp, lbcspsp, bu_facsp = qgpv_pert.qgpv_surf(dqgpv_pert, dtheta_pert,\
        xsp_pert, ysp_pert, azsp, axsp, aynp, aysp, x, y, z, xg, yg, nx, ny, nz, nsp_pert, facz, dz, dtdzsp_pert)
    fbspsp, ftspsp, fzbspsp, fztspsp, fspsp = qgpv_pert.qgpv_inversion(nx, ny, nz, bu_facsp, \
        facx, facy, facz, kmax, lmax, pvspsp, ubcspsp, lbcspsp, dz)
    usp_pert, vsp_pert, thetasp_pert, rhosp_pert, fxysp = qgpv_pert.qgpv_solver(fspsp, nx, ny, nz, \
        bu_facsp, dxg, dyg, dz, lbcsp, ubcsp, dtdzsp_pert)

    u_pert = u_pert + usp_pert
    v_pert = v_pert + vsp_pert
    theta_pert = theta_pert + thetasp_pert
    rho_pert = rho_pert + rhosp_pert      

#%%%%%%%%%% MOISTURE AND WRF GRID

### Set up 
dnw, rdnw, dn, rdn, fnp, fnm, cof1, cof2, cf1, cf2, cf3, cfn, cfn1, rdx,\
     rdy = wrf_fields.wrf_grid(nx, ny, nz, hres, znw)

### Add QGPV perturbation to initial jet    
u, v, theta, rho, z = wrf_fields.add_pert(u_jet, theta_jet, rho_jet, \
    z_jet, u_pert, v_pert, rho_pert, theta_pert, nx, ny, nz)

### Compute base fields in middle of domain
mub, pb, t_init, alb, phb = wrf_fields.middle_sounding(nx, ny, nz, p_surf_jet, \
    ptop, znu, dnw, z, theta, rho, u, v, avg_rh)

### Compute fields in entire domain
u_field, v_field, t, ph, mu, p, moist, t_base, qv_base, u_base, v_base, \
    tsk, tmn = wrf_fields.full_domain(nx, ny, nz, p_surf_jet, ptop, znu, \
                dnw, mub, pb, alb, dn, z, theta, rho, u, v, avg_rh)

#%%%%%%%%%% WRITE NETCDF FILE

hres_m = hres*1000.
ncfile = nc.Dataset(file_name,'w',format=netcdf_type)
ncfile = write_wrfinputd01.write(ncfile, nx, ny, nz, hres_m, title_str, \
        time_str, u_field, v_field, t, ph, phb, t_init, mu, mub, p, pb, \
        fnm, fnp, rdnw, rdn, dnw, dn, t_base, cfn, cfn1, rdx, rdy, cf1, \
        cf2, cf3, tsk, u_base, v_base, qv_base, tmn, moist, znw, znu, \
        diff_opt, km_opt, damp_opt, dampcoef, khdif, kvdif, mp_physics, \
        ra_lw_physics, ra_sw_physics, sf_sfclay_physics, sf_surface_physics, \
        bl_pbl_physics, cu_physics, sf_lake_physics, surface_input_source, \
        hypsometric_opt, dt, num_land_cat, num_soil_layers, num_soil_cat, \
        spec_bdy_width, ptop)
ncfile.close()
