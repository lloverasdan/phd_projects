#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Daniel Lloveras
"""

#%%%%%%%%%% IMPORT MODULES AND DEFINE PATHS

### Modules
import numpy as np
import netCDF4 as nc
from bwave_ideal_wrf import epv_jet, qgpv_pert, wrf_fields, write_wrfinputd01

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
file_name = ''
netcdf_type = 'NETCDF3_64BIT_OFFSET'
title_str = "OUTPUT FROM IDEAL V3.6.1 PREPROCESSOR"
time_str = "2021-01-01_00:00:00"

#%%%%%%%%%% USER OPTIONS

### Grid parameters
nx = 250 # number of grid points in x direction
ny = 225 # number of grid points in y direction
nz = 100 # number of grid points in z direction
hres = 32. # horizontal grid resolution in km
zl = 20. # model top in km
pbot = 101200. # bottom pressure in Pa
ptop = 5000. # top pressure in Pa
pibot = 1004. # bottom height for EPV inversion in pi
pitop = 424. # top height for EPV inversion in pi
npi_h = 180 # number of grid points in pi for high-resolution EPV inversion
npi_l = 40 # number of grid points in pi for low-resolution EPV inversion
ny_l = 30 # number of grid points in y for low-resolution EPV inversion
nit = 50000  # number of iterations for EPV inversion
om = 1.8  # successive over-relaxation coefficient for EPV inversion

### Background jet parameters
pv_dist = 1 # PV distribution option, CTSS=1, 2LPV=2
pvt = 0.3e-6 # tropospheric PV in (m^2*K)/(s*kg), CTSS=0.3e-6, 2LPV=0.4e-6
pvs = 5.0e-6 # stratospheric PV in (m^2*K)/(s*kg), CTSS=5.0e-6, 2LPV=4.0e-6
pim = 660. # mean height of tropopause in pi, CTSS=660, 2LPV=720
dpitr = 30. # max displacement from mean tropopause height in pi, CTSS=30, 2LPV=50
dpipv = 15. # tropopause depth in pi over which strongest PV change occurs, CTSS=15, 2LPV=15
thtop = 526. # mean top potential temperature in K, CTSS=526, 2LPV=417
dth = 10. # max displacement from mean top potential temperature in K, CTSS=10, 2LPV=0
dyth = 1.0e6 # meridional scale for top potential temperature transition in m, CTSS=1.0e6, 2LPV=1.0e6
dytr = 5.0e5  # meridional scale for tropopause height transition in m, CTSS=5.0e5, 2LPV=5.0e5
avg_rh = 0.72 # average relative humidity in the domain for moist soundings CTSS=0.72, 2LPV=0.75

### Upper-level QGPV-anomaly parameters
up_pert = True # upper-level nomaly option, Yes=True, No=False
up_pert_type = 2 # upper-level anomaly structure, Gaussian=1, Cosine=2, Theta=3
up_pert_mag = 1.0e-4 # magnitude of upper-level anomaly in s^-1
x_up_pert = 94 # x center gridpoint for upper-level anomaly
y_up_pert = 112 # y center gridpoint for upper-level anomaly
z_up_pert = 50  # z center gridpoint for upper-level anomaly, only if up_pert_type=1,2
th_up_pert = 310 # theta center value in K for upper-level anomaly, only if up_pert_type=3
ax_up_pert = 200. # x decay scale of upper-level anomaly in km
ay_up_pert = 600. # y decay scale of upper-level anomaly in km
az_up_pert = 1.5 # z decay scale of upper-level anomaly in km, only if up_pert_type=1,2
ath_up_pert = 2.0 # theta decay scale of upper-level anomaly in K, only if up_pert_type=3

### Surface potential-temperature anomaly parameters
surf_pert = True # surface anomaly option, Yes=True, No=False
surf_pert_mag = 5.0 # magnitude of surface theta anomaly in K
x_surf_pert = 130 # x center gridpoint for surface anomaly
y_surf_pert = 85 # y center gridpoint for surface anomaly
ax_surf_pert = 600. # x decay scale for surface anomaly in km
ay_surf_pert = 200. # y decay scale for surface anomaly in km

#%%%%%%%%%% INITIAL JET

### Set up grids
ny_h = ny
ly = ny_h*hres*1000
dy_l, ylvls_l, dy_h, ylvls_h = epv_jet.y_grid(ly, ny_l, ny_h)
dpi_l, pilvls_l, dpi_h, pilvls_h, p_pi_l, p_pi_h = epv_jet.pi_grid(pbot, pibot, pitop, npi_l, npi_h)
znw, znu = epv_jet.eta_grid(nz)

### Low-resolution run
f_l = np.zeros((npi_l,ny_l))
pv_l, f_l, u_l, theta_l, pv_out_l = epv_jet.solve_PV_inversion(ly, ny_l, \
    dy_l, dytr, pim, dpitr, pvt, pvs, dpipv, npi_l, pibot, pitop, dpi_l, \
    dyth, thtop, dth, f_l, om, nit, pv_dist)

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
    dyth, thtop, dth, f_interp, om, nit, pv_dist)

### Compute on eta grid
p_jet, theta_jet, u_jet, rho_jet, z_jet, dtdz_jet, n_jet = \
    epv_jet.eta_calc(ny_h, nz, pbot, ptop, znu, npi_h, p_pi_h, u_h, theta_h, f_h)

### Surface pressure
p_surf_jet = np.ones(ny)
p_surf_jet = p_surf_jet*pbot

#%%%%%%%%%% QGPV ANOMALIES

### Set up grids
xl, yl, x, y, xg, yg, dz, zp, facz = qgpv_pert.cartesian_mesh(nx, ny, nz, hres, zl)
kmax, lmax, facx, facy, dxg, dyg = qgpv_pert.spectral_mesh(nx, ny, xl, yl)

n_up_pert = n_jet[:,y_up_pert]
n_surf_pert = n_jet[:,y_surf_pert]
dtdz_up_pert = dtdz_jet[:,y_up_pert]
dtdz_surf_pert = dtdz_jet[:,y_surf_pert]


### Upper-level anomaly
if up_pert:
    
    ### Initialize anomaly
    if up_pert_type==1:
        pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac = qgpv_pert.qgpv_anom_gauss(up_pert_mag,\
            x_up_pert, y_up_pert, z_up_pert, az_up_pert, ax_up_pert, ay_up_pert, x, y, zp, xg, yg, nx, ny, \
            nz, n_up_pert, facz, dz, dtdz_up_pert)
    elif up_pert_type==2:
        pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac = qgpv_pert.qgpv_anom_cos(up_pert_mag,\
            x_up_pert, y_up_pert, z_up_pert, az_up_pert, ax_up_pert, ay_up_pert, x, y, zp, xg, yg, nx, ny, \
            nz, n_up_pert, facz, dz, dtdz_up_pert)
    elif up_pert_type==3:
        pvxy, ubcxy, lbcxy, pvsp, ubcsp, lbcsp, bu_fac = qgpv_pert.qgpv_theta(qgpv_mag,\
            x_up_pert, y_up_pert, th_up_pert, ath_up_pert, ax_up_pert, ay_up_pert, x, y, zp, xg, yg, nx, ny, \
            nz, n_up_pert, facz, dz, theta_jet, z_jet, dtdz_up_pert)
    
    ### Invert anomaly
    fbsp, ftsp, fzbsp, fztsp, fsp = qgpv_pert.qgpv_inversion(nx, ny, nz, bu_fac, \
            facx, facy, facz, kmax, lmax, pvsp, ubcsp, lbcsp, dz)

    ### Compute u, v, rho, theta, p perturbation fields
    u_up_pert, v_up_pert, theta_up_pert, rho_up_pert, fxy_up = qgpv_pert.qgpv_solver(fsp, nx, ny, nz, \
            bu_fac, dxg, dyg, dz, lbcxy, ubcxy, dtdz_up_pert)

    p_up_pert = fxy_up*F0
    
elif up_pert == False:
    u_up_pert = np.zeros((nz,ny,nx))
    v_up_pert = np.zeros((nz,ny,nx))
    theta_up_pert = np.zeros((nz,ny,nx))
    rho_up_pert = np.zeros((nz,ny,nx))
    p_up_pert = np.zeros((nz,ny,nx))
    
if surf_pert:
    
    ### Initialize anomaly
    pvs, ubcs, lbcs, pvsps, ubcsps, lbcsps, bu_facs = qgpv_pert.qgpv_surf(surf_pert_mag,\
        x_surf_pert, y_surf_pert, ax_surf_pert, ay_surf_pert, x, y, zp, xg, yg, nx, ny, nz,\
        n_surf_pert, facz, dz, dtdz_surf_pert)
    
    ### Invert anomaly
    fbsps, ftsps, fzbsps, fztsps, fsps = qgpv_pert.qgpv_inversion(nx, ny, nz, bu_facs, \
        facx, facy, facz, kmax, lmax, pvsps, ubcsps, lbcsps, dz)
    
    ### Compute u, v, rho, theta, p perturbation fields
    u_surf_pert, v_surf_pert, theta_surf_pert, rho_surf_pert, fxy_surf = qgpv_pert.qgpv_solver(fsps, nx, ny, nz, \
        bu_facs, dxg, dyg, dz, lbcs, ubcs, dtdz_surf_pert)

    p_surf_pert = fxy_surf*F0

elif surf_pert == False:
    u_surf_pert = np.zeros((nz,ny,nx))
    v_surf_pert = np.zeros((nz,ny,nx))
    theta_surf_pert = np.zeros((nz,ny,nx))
    rho_surf_pert = np.zeros((nz,ny,nx))
    p_surf_pert = np.zeros((nz,ny,nx))

### Add perturbations
u_pert = u_up_pert + u_surf_pert
v_pert = v_up_pert + v_surf_pert
theta_pert = theta_up_pert + theta_surf_pert
rho_pert = rho_up_pert + rho_surf_pert
p_pert = p_up_pert + p_surf_pert

#%%%%%%%%%% MOISTURE AND WRF GRID

### Set up 
dnw, rdnw, dn, rdn, fnp, fnm, cof1, cof2, cf1, cf2, cf3, cfn, cfn1, rdx,\
     rdy = wrf_fields.wrf_grid(nx, ny, nz, hres, znw)

### Add QGPV perturbation
u_eta, v_eta, theta_eta, rho_eta, pres_eta, z_eta, up_eta, vp_eta, \
    thetap_eta, rhop_eta, presp_eta = wrf_fields.add_pert(u_jet, theta_jet, \
    rho_jet, z_jet, p_jet, zp, u_pert, v_pert, rho_pert, theta_pert, p_pert, nx, ny, nz)

p_surf_eta = np.zeros((ny,nx))
for i in range(nx):
    p_surf_eta[:,i] = p_surf_jet

p_bot = p_surf_eta + p_pert[0,:,:]
p_top = np.ones((ny,nx))*ptop + p_pert[-1,:,:]

### Compute base fields in middle of domain
mub, pb, t_init, alb, phb = wrf_fields.middle_sounding(nx, ny, nz, p_bot, \
    p_top, znu, dnw, z_eta, theta_eta, rho_eta, u_eta, v_eta, pres_eta, avg_rh)             
        
### Compute fields in entire domain
u, v, t, ph, mu, p, moist, tsk = wrf_fields.full_domain(nx, ny, nz, p_bot, p_top, znu, \
    dnw, mub, pb, alb, dn, z_eta, theta_eta, rho_eta, u_eta, v_eta, pres_eta, avg_rh)

#%%%%%%%%%% WRITE NETCDF FILE

hres_m = hres*1000.
ncfile = nc.Dataset(file_name,'w',format=netcdf_type)
ncfile = write_wrfinputd01.write(ncfile, nx, ny, nz, hres_m, title_str, \
        time_str, u, v, t, ph, phb, t_init, mu, mub, p, pb, \
        fnm, fnp, rdnw, rdn, dnw, dn, cfn, cfn1, rdx, rdy, cf1, \
        cf2, cf3, moist, znw, znu, ptop, tsk)
ncfile.close()
