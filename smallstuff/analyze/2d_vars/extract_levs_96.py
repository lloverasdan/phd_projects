#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

from netCDF4 import Dataset
import numpy as np
from wrf import getvar, interplevel

### Times
ti = np.asarray([12,30,48])
nx = 2000
ny = 1800

### Load the netCDF files
nc_ctl_96 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_96h/wrfout_d01_2021-01-05_00_00_00')
nc_adj_96 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/adj_96h/wrfout_d01_2021-01-05_00_00_00')
nc_adj10_96 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/adj_96h_10/wrfout_d01_2021-01-05_00_00_00')
nc_adj100_96 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/adj_96h_100/wrfout_d01_2021-01-05_00_00_00')
nc_wave_96 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/wave_96h/wrfout_d01_2021-01-05_00_00_00')
nc_wave10_96 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/wave_96h_10/wrfout_d01_2021-01-05_00_00_00')
nc_wave100_96 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/wave_96h_100/wrfout_d01_2021-01-05_00_00_00')

### Read in the data
plev = 325.
ulev_ctl = np.zeros((len(ti),ny,nx))
# vlev_ctl = np.zeros((len(ti),ny,nx))
# zlev_ctl = np.zeros((len(ti),ny,nx))

ulev_adj = np.zeros((len(ti),ny,nx))
# vlev_adj = np.zeros((len(ti),ny,nx))
# zlev_adj = np.zeros((len(ti),ny,nx))

ulev_adj10 = np.zeros((len(ti),ny,nx))
# vlev_adj10 = np.zeros((len(ti),ny,nx))
# zlev_adj10 = np.zeros((len(ti),ny,nx))

ulev_adj100 = np.zeros((len(ti),ny,nx))
# vlev_adj100 = np.zeros((len(ti),ny,nx))
# zlev_adj100 = np.zeros((len(ti),ny,nx))

ulev_wave = np.zeros((len(ti),ny,nx))
# vlev_wave = np.zeros((len(ti),ny,nx))
# zlev_wave = np.zeros((len(ti),ny,nx))

ulev_wave10 = np.zeros((len(ti),ny,nx))
# vlev_wave10 = np.zeros((len(ti),ny,nx))
# zlev_wave10 = np.zeros((len(ti),ny,nx))

ulev_wave100 = np.zeros((len(ti),ny,nx))
# vlev_wave100 = np.zeros((len(ti),ny,nx))
# zlev_wave100 = np.zeros((len(ti),ny,nx))

for i in range(len(ti)):
    p_ctl = np.asarray(getvar(nc_ctl_96,'th',timeidx=int(ti[i]/3)))
    p_adj = np.asarray(getvar(nc_adj_96,'th',timeidx=int(ti[i]/3)))
    p_adj10 = np.asarray(getvar(nc_adj10_96,'th',timeidx=int(ti[i]/3)))
    p_adj100 = np.asarray(getvar(nc_adj100_96,'th',timeidx=int(ti[i]/3)))
    p_wave = np.asarray(getvar(nc_wave_96,'th',timeidx=int(ti[i]/3)))
    p_wave10 = np.asarray(getvar(nc_wave10_96,'th',timeidx=int(ti[i]/3)))
    p_wave100 = np.asarray(getvar(nc_wave100_96,'th',timeidx=int(ti[i]/3)))

    u_ctl = np.asarray(getvar(nc_ctl_96,'pvo',timeidx=int(ti[i]/3)))
    u_adj = np.asarray(getvar(nc_adj_96,'pvo',timeidx=int(ti[i]/3)))
    u_adj10 = np.asarray(getvar(nc_adj10_96,'pvo',timeidx=int(ti[i]/3)))
    u_adj100 = np.asarray(getvar(nc_adj100_96,'pvo',timeidx=int(ti[i]/3)))
    u_wave = np.asarray(getvar(nc_wave_96,'pvo',timeidx=int(ti[i]/3)))
    u_wave10 = np.asarray(getvar(nc_wave10_96,'pvo',timeidx=int(ti[i]/3)))
    u_wave100 = np.asarray(getvar(nc_wave100_96,'pvo',timeidx=int(ti[i]/3)))

#     v_ctl = np.asarray(getvar(nc_ctl_96,'va',timeidx=int(ti[i]/3)))
#     v_adj = np.asarray(getvar(nc_adj_96,'va',timeidx=int(ti[i]/3)))
#     v_adj10 = np.asarray(getvar(nc_adj10_96,'va',timeidx=int(ti[i]/3)))
#     v_adj100 = np.asarray(getvar(nc_adj100_96,'va',timeidx=int(ti[i]/3)))
#     v_wave = np.asarray(getvar(nc_wave_96,'va',timeidx=int(ti[i]/3)))
#     v_wave10 = np.asarray(getvar(nc_wave10_96,'va',timeidx=int(ti[i]/3)))
#     v_wave100 = np.asarray(getvar(nc_wave100_96,'va',timeidx=int(ti[i]/3)))

#     z_ctl = np.asarray(getvar(nc_ctl_96,'z',timeidx=int(ti[i]/3),units='dm'))
#     z_adj = np.asarray(getvar(nc_adj_96,'z',timeidx=int(ti[i]/3),units='dm'))
#     z_adj10 = np.asarray(getvar(nc_adj10_96,'z',timeidx=int(ti[i]/3),units='dm'))
#     z_adj100 = np.asarray(getvar(nc_adj100_96,'z',timeidx=int(ti[i]/3),units='dm'))
#     z_wave = np.asarray(getvar(nc_wave_96,'z',timeidx=int(ti[i]/3),units='dm'))
#     z_wave10 = np.asarray(getvar(nc_wave10_96,'z',timeidx=int(ti[i]/3),units='dm'))
#     z_wave100 = np.asarray(getvar(nc_wave100_96,'z',timeidx=int(ti[i]/3),units='dm'))

    ulev_ctl[i,:,:] = interplevel(u_ctl,p_ctl,plev)
    ulev_adj[i,:,:] = interplevel(u_adj,p_adj,plev)
    ulev_adj10[i,:,:] = interplevel(u_adj10,p_adj10,plev)
    ulev_adj100[i,:,:] = interplevel(u_adj100,p_adj100,plev)
    ulev_wave[i,:,:] = interplevel(u_wave,p_wave,plev)
    ulev_wave10[i,:,:] = interplevel(u_wave10,p_wave10,plev)
    ulev_wave100[i,:,:] = interplevel(u_wave100,p_wave100,plev)

#     vlev_ctl[i,:,:] = interplevel(v_ctl,p_ctl,plev)
#     vlev_adj[i,:,:] = interplevel(v_adj,p_adj,plev)
#     vlev_adj10[i,:,:] = interplevel(v_adj10,p_adj10,plev)
#     vlev_adj100[i,:,:] = interplevel(v_adj100,p_adj100,plev)
#     vlev_wave[i,:,:] = interplevel(v_wave,p_wave,plev)
#     vlev_wave10[i,:,:] = interplevel(v_wave10,p_wave10,plev)
#     vlev_wave100[i,:,:] = interplevel(v_wave100,p_wave100,plev)

#     zlev_ctl[i,:,:] = interplevel(z_ctl,p_ctl,plev)
#     zlev_adj[i,:,:] = interplevel(z_adj,p_adj,plev)
#     zlev_adj10[i,:,:] = interplevel(z_adj10,p_adj10,plev)
#     zlev_adj100[i,:,:] = interplevel(z_adj100,p_adj100,plev)
#     zlev_wave[i,:,:] = interplevel(z_wave,p_wave,plev)
#     zlev_wave10[i,:,:] = interplevel(z_wave10,p_wave10,plev)
#     zlev_wave100[i,:,:] = interplevel(z_wave100,p_wave100,plev)
        
### Save the output
np.save('/p/work1/lloveras/adj_4km/processed/96h/325k_pv/pv325_ctl',ulev_ctl)
np.save('/p/work1/lloveras/adj_4km/processed/96h/325k_pv/pv325_adj',ulev_adj)
np.save('/p/work1/lloveras/adj_4km/processed/96h/325k_pv/pv325_adj10',ulev_adj10)
np.save('/p/work1/lloveras/adj_4km/processed/96h/325k_pv/pv325_adj100',ulev_adj100)
np.save('/p/work1/lloveras/adj_4km/processed/96h/325k_pv/pv325_wave',ulev_wave)
np.save('/p/work1/lloveras/adj_4km/processed/96h/325k_pv/pv325_wave10',ulev_wave10)
np.save('/p/work1/lloveras/adj_4km/processed/96h/325k_pv/pv325_wave100',ulev_wave100)

# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/v900_ctl',vlev_ctl)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/v900_adj',vlev_adj)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/v900_adj10',vlev_adj10)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/v900_adj100',vlev_adj100)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/v900_wave',vlev_wave)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/v900_wave10',vlev_wave10)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/v900_wave100',vlev_wave100)

# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/z900_ctl',zlev_ctl)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/z900_adj',zlev_adj)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/z900_adj10',zlev_adj10)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/z900_adj100',zlev_adj100)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/z900_wave',zlev_wave)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/z900_wave10',zlev_wave10)
# np.save('/p/work1/lloveras/adj_4km/processed/96h/900/z900_wave100',zlev_wave100)
