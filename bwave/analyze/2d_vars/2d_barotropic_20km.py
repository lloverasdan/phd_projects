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
ti = np.asarray([0,1,2,3,4,5,6,7,8]) # time indices in days
nx = 400
ny = 360

### Load the netCDF files
nc1 = Dataset('/p/work1/lloveras/bwave_nov/20km_files/barotropic/output/wrfout_moist')
nc2 = Dataset('/p/work1/lloveras/bwave_nov/20km_files/barotropic/output/wrfout_dry')

### Read in the data
slp1 = np.zeros((len(ti),ny,nx))
dbz1 = np.zeros((len(ti),ny,nx))
tsfc1 = np.zeros((len(ti),ny,nx))
z2501 = np.zeros((len(ti),ny,nx))

slp2 = np.zeros((len(ti),ny,nx))
tsfc2 = np.zeros((len(ti),ny,nx))
z2502 = np.zeros((len(ti),ny,nx))

for i in range(len(ti)):
    slp1[i,:,:] = np.asarray(getvar(nc1,'slp',timeidx=int(ti[i]*24/6)))
    dbz1[i,:,:] = np.asarray(getvar(nc1,'mdbz',timeidx=int(ti[i]*24/6)))
    tsfc1[i,:,:] = np.asarray(getvar(nc1,'tc',timeidx=int(ti[i]*24/6)))[0,:,:]
    p1 = np.asarray(getvar(nc1,'pressure',timeidx=int(ti[i]*24/6)))
    z1 = np.asarray(getvar(nc1,'z',timeidx=int(ti[i]*24/6)))
    z2501[i,:,:] = interplevel(z1,p1,250.)
    
    slp2[i,:,:] = np.asarray(getvar(nc2,'slp',timeidx=int(ti[i]*24/6)))
    tsfc2[i,:,:] = np.asarray(getvar(nc2,'tc',timeidx=int(ti[i]*24/6)))[0,:,:]
    p2 = np.asarray(getvar(nc2,'pressure',timeidx=int(ti[i]*24/6)))
    z2 = np.asarray(getvar(nc2,'z',timeidx=int(ti[i]*24/6)))
    z2502[i,:,:] = interplevel(z2,p2,250.)

### Save the output
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/slp_barotropic_moist_20km',slp1)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/dbz_barotropic_20km',dbz1)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/tsfc_barotropic_moist_20km',tsfc1)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/z250_barotropic_moist_20km',z2501)

np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/slp_barotropic_dry_20km',slp2)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/tsfc_barotropic_dry_20km',tsfc2)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/z250_barotropic_dry_20km',z2502)
