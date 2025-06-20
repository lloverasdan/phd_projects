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
nx = 2000
ny = 1800

### Load the netCDF files
nc1 = Dataset('/p/work1/lloveras/bwave_nov/4km_files/output/wrfout_barotropic')

### Read in the data
slp1 = np.zeros((len(ti),ny,nx))
dbz1 = np.zeros((len(ti),ny,nx))
tsfc1 = np.zeros((len(ti),ny,nx))
z2501 = np.zeros((len(ti),ny,nx))

for i in range(len(ti)):
    slp1[i,:,:] = np.asarray(getvar(nc1,'slp',timeidx=int(ti[i]*24/6)))
    dbz1[i,:,:] = np.asarray(getvar(nc1,'mdbz',timeidx=int(ti[i]*24/6)))
    tsfc1[i,:,:] = np.asarray(getvar(nc1,'tc',timeidx=int(ti[i]*24/6)))[0,:,:]

    p1 = np.asarray(getvar(nc1,'pressure',timeidx=int(ti[i]*24/6)))
    z1 = np.asarray(getvar(nc1,'z',timeidx=int(ti[i]*24/6)))
    z2501[i,:,:] = interplevel(z1,p1,250.)

### Save the output
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/slp_barotropic_4km',slp1)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/dbz_barotropic_4km',dbz1)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/tsfc_barotropic_4km',tsfc1)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/z250_barotropic_4km',z2501)
