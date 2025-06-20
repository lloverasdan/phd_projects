#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

from netCDF4 import Dataset
from wrf import getvar
import numpy as np

### Input
nc1 = Dataset('/p/work1/lloveras/bwave/20km_files/barotropic/output/wrfout_moist_rh75')
nc2 = Dataset('/p/work1/lloveras/bwave/20km_files/standard/output/wrfout_moist_rh75')
nc3 = Dataset('/p/work1/lloveras/bwave/20km_files/surface/output/wrfout_moist_rh75')
# nc4 = Dataset('/p/work1/lloveras/bwave_nov/20km_files/surface/output/wrfout_dry')
# nc5 = Dataset('/p/work1/lloveras/bwave_nov/20km_files/barotropic/output/wrfout_moist')
# nc6 = Dataset('/p/work1/lloveras/bwave_nov/20km_files/barotropic/output/wrfout_dry')

### Output
file1 = '/p/work1/lloveras/bwave/processed/slp_time/barotropic_moist_rh75_20km'
file2 = '/p/work1/lloveras/bwave/processed/slp_time/standard_moist_rh75_20km'
file3 = '/p/work1/lloveras/bwave/processed/slp_time/surface_moist_rh75_20km'
# file4 = '/p/work1/lloveras/bwave_nov/processed/slp_time/surface_dry_20km'
# file5 = '/p/work1/lloveras/bwave_nov/processed/slp_time/barotropic_moist_20km'
# file6 = '/p/work1/lloveras/bwave_nov/processed/slp_time/barotropic_dry_20km'

### Analysis
nt = 33
nx = 400
ny = 360

pb1 = np.tile(np.expand_dims(np.asarray(getvar(nc1,'slp'))[:,0],-1),nx)
pb2 = np.tile(np.expand_dims(np.asarray(getvar(nc2,'slp'))[:,0],-1),nx)
pb3 = np.tile(np.expand_dims(np.asarray(getvar(nc3,'slp'))[:,0],-1),nx)
# pb4 = np.tile(np.expand_dims(np.asarray(getvar(nc4,'slp'))[:,0],-1),nx)
# pb5 = np.tile(np.expand_dims(np.asarray(getvar(nc5,'slp'))[:,0],-1),nx)
# pb6 = np.tile(np.expand_dims(np.asarray(getvar(nc6,'slp'))[:,0],-1),nx)

slp1 = np.zeros(nt)
slp2 = np.zeros(nt)
slp3 = np.zeros(nt)
# slp4 = np.zeros(nt)
# slp5 = np.zeros(nt)
# slp6 = np.zeros(nt)
for ti in range(nt):
    slp1[ti] = np.amin(np.asarray(getvar(nc1,'slp',timeidx=ti)) - pb1)
    slp2[ti] = np.amin(np.asarray(getvar(nc2,'slp',timeidx=ti)) - pb2)
    slp3[ti] = np.amin(np.asarray(getvar(nc3,'slp',timeidx=ti)) - pb3)
    # slp4[ti] = np.amin(np.asarray(getvar(nc4,'slp',timeidx=ti)) - pb4)
    # slp5[ti] = np.amin(np.asarray(getvar(nc5,'slp',timeidx=ti)) - pb5)
    # slp6[ti] = np.amin(np.asarray(getvar(nc6,'slp',timeidx=ti)) - pb6)    

np.save(file1,slp1)
np.save(file2,slp2)
np.save(file3,slp3)
# np.save(file4,slp4)
# np.save(file5,slp5)
# np.save(file6,slp6)