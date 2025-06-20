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
nc1 = Dataset('/p/work1/lloveras/bwave_nov/4km_files/output/wrfout_surface')

### Output
file1 = '/p/work1/lloveras/bwave_nov/processed/slp_time/surface_4km'

### Analysis
nt = 33
nx = 2000
ny = 1800

pb1 = np.tile(np.expand_dims(np.asarray(getvar(nc1,'slp'))[:,0],-1),nx)

slp1 = np.zeros(nt)
for ti in range(nt):
    slp1[ti] = np.amin(np.asarray(getvar(nc1,'slp',timeidx=ti)) - pb1)

np.save(file1,slp1)
