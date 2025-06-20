#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

from netCDF4 import Dataset
import numpy as np
from wrf import getvar, interplevel

### Load the netCDF files
nc1 = Dataset('/p/work1/lloveras/bwave_nov/4km_files/output/wrfout_surface')

y = 800
x1 = 25
x2 = 275

p1 = np.asarray(getvar(nc1,'pressure',timeidx=24))[:,y,x1:x2]
z1 = np.asarray(getvar(nc1,'z',timeidx=24))[:,y,x1:x2]
dbz1 = np.asarray(getvar(nc1,'dbz',timeidx=24))[:,y,x1:x2]
t1 = np.asarray(getvar(nc1,'th',timeidx=24))[:,y,x1:x2]

np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/p_cross_surface_4km',p1)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/z_cross_surface_4km',z1)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/dbz_cross_surface_4km',dbz1)
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/th_cross_surface_4km',t1)