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
nc1 = Dataset('/p/work1/lloveras/bwave_nov/4km_files/output/wrfout_standard')

### Read in the data
cape1 = np.zeros((len(ti),ny,nx))

for i in range(len(ti)):
    cape1[i,:,:] = np.asarray(getvar(nc1,'cape_2d',timeidx=int(ti[i]*24/6)))[0]

### Save the output
np.save('/p/work1/lloveras/bwave_nov/processed/2d_vars/cape_standard_4km',cape1)
