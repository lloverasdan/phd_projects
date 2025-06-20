#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from functions import calc
from netCDF4 import Dataset
import numpy as np
from wrf import getvar

### Input
nc1 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_72h/wrfout_d01_2021-01-04_00_00_00')
nc2 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/adj_72h_100/wrfout_d01_2021-01-04_00_00_00')

### Output
slp_file = '/p/work1/lloveras/adj_4km/processed/slp/mag_diffs/mag_slp_adj_72h_100'

### Analysis
nt = 17
slp_vals = np.zeros(nt)
for ti in range(nt):
    slp1 = np.amin(np.asarray(getvar(nc1,'slp',timeidx=ti)))
    slp2 = np.amin(np.asarray(getvar(nc2,'slp',timeidx=ti)))
    slp_vals[ti] = slp2 - slp1

np.save(slp_file,slp_vals)
