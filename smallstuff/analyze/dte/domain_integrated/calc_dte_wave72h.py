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

### Input
nc1 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_72h/wrfout_d01_2021-01-04_00_00_00')
nc2 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/wave_72h/wrfout_d01_2021-01-04_00_00_00')
nc_in = Dataset('/p/work1/lloveras/adj_4km/wrf_output/wrfin_ctl')

### Output
dte_file = '/p/work1/lloveras/adj_4km/processed/dte/domain_integrated/normalized/long_runs/dte_wave_72h'

### Analysis
nt = 17
dte_vals = np.zeros(nt)
for ti in range(nt):
    dte_vals[ti] = calc.dte(nc1,nc2,nc_in,ti)

np.save(dte_file,dte_vals)