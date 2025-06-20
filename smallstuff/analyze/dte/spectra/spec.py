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
nc1 = Dataset('/p/work1/lloveras/adj_4km/long_runs/ctl_48h/wrfout_d01_2021-01-05_00_00_00')
nx = 2000
ny = 3600
nz = 90
res = 4000
nmax = int(np.ceil(np.sqrt(2)*np.maximum(nx/2,ny/2)))

### Output
spec_file = '/p/work1/lloveras/adj_4km/processed_data/long_runs/spec_ctl_48h_144h'

### Analysis
nt = 2
spec = np.zeros((nt,nz,nmax))
for ti in range(nt):
    spec[ti,:,:] = calc.te_back_2d(nc1,res,timeid=int((ti+1)*4))

np.save(spec_file,spec)
