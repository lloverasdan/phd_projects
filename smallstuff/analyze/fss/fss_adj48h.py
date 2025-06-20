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
nc1 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_48h/wrfout_d01_2021-01-03_00_00_00')
nc2 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/adj_48h_neg_10/wrfout_d01_2021-01-03_00_00_00')
threshold = 1.0
neighborhood = 4

### Output
fss_file = '/p/work1/lloveras/adj_4km/processed/fss/fss_adj_48h_neg_10'

### Analysis
nt = 16
fss_vals = np.zeros(nt)
for ti in range(1,nt+1):
    
    observed = ((np.asarray(getvar(nc1,'RAINC',timeidx=ti)) +\
                 np.asarray(getvar(nc1,'RAINNC',timeidx=ti)) +\
                 np.asarray(getvar(nc1,'SNOWNC',timeidx=ti)) +\
                 np.asarray(getvar(nc1,'HAILNC',timeidx=ti)) +\
                 np.asarray(getvar(nc1,'GRAUPELNC',timeidx=ti))) -\
                (np.asarray(getvar(nc1,'RAINC',timeidx=ti-1)) +\
                 np.asarray(getvar(nc1,'RAINNC',timeidx=ti-1)) +\
                 np.asarray(getvar(nc1,'SNOWNC',timeidx=ti-1)) +\
                 np.asarray(getvar(nc1,'HAILNC',timeidx=ti-1)) +\
                 np.asarray(getvar(nc1,'GRAUPELNC',timeidx=ti-1))))/3
    
    modeled = ((np.asarray(getvar(nc2,'RAINC',timeidx=ti)) +\
                 np.asarray(getvar(nc2,'RAINNC',timeidx=ti)) +\
                 np.asarray(getvar(nc2,'SNOWNC',timeidx=ti)) +\
                 np.asarray(getvar(nc2,'HAILNC',timeidx=ti)) +\
                 np.asarray(getvar(nc2,'GRAUPELNC',timeidx=ti))) -\
                (np.asarray(getvar(nc2,'RAINC',timeidx=ti-1)) +\
                 np.asarray(getvar(nc2,'RAINNC',timeidx=ti-1)) +\
                 np.asarray(getvar(nc2,'SNOWNC',timeidx=ti-1)) +\
                 np.asarray(getvar(nc2,'HAILNC',timeidx=ti-1)) +\
                 np.asarray(getvar(nc2,'GRAUPELNC',timeidx=ti-1))))/3
    
    fss_vals[ti-1] = calc.fss(modeled, observed, threshold, neighborhood)

np.save(fss_file,fss_vals)
