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
nc1 = Dataset('/p/work1/lloveras/bwave_nov/4km_files/output/wrfout_standard')

### Output
avgfile1 = '/p/work1/lloveras/bwave_nov/processed/precip_avg/standard_4km'
areafile1 = '/p/work1/lloveras/bwave_nov/processed/precip_area/standard_4km'

### Analysis
nt = 32
avg1 = np.zeros(nt)
area1 = np.zeros(nt)
for ti in range(1,nt+1):
    temp1 = ((np.asarray(getvar(nc1,'RAINC',timeidx=ti)) +\
             np.asarray(getvar(nc1,'RAINNC',timeidx=ti)) +\
             np.asarray(getvar(nc1,'SNOWNC',timeidx=ti)) +\
             np.asarray(getvar(nc1,'HAILNC',timeidx=ti)) +\
             np.asarray(getvar(nc1,'GRAUPELNC',timeidx=ti))) -\
            (np.asarray(getvar(nc1,'RAINC',timeidx=ti-1)) +\
             np.asarray(getvar(nc1,'RAINNC',timeidx=ti-1)) +\
             np.asarray(getvar(nc1,'SNOWNC',timeidx=ti-1)) +\
             np.asarray(getvar(nc1,'HAILNC',timeidx=ti-1)) +\
             np.asarray(getvar(nc1,'GRAUPELNC',timeidx=ti-1))))/6
    
    area1[ti-1] = np.count_nonzero(temp1 > 0.)
    temp1[temp1 <= 0] = np.nan
    avg1[ti-1] = np.nan_to_num(np.nanmean(temp1))

np.save(avgfile1,avg1)
np.save(areafile1,area1)