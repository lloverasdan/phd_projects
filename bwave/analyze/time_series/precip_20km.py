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
nc1 = Dataset('/p/work1/lloveras/bwave_nov/20km_files/standard/output/wrfout_moist')
nc2 = Dataset('/p/work1/lloveras/bwave_nov/20km_files/surface/output/wrfout_moist')
nc3 = Dataset('/p/work1/lloveras/bwave_nov/20km_files/barotropic/output/wrfout_moist')

### Output
avgfile1 = '/p/work1/lloveras/bwave_nov/processed/precip_avg/standard_20km'
areafile1 = '/p/work1/lloveras/bwave_nov/processed/precip_area/standard_20km'

avgfile2 = '/p/work1/lloveras/bwave_nov/processed/precip_avg/surface_20km'
areafile2 = '/p/work1/lloveras/bwave_nov/processed/precip_area/surface_20km'

avgfile3 = '/p/work1/lloveras/bwave_nov/processed/precip_avg/barotropic_20km'
areafile3 = '/p/work1/lloveras/bwave_nov/processed/precip_area/barotropic_20km'

### Analysis
nt = 32
avg1 = np.zeros(nt)
area1 = np.zeros(nt)
avg2 = np.zeros(nt)
area2 = np.zeros(nt)
avg3 = np.zeros(nt)
area3 = np.zeros(nt)
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
    
    temp2 = ((np.asarray(getvar(nc2,'RAINC',timeidx=ti)) +\
             np.asarray(getvar(nc2,'RAINNC',timeidx=ti)) +\
             np.asarray(getvar(nc2,'SNOWNC',timeidx=ti)) +\
             np.asarray(getvar(nc2,'HAILNC',timeidx=ti)) +\
             np.asarray(getvar(nc2,'GRAUPELNC',timeidx=ti))) -\
            (np.asarray(getvar(nc2,'RAINC',timeidx=ti-1)) +\
             np.asarray(getvar(nc2,'RAINNC',timeidx=ti-1)) +\
             np.asarray(getvar(nc2,'SNOWNC',timeidx=ti-1)) +\
             np.asarray(getvar(nc2,'HAILNC',timeidx=ti-1)) +\
             np.asarray(getvar(nc2,'GRAUPELNC',timeidx=ti-1))))/6
    
    area2[ti-1] = np.count_nonzero(temp2 > 0.)
    temp2[temp2 <= 0] = np.nan
    avg2[ti-1] = np.nan_to_num(np.nanmean(temp2))
    
    
    temp3 = ((np.asarray(getvar(nc3,'RAINC',timeidx=ti)) +\
             np.asarray(getvar(nc3,'RAINNC',timeidx=ti)) +\
             np.asarray(getvar(nc3,'SNOWNC',timeidx=ti)) +\
             np.asarray(getvar(nc3,'HAILNC',timeidx=ti)) +\
             np.asarray(getvar(nc3,'GRAUPELNC',timeidx=ti))) -\
            (np.asarray(getvar(nc3,'RAINC',timeidx=ti-1)) +\
             np.asarray(getvar(nc3,'RAINNC',timeidx=ti-1)) +\
             np.asarray(getvar(nc3,'SNOWNC',timeidx=ti-1)) +\
             np.asarray(getvar(nc3,'HAILNC',timeidx=ti-1)) +\
             np.asarray(getvar(nc3,'GRAUPELNC',timeidx=ti-1))))/6
    
    area3[ti-1] = np.count_nonzero(temp3 > 0.)
    temp3[temp3 <= 0] = np.nan
    avg3[ti-1] = np.nan_to_num(np.nanmean(temp3))

np.save(avgfile1,avg1)
np.save(areafile1,area1)

np.save(avgfile2,avg2)
np.save(areafile2,area2)

np.save(avgfile3,avg3)
np.save(areafile3,area3)
