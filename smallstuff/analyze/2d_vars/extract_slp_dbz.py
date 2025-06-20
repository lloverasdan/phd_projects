#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 14:00:02 2020

@author: lloverasdan
"""

from netCDF4 import Dataset
import numpy as np
from wrf import getvar

### Load the netCDF files
nc_ctl = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/ctl_48h/wrfout_d01_2021-01-03_00_00_00')
# nc_adj10_48 = Dataset('/p/work1/lloveras/adj_4km/wrf_output/long_runs/adj_48h_neg_10/wrfout_d01_2021-01-03_00_00_00')

### Read in the data
slp_ctl = np.asarray(getvar(nc_ctl,'slp',timeidx=0))
# dbz_adj = np.asarray(getvar(nc_adj_48,'mdbz',timeidx=16))
# slp_adj10 = np.asarray(getvar(nc_adj10_48,'slp',timeidx=16))
# dbz_adj10 = np.asarray(getvar(nc_adj10_48,'mdbz',timeidx=16))

### Save the output
np.save('/p/work1/lloveras/adj_4km/processed/2d_fields/data_48h/slp_ctl_0h',slp_ctl)
# np.save('/p/work1/lloveras/adj_4km/processed/2d_fields/data_48h/mdbz_adj_neg_48h',dbz_adj)
# np.save('/p/work1/lloveras/adj_4km/processed/2d_fields/data_48h/slp_adj10_neg_48h',slp_adj10)
# np.save('/p/work1/lloveras/adj_4km/processed/2d_fields/data_48h/mdbz_adj10_neg_48h',dbz_adj10)