import numpy as np
import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from bwave_ideal_wrf import wrf_fields


CP = 1005.7
RD = 287.04
P0 = 1000.
TR = 300.
LV = 2.501e6
EPS = 1.

nx_coamps = 250
ny_coamps = 225
dx_coamps = 32000.
nz_coamps = 100
nx_wrf = 2000
ny_wrf = 1800
nz_wrf = 100
dx_wrf = 4000.
cent_length = 2500000.
shape_coamps = [nz_coamps, ny_coamps, nx_coamps]
shape_wrf = [nz_wrf, ny_wrf, nx_wrf]

slp_ctl_0h = np.load('/p/work1/lloveras/adj_4km/processed/2d_fields/data_48h/slp_ctl_0h.npy')
tlm_p = np.reshape(np.fromfile('/p/work1/lloveras/adj_4km/tlm_files/2207070000ZPP1.GRD1', dtype='>f4')[:-200],shape_coamps).astype(float)[0,:,:]

x_coamps = np.linspace(int(dx_coamps),int(nx_coamps*dx_coamps),int(nx_coamps))
x_wrf = np.linspace(int(dx_wrf),int(nx_wrf*dx_wrf),int(nx_wrf))

coamps_temp = np.zeros([ny_coamps, nx_wrf])
for j in range(ny_coamps):
    for i in range(nx_wrf):
        coamps_temp[j,i] = wrf_fields.interp_0(tlm_p[j,:], x_coamps, x_wrf[i], nx_coamps)
        
y_coamps = np.linspace(int(dx_coamps),int(ny_coamps*dx_coamps),int(ny_coamps))
y_wrf = np.linspace(int(dx_wrf),int(ny_wrf*dx_wrf),int(ny_wrf))

coamps_4km = np.zeros([ny_wrf, nx_wrf])
for j in range(ny_wrf):
    for i in range(nx_wrf):
        coamps_4km[j,i] = wrf_fields.interp_0(coamps_temp[:,i], y_coamps, y_wrf[j], ny_coamps)

minp_ind = np.unravel_index(slp_ctl_0h.argmin(), slp_ctl_0h.shape)
cent_ind =  int(cent_length/dx_wrf)
rol_val = -(cent_ind - minp_ind[-1])

coamps_4km_roll = np.roll(coamps_4km, rol_val, axis=-1)

np.save('/p/work1/lloveras/adj_4km/tlm_files/tlm_p_4km',coamps_4km_roll)
