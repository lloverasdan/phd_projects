#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 13:49:29 2019

@author: Daniel Lloveras

Defines the function for writing the wrfinput_d01 netCDF file to be used by
WRF at runtime
"""

import numpy as np

#%% Constants

RD = 287
FF = 1e-4

def write(ncfile, nx, ny, nz, hres_m, title_str, time_str, u, v, t, ph, phb,\
          t_init, mu, mub, p, pb, fnm, fnp, rdnw, rdn, dnw, dn, t_base,\
          cfn, cfn1, rdx, rdy, cf1, cf2, cf3, tsk, u_base, v_base, qv_base,\
          tmn, moist, znw, znu, diff_opt, km_opt, damp_opt, dampcoef, khdif,\
          kvdif, mp_physics, ra_lw_physics, ra_sw_physics, sf_sfclay_physics,\
          sf_surface_physics, bl_pbl_physics, cu_physics, sf_lake_physics,\
          surface_input_source, hypsometric_opt, dt, num_land_cat,\
          num_soil_layers, num_soil_cat, spec_bdy_width, ptop):
    
    ncfile.createDimension('Time',None)
    ncfile.createDimension('DateStrLen',19)
    ncfile.createDimension('west_east',nx)
    ncfile.createDimension('south_north',ny)
    ncfile.createDimension('bottom_top',nz)
    ncfile.createDimension('bottom_top_stag',nz+1)
    ncfile.createDimension('soil_layers_stag',num_soil_layers)
    ncfile.createDimension('west_east_stag',nx+1)
    ncfile.createDimension('south_north_stag',ny+1)
    ncfile.createDimension('DIM0009',spec_bdy_width)
    ncfile.createDimension('land_cat_stag',num_land_cat)
    ncfile.createDimension('soil_cat_stag',num_soil_cat)
    ncfile.createDimension('num_ext_model_couple_dom_stag',1)
    
    ncfile.TITLE = title_str
    ncfile.START_DATE = time_str
    ncfile.SIMULATION_START_DATE = time_str
    ncfile.WEST_EAST_GRID_DIMENSION = nx+1
    ncfile.SOUTH_NORTH_GRID_DIMENSION = ny+1
    ncfile.BOTTOM_TOP_GRID_DIMENSION = nz+1
    ncfile.DX = hres_m
    ncfile.DY = hres_m
    ncfile.GRIDTYPE = "C"
    ncfile.DIFF_OPT = diff_opt
    ncfile.KM_OPT = km_opt
    ncfile.DAMP_OPT = damp_opt
    ncfile.DAMPCOEF = dampcoef
    ncfile.KHDIF = khdif
    ncfile.KVDIF = kvdif
    ncfile.MP_PHYSICS = mp_physics
    ncfile.RA_LW_PHYSICS = ra_lw_physics
    ncfile.RA_SW_PHYSICS = ra_sw_physics
    ncfile.SF_SFCLAY_PHYSICS = sf_sfclay_physics
    ncfile.SF_SURFACE_PHYSICS = sf_surface_physics
    ncfile.BL_PBL_PHYSICS = bl_pbl_physics
    ncfile.CU_PHYSICS = cu_physics
    ncfile.SF_LAKE_PHYSICS = sf_lake_physics
    ncfile.SURFACE_INPUT_SOURCE = surface_input_source
    ncfile.SST_UPDATE = 0
    ncfile.GRID_FDDA = 0
    ncfile.GFDDA_INTERVAL_M = 0
    ncfile.GFDDA_END_H = 0
    ncfile.GRID_SFDDA = 0
    ncfile.SGFDDA_INTERVAL_M = 0
    ncfile.SGFDDA_END_H = 0
    ncfile.HYPSOMETRIC_OPT = hypsometric_opt
    ncfile.WEST_EAST_PATCH_START_UNSTAG = 1
    ncfile.WEST_EAST_PATCH_END_UNSTAG = nx
    ncfile.WEST_EAST_PATCH_START_STAG = 1
    ncfile.WEST_EAST_PATCH_END_STAG = nx+1
    ncfile.SOUTH_NORTH_PATCH_START_UNSTAG = 1
    ncfile.SOUTH_NORTH_PATCH_END_UNSTAG = ny
    ncfile.SOUTH_NORTH_PATCH_START_STAG = 1
    ncfile.SOUTH_NORTH_PATCH_END_STAG = ny+1
    ncfile.BOTTOM_TOP_PATCH_START_UNSTAG = 1
    ncfile.BOTTOM_TOP_PATCH_END_UNSTAG = nz
    ncfile.BOTTOM_TOP_PATCH_START_STAG = 1
    ncfile.BOTTOM_TOP_PATCH_END_STAG = nz+1
    ncfile.GRID_ID = 1
    ncfile.PARENT_ID = 0
    ncfile.I_PARENT_START = 0
    ncfile.J_PARENT_START = 0
    ncfile.PARENT_GRID_RATIO = 1
    ncfile.DT = dt
    ncfile.CEN_LAT = 0.
    ncfile.CEN_LON = 0.
    ncfile.TRUELAT1 = 0.
    ncfile.TRUELAT2 = 0.
    ncfile.MOAD_CEN_LAT = 0.
    ncfile.STAND_LON = 0.
    ncfile.POLE_LAT = 0.
    ncfile.POLE_LON = 0.
    ncfile.GMT = 0.
    ncfile.JULYR = 0
    ncfile.JULDAY = 1
    ncfile.MAP_PROJ = 0
    ncfile.MAP_PROJ_CHAR = "Cartesian"
    ncfile.MMINLU = "USGS"
    ncfile.NUM_LAND_CAT = num_land_cat
    ncfile.ISWATER = 0
    ncfile.ISLAKE = 0
    ncfile.ISICE = 0
    ncfile.ISURBAN = 0
    ncfile.ISOILWATER = 0
    
    ### Rename attributes to include hyphens
    ncfile.renameAttribute('WEST_EAST_GRID_DIMENSION','WEST-EAST_GRID_DIMENSION')
    ncfile.renameAttribute('SOUTH_NORTH_GRID_DIMENSION','SOUTH-NORTH_GRID_DIMENSION')
    ncfile.renameAttribute('BOTTOM_TOP_GRID_DIMENSION','BOTTOM-TOP_GRID_DIMENSION')
    ncfile.renameAttribute('WEST_EAST_PATCH_START_UNSTAG','WEST-EAST_PATCH_START_UNSTAG')
    ncfile.renameAttribute('WEST_EAST_PATCH_END_UNSTAG','WEST-EAST_PATCH_END_UNSTAG')
    ncfile.renameAttribute('WEST_EAST_PATCH_START_STAG','WEST-EAST_PATCH_START_STAG')
    ncfile.renameAttribute('WEST_EAST_PATCH_END_STAG','WEST-EAST_PATCH_END_STAG')
    ncfile.renameAttribute('SOUTH_NORTH_PATCH_START_UNSTAG','SOUTH-NORTH_PATCH_START_UNSTAG')
    ncfile.renameAttribute('SOUTH_NORTH_PATCH_END_UNSTAG','SOUTH-NORTH_PATCH_END_UNSTAG')
    ncfile.renameAttribute('SOUTH_NORTH_PATCH_START_STAG','SOUTH-NORTH_PATCH_START_STAG')
    ncfile.renameAttribute('SOUTH_NORTH_PATCH_END_STAG','SOUTH-NORTH_PATCH_END_STAG')
    ncfile.renameAttribute('BOTTOM_TOP_PATCH_START_UNSTAG','BOTTOM-TOP_PATCH_START_UNSTAG')
    ncfile.renameAttribute('BOTTOM_TOP_PATCH_END_UNSTAG','BOTTOM-TOP_PATCH_END_UNSTAG')
    ncfile.renameAttribute('BOTTOM_TOP_PATCH_START_STAG','BOTTOM-TOP_PATCH_START_STAG')
    ncfile.renameAttribute('BOTTOM_TOP_PATCH_END_STAG','BOTTOM-TOP_PATCH_END_STAG')
    
    ### Time
    Times = ncfile.createVariable('Times','S1',('Time','DateStrLen'))
    
    ### Latitude
    XLAT = ncfile.createVariable('XLAT','f4',('Time','south_north','west_east'))
    XLAT.FieldType = 104
    XLAT.MemoryOrder = "XY "
    XLAT.description = "LATITUDE, SOUTH IS NEGATIVE"
    XLAT.units = "degree_north"
    XLAT.stagger = ""
    
    ### Longitude
    XLONG = ncfile.createVariable('XLONG','f4',('Time','south_north','west_east'))
    XLONG.FieldType = 104
    XLONG.MemoryOrder = "XY "
    XLONG.description = "LONGITUDE, WEST IS NEGATIVE"
    XLONG.units = "degree_east"
    XLONG.stagger = ""
    
    ### Land use category
    LU_INDEX = ncfile.createVariable('LU_INDEX','f4',('Time','south_north','west_east'))
    LU_INDEX.FieldType = 104
    LU_INDEX.MemoryOrder = "XY "
    LU_INDEX.description = "LAND USE CATEGORY"
    LU_INDEX.units = ""
    LU_INDEX.stagger = ""
    LU_INDEX.coordinates = "XLONG XLAT"
    
    ### Eta values on half (mass) levels
    ZNU = ncfile.createVariable('ZNU','f4',('Time','bottom_top'))
    ZNU.FieldType = 104
    ZNU.MemoryOrder = "Z  "
    ZNU.description = "eta values on half (mass) levels"
    ZNU.units = ""
    ZNU.stagger = ""
    
    ### Eta values on full (w) levels
    ZNW = ncfile.createVariable('ZNW','f4',('Time','bottom_top_stag'))
    ZNW.FieldType = 104
    ZNW.MemoryOrder = "Z  "
    ZNW.description = "eta values on full (w) levels"
    ZNW.units = ""
    ZNW.stagger = "Z"
    
    ### Depths of centers of soil layers
    ZS = ncfile.createVariable('ZS','f4',('Time','soil_layers_stag'))
    ZS.FieldType = 104
    ZS.MemoryOrder = "Z  "
    ZS.description = "DEPTHS OF CENTERS OF SOIL LAYERS"
    ZS.units = "m"
    ZS.stagger = "Z"
    
    ### Thickness of soil layers
    DZS = ncfile.createVariable('DZS','f4',('Time','soil_layers_stag'))
    DZS.FieldType = 104
    DZS.MemoryOrder = "Z  "
    DZS.description = "THICKNESS OF SOIL LAYERS"
    DZS.units = "m"
    DZS.stagger = "Z"
    
    ### Variance of subgrid-scale orography
    VAR_SSO = ncfile.createVariable('VAR_SSO','f4',('Time','south_north','west_east'))
    VAR_SSO.FieldType = 104
    VAR_SSO.MemoryOrder = "XY "
    VAR_SSO.description = "variance of subgrid-scale orography"
    VAR_SSO.units = "m2"
    VAR_SSO.stagger = ""
    VAR_SSO.coordinates = "XLONG XLAT"
    
    ### X-wind component
    U = ncfile.createVariable('U','f4',('Time','bottom_top','south_north','west_east_stag'))
    U.FieldType = 104
    U.MemoryOrder = "XYZ"
    U.description = "x-wind component"
    U.units = "m s-1"
    U.stagger = "X"
    U.coordinates = "XLONG_U XLAT_U"
    
    ### Y-wind component
    V = ncfile.createVariable('V','f4',('Time','bottom_top','south_north_stag','west_east'))
    V.FieldType = 104
    V.MemoryOrder = "XYZ"
    V.description = "y-wind component"
    V.units = "m s-1"
    V.stagger = "Y"
    V.coordinates = "XLONG_V XLAT_V"
    
    ### Z-wind component
    W = ncfile.createVariable('W','f4',('Time','bottom_top_stag','south_north','west_east'))
    W.FieldType = 104
    W.MemoryOrder = "XYZ"
    W.description = "z-wind component"
    W.units = "m s-1"
    W.stagger = "Z"
    W.coordinates = "XLONG XLAT"
    
    ### Perturbation geopotential
    PH = ncfile.createVariable('PH','f4',('Time','bottom_top_stag','south_north','west_east'))
    PH.FieldType = 104
    PH.MemoryOrder = "XYZ"
    PH.description = "perturbation geopotential"
    PH.units = "m2 s-2"
    PH.stagger = "Z"
    PH.coordinates = "XLONG XLAT"
    
    ### Base-state geopotential
    PHB = ncfile.createVariable('PHB','f4',('Time','bottom_top_stag','south_north','west_east'))
    PHB.FieldType = 104
    PHB.MemoryOrder = "XYZ"
    PHB.description = "base-state geopotential"
    PHB.units = "m2 s-2"
    PHB.stagger = "Z"
    PHB.coordinates = "XLONG XLAT"
    
    ### Perturbation potential temperature (theta - t0)
    T = ncfile.createVariable('T','f4',('Time','bottom_top','south_north','west_east'))
    T.FieldType = 104
    T.MemoryOrder = "XYZ"
    T.description = "perturbation potential temperature (theta-t0)"
    T.units = "K"
    T.stagger = ""
    T.coordinates = "XLONG XLAT"
    
    ### Initial potential temperature
    T_INIT = ncfile.createVariable('T_INIT','f4',('Time','bottom_top','south_north','west_east'))
    T_INIT.FieldType = 104
    T_INIT.MemoryOrder = "XYZ"
    T_INIT.description = "initial potential temperature"
    T_INIT.units = "K"
    T_INIT.stagger = ""
    T_INIT.coordinates = "XLONG XLAT"
    
    ### Perturbation dy air mass in column
    MU = ncfile.createVariable('MU','f4',('Time','south_north','west_east'))
    MU.FieldType = 104
    MU.MemoryOrder = "XY "
    MU.description = "perturbation dry air mass in column"
    MU.units = "Pa"
    MU.stagger = ""
    MU.coordinates = "XLONG XLAT"
    
    ### Base state dry air mass in column
    MUB = ncfile.createVariable('MUB','f4',('Time','south_north','west_east'))
    MUB.FieldType = 104
    MUB.MemoryOrder = "XY "
    MUB.description = "base state dry air mass in column"
    MUB.units = "Pa"
    MUB.stagger = ""
    MUB.coordinates = "XLONG XLAT"
    
    ### Perturbation pressure
    P = ncfile.createVariable('P','f4',('Time','bottom_top','south_north','west_east'))
    P.FieldType = 104
    P.MemoryOrder = "XYZ"
    P.description = "perturbation pressure"
    P.units = "Pa"
    P.stagger = ""
    P.coordinates = "XLONG XLAT"
    
    ### Base state pressure
    PB = ncfile.createVariable('PB','f4',('Time','bottom_top','south_north','west_east'))
    PB.FieldType = 104
    PB.MemoryOrder = "XYZ"
    PB.description = "BASE STATE PRESSURE"
    PB.units = "Pa"
    PB.stagger = ""
    PB.coordinates = "XLONG XLAT"
    
    ### Upper weight for vertical stretching
    FNM = ncfile.createVariable('FNM','f4',('Time','bottom_top'))
    FNM.FieldType = 104
    FNM.MemoryOrder = "Z  "
    FNM.description = "upper weight for vertical stretching"
    FNM.units = ""
    FNM.stagger = ""
    
    ### Lower weight for vertical stretching
    FNP = ncfile.createVariable('FNP','f4',('Time','bottom_top'))
    FNP.FieldType = 104
    FNP.MemoryOrder = "Z  "
    FNP.description = "lower weight for vertical stretching"
    FNP.units = ""
    FNP.stagger = ""
    
    ### Inverse d(eta) values between full (w) levels
    RDNW = ncfile.createVariable('RDNW','f4',('Time','bottom_top'))
    RDNW.FieldType = 104
    RDNW.MemoryOrder = "Z  "
    RDNW.description = "inverse d(eta) values between full (w) levels"
    RDNW.units = ""
    RDNW.stagger = ""
    
    ### Inverse d(eta) values between half (mass) levels
    RDN = ncfile.createVariable('RDN','f4',('Time','bottom_top'))
    RDN.FieldType = 104
    RDN.MemoryOrder = "Z  "
    RDN.description = "inverse d(eta) values between half (mass) levels"
    RDN.units = ""
    RDN.stagger = ""
    
    ### d(eta) values between full (w) levels
    DNW = ncfile.createVariable('DNW','f4',('Time','bottom_top'))
    DNW.FieldType = 104
    DNW.MemoryOrder = "Z  "
    DNW.description = "d(eta) values between full (w) levels"
    DNW.units = ""
    DNW.stagger = ""
    
    ### d(eta) values between half (mass) levels
    DN = ncfile.createVariable('DN','f4',('Time','bottom_top'))
    DN.FieldType = 104
    DN.MemoryOrder = "Z  "
    DN.description = "d(eta) values between half (mass) levels"
    DN.units = ""
    DN.stagger = ""
    
    ### Base state T in idealized cases
    T_BASE = ncfile.createVariable('T_BASE','f4',('Time','bottom_top'))
    T_BASE.FieldType = 104
    T_BASE.MemoryOrder = "Z  "
    T_BASE.description = "BASE STATE T IN IDEALIZED CASES"
    T_BASE.units = "K"
    T_BASE.stagger = ""
    
    ### Extrapolation constant
    CFN = ncfile.createVariable('CFN','f4','Time')
    CFN.FieldType = 104
    CFN.MemoryOrder = "0  "
    CFN.description = "extrapolation constant"
    CFN.units = ""
    CFN.stagger = ""
    
    ### 1 - extrapolation constant
    CFN1 = ncfile.createVariable('CFN1','f4','Time')
    CFN1.FieldType = 104
    CFN1.MemoryOrder = "0  "
    CFN1.description = "extrapolation constant"
    CFN1.units = ""
    CFN1.stagger = ""
    
    ### Step number
    STEP_NUMBER = ncfile.createVariable('STEP_NUMBER','i4','Time')
    STEP_NUMBER.FieldType = 106
    STEP_NUMBER.MemoryOrder = "0  "
    STEP_NUMBER.description = ""
    STEP_NUMBER.units = "-"
    STEP_NUMBER.stagger = ""
    
    ### Hydrostatic pressure
    P_HYD = ncfile.createVariable('P_HYD','f4',('Time','bottom_top','south_north','west_east'))
    P_HYD.FieldType = 104
    P_HYD.MemoryOrder = "XYZ"
    P_HYD.description = "hydrostatic pressure"
    P_HYD.units = "Pa"
    P_HYD.stagger = ""
    P_HYD.coordinates = "XLONG XLAT"
    
    ### 2 meter qv
    Q2 = ncfile.createVariable('Q2','f4',('Time','south_north','west_east'))
    Q2.FieldType = 104
    Q2.MemoryOrder = "XY "
    Q2.description = "QV at 2 M"
    Q2.units = "kg kg-1"
    Q2.stagger = ""
    Q2.coordinates = "XLONG XLAT"
    
    ### 2 meter temperature
    T2 = ncfile.createVariable('T2','f4',('Time','south_north','west_east'))
    T2.FieldType = 104
    T2.MemoryOrder = "XY "
    T2.description = "TEMP at 2 M"
    T2.units = "K"
    T2.stagger = ""
    T2.coordinates = "XLONG XLAT"
    
    ### 2 meter potential temperature
    TH2 = ncfile.createVariable('TH2','f4',('Time','south_north','west_east'))
    TH2.FieldType = 104
    TH2.MemoryOrder = "XY "
    TH2.description = "POT TEMP at 2 M"
    TH2.units = "K"
    TH2.stagger = ""
    TH2.coordinates = "XLONG XLAT"
    
    ### Surface pressure
    PSFC = ncfile.createVariable('PSFC','f4',('Time','south_north','west_east'))
    PSFC.FieldType = 104
    PSFC.MemoryOrder = "XY "
    PSFC.description = "SFC PRESSURE"
    PSFC.units = "Pa"
    PSFC.stagger = ""
    PSFC.coordinates = "XLONG XLAT"
    
    ### 10 meter U
    U10 = ncfile.createVariable('U10','f4',('Time','south_north','west_east'))
    U10.FieldType = 104
    U10.MemoryOrder = "XY "
    U10.description = "U at 10 M"
    U10.units = "m s-1"
    U10.stagger = ""
    U10.coordinates = "XLONG XLAT"
    
    ### 10 meter V
    V10 = ncfile.createVariable('V10','f4',('Time','south_north','west_east'))
    V10.FieldType = 104
    V10.MemoryOrder = "XY "
    V10.description = "V at 10 M"
    V10.units = "m s-1"
    V10.stagger = ""
    V10.coordinates = "XLONG XLAT"
    
    ### Inverse x grid length
    RDX = ncfile.createVariable('RDX','f4','Time')
    RDX.FieldType = 104
    RDX.MemoryOrder = "0  "
    RDX.description = "INVERSE X GRID LENGTH"
    RDX.units = ""
    RDX.stagger = ""
    
    ### Inverse y grid length
    RDY = ncfile.createVariable('RDY','f4','Time')
    RDY.FieldType = 104
    RDY.MemoryOrder = "0  "
    RDY.description = "INVERSE Y GRID LENGTH"
    RDY.units = ""
    RDY.stagger = ""
    
    ### Small timestep
    DTS = ncfile.createVariable('DTS','f4','Time')
    DTS.FieldType = 104
    DTS.MemoryOrder = "0  "
    DTS.description = "SMALL TIMESTEP"
    DTS.units = ""
    DTS.stagger = ""
    
    ### Time weight constant for small steps
    DTSEPS = ncfile.createVariable('DTSEPS','f4','Time')
    DTSEPS.FieldType = 104
    DTSEPS.MemoryOrder = "0  "
    DTSEPS.description = "TIME WEIGHT CONSTANT FOR SMALL STEPS"
    DTSEPS.units = ""
    DTSEPS.stagger = ""
    
    ### Time weight constant for small steps
    RESM = ncfile.createVariable('RESM','f4','Time')
    RESM.FieldType = 104
    RESM.MemoryOrder = "0  "
    RESM.description = "TIME WEIGHT CONSTANT FOR SMALL STEPS"
    RESM.units = ""
    RESM.stagger = ""
    
    ### Zeta at model top
    ZETATOP = ncfile.createVariable('ZETATOP','f4','Time')
    ZETATOP.FieldType = 104
    ZETATOP.MemoryOrder = "0  "
    ZETATOP.description = "ZETA AT MODEL TOP"
    ZETATOP.units = ""
    ZETATOP.stagger = ""
    
    ### 2nd order extrapolation constant
    CF1 = ncfile.createVariable('CF1','f4','Time')
    CF1.FieldType = 104
    CF1.MemoryOrder = "0  "
    CF1.description = "2nd order extrapolation constant"
    CF1.units = ""
    CF1.stagger = ""
    
    ### 2nd order extrapolation constant
    CF2 = ncfile.createVariable('CF2','f4','Time')
    CF2.FieldType = 104
    CF2.MemoryOrder = "0  "
    CF2.description = "2nd order extrapolation constant"
    CF2.units = ""
    CF2.stagger = ""
    
    ### 2nd order extrapolation constant
    CF3 = ncfile.createVariable('CF3','f4','Time')
    CF3.FieldType = 104
    CF3.MemoryOrder = "0  "
    CF3.description = "2nd order extrapolation constant"
    CF3.units = ""
    CF3.stagger = ""
    
    ### Water vapor mixing ratio
    QVAPOR = ncfile.createVariable('QVAPOR','f4',('Time','bottom_top','south_north','west_east'))
    QVAPOR.FieldType = 104
    QVAPOR.MemoryOrder = "XYZ"
    QVAPOR.description = "Water vapor mixing ratio"
    QVAPOR.units = "kg kg-1"
    QVAPOR.stagger = ""
    QVAPOR.coordinates = "XLONG XLAT"
    
    ### Cloud water mixing ratio
    QCLOUD = ncfile.createVariable('QCLOUD','f4',('Time','bottom_top','south_north','west_east'))
    QCLOUD.FieldType = 104
    QCLOUD.MemoryOrder = "XYZ"
    QCLOUD.description = "Cloud water mixing ratio"
    QCLOUD.units = "kg kg-1"
    QCLOUD.stagger = ""
    QCLOUD.coordinates = "XLONG XLAT"
    
    ### Rain water mixing ratio
    QRAIN = ncfile.createVariable('QRAIN','f4',('Time','bottom_top','south_north','west_east'))
    QRAIN.FieldType = 104
    QRAIN.MemoryOrder = "XYZ"
    QRAIN.description = "Rain water mixing ratio"
    QRAIN.units = "kg kg-1"
    QRAIN.stagger = ""
    QRAIN.coordinates = "XLONG XLAT"
    
    ### Ice mixing ratio
    QICE = ncfile.createVariable('QICE','f4',('Time','bottom_top','south_north','west_east'))
    QICE.FieldType = 104
    QICE.MemoryOrder = "XYZ"
    QICE.description = "Ice mixing ratio"
    QICE.units = "kg kg-1"
    QICE.stagger = ""
    QICE.coordinates = "XLONG XLAT"
    
    ### Snow mixing ratio
    QSNOW = ncfile.createVariable('QSNOW','f4',('Time','bottom_top','south_north','west_east'))
    QSNOW.FieldType = 104
    QSNOW.MemoryOrder = "XYZ"
    QSNOW.description = "Snow mixing ratio"
    QSNOW.units = "kg kg-1"
    QSNOW.stagger = ""
    QSNOW.coordinates = "XLONG XLAT"
    
    ### Graupel mixing ratio
    QGRAUP = ncfile.createVariable('QGRAUP','f4',('Time','bottom_top','south_north','west_east'))
    QGRAUP.FieldType = 104
    QGRAUP.MemoryOrder = "XYZ"
    QGRAUP.description = "Graupel mixing ratio"
    QGRAUP.units = "kg kg-1"
    QGRAUP.stagger = ""
    QGRAUP.coordinates = "XLONG XLAT"
    
    ### Hail mixing ratio
    QHAIL = ncfile.createVariable('QHAIL','f4',('Time','bottom_top','south_north','west_east'))
    QHAIL.FieldType = 104
    QHAIL.MemoryOrder = "XYZ"
    QHAIL.description = "Hail mixing ratio"
    QHAIL.units = "kg kg-1"
    QHAIL.stagger = ""
    QHAIL.coordinates = "XLONG XLAT"
    
    ### Cloud water number concentration
    QNCLOUD = ncfile.createVariable('QNCLOUD','f4',('Time','bottom_top','south_north','west_east'))
    QNCLOUD.FieldType = 104
    QNCLOUD.MemoryOrder = "XYZ"
    QNCLOUD.description = "cloud water Number concentration"
    QNCLOUD.units = "  kg(-1)"
    QNCLOUD.stagger = ""
    QNCLOUD.coordinates = "XLONG XLAT"
    
    ### Rain number concentration
    QNRAIN = ncfile.createVariable('QNRAIN','f4',('Time','bottom_top','south_north','west_east'))
    QNRAIN.FieldType = 104
    QNRAIN.MemoryOrder = "XYZ"
    QNRAIN.description = "Rain Number concentration"
    QNRAIN.units = "  kg(-1)"
    QNRAIN.stagger = ""
    QNRAIN.coordinates = "XLONG XLAT"
    
    ### Ice number concentration
    QNICE = ncfile.createVariable('QNICE','f4',('Time','bottom_top','south_north','west_east'))
    QNICE.FieldType = 104
    QNICE.MemoryOrder = "XYZ"
    QNICE.description = "Ice Number concentration"
    QNICE.units = "  kg(-1)"
    QNICE.stagger = ""
    QNICE.coordinates = "XLONG XLAT"
    
    ### Snow number concentration
    QNSNOW = ncfile.createVariable('QNSNOW','f4',('Time','bottom_top','south_north','west_east'))
    QNSNOW.FieldType = 104
    QNSNOW.MemoryOrder = "XYZ"
    QNSNOW.description = "Snow Number concentration"
    QNSNOW.units = "  kg(-1)"
    QNSNOW.stagger = ""
    QNSNOW.coordinates = "XLONG XLAT"
    
    ### Graupel number concentration
    QNGRAUPEL = ncfile.createVariable('QNGRAUPEL','f4',('Time','bottom_top','south_north','west_east'))
    QNGRAUPEL.FieldType = 104
    QNGRAUPEL.MemoryOrder = "XYZ"
    QNGRAUPEL.description = "Graupel Number concentration"
    QNGRAUPEL.units = "  kg(-1)"
    QNGRAUPEL.stagger = ""
    QNGRAUPEL.coordinates = "XLONG XLAT"
    
    ### Hail number concentration
    QNHAIL = ncfile.createVariable('QNHAIL','f4',('Time','bottom_top','south_north','west_east'))
    QNHAIL.FieldType = 104
    QNHAIL.MemoryOrder = "XYZ"
    QNHAIL.description = "Hail Number concentration"
    QNHAIL.units = "  kg(-1)"
    QNHAIL.stagger = ""
    QNHAIL.coordinates = "XLONG XLAT"
    
    ### Graupel particle volume
    QVGRAUPEL = ncfile.createVariable('QVGRAUPEL','f4',('Time','bottom_top','south_north','west_east'))
    QVGRAUPEL.FieldType = 104
    QVGRAUPEL.MemoryOrder = "XYZ"
    QVGRAUPEL.description = "Hail Number concentration"
    QVGRAUPEL.units = "m(3)  kg(-1)"
    QVGRAUPEL.stagger = ""
    QVGRAUPEL.coordinates = "XLONG XLAT"
    
    ### Relaxation term for boundary zone
    FCX = ncfile.createVariable('FCX','f4',('Time','DIM0009'))
    FCX.FieldType = 104
    FCX.MemoryOrder = "C  "
    FCX.description = "RELAXATION TERM FOR BOUNDARY ZONE"
    FCX.units = ""
    FCX.stagger = ""
    
    ### 2nd relaxation term for boundary zone
    GCX = ncfile.createVariable('GCX','f4',('Time','DIM0009'))
    GCX.FieldType = 104
    GCX.MemoryOrder = "C  "
    GCX.description = "2ND RELAXATION TERM FOR BOUNDARY ZONE"
    GCX.units = ""
    GCX.stagger = ""
    
    ### Time since boundary read
    DTBC = ncfile.createVariable('DTBC','f4','Time')
    DTBC.FieldType = 104
    DTBC.MemoryOrder = "0  "
    DTBC.description = "TIME SINCE BOUNDARY READ"
    DTBC.units = ""
    DTBC.stagger = ""
    
    ### Elevation x slope
    TOPOSLPX = ncfile.createVariable('TOPOSLPX','f4',('Time','south_north','west_east'))
    TOPOSLPX.FieldType = 104
    TOPOSLPX.MemoryOrder = "XY "
    TOPOSLPX.description = "ELEVATION X SLOPE"
    TOPOSLPX.units = ""
    TOPOSLPX.stagger = ""
    TOPOSLPX.coordinates = "XLONG XLAT"
    
    ### Elevation y slope
    TOPOSLPY = ncfile.createVariable('TOPOSLPY','f4',('Time','south_north','west_east'))
    TOPOSLPY.FieldType = 104
    TOPOSLPY.MemoryOrder = "XY "
    TOPOSLPY.description = "ELEVATION Y SLOPE"
    TOPOSLPY.units = ""
    TOPOSLPY.stagger = ""
    TOPOSLPY.coordinates = "XLONG XLAT"
    
    ### Annual max vegetation fraction
    SHDMAX = ncfile.createVariable('SHDMAX','f4',('Time','south_north','west_east'))
    SHDMAX.FieldType = 104
    SHDMAX.MemoryOrder = "XY "
    SHDMAX.description = "ANNUAL MAX VEG FRACTION"
    SHDMAX.units = ""
    SHDMAX.stagger = ""
    SHDMAX.coordinates = "XLONG XLAT"
    
    ### Annual min vegetation fraction
    SHDMIN = ncfile.createVariable('SHDMIN','f4',('Time','south_north','west_east'))
    SHDMIN.FieldType = 104
    SHDMIN.MemoryOrder = "XY "
    SHDMIN.description = "ANNUAL MIN VEG FRACTION"
    SHDMIN.units = ""
    SHDMIN.stagger = ""
    SHDMIN.coordinates = "XLONG XLAT"
    
    ### Annual max snow albedo in fraction
    SNOALB = ncfile.createVariable('SNOALB','f4',('Time','south_north','west_east'))
    SNOALB.FieldType = 104
    SNOALB.MemoryOrder = "XY "
    SNOALB.description = "ANNUAL MAX SNOW ALBEDO IN FRACTION"
    SNOALB.units = ""
    SNOALB.stagger = ""
    SNOALB.coordinates = "XLONG XLAT"
    
    ### Land use fraction by category
    LANDUSEF = ncfile.createVariable('LANDUSEF','f4',('Time','land_cat_stag','south_north','west_east'))
    LANDUSEF.FieldType = 104
    LANDUSEF.MemoryOrder = "XYZ"
    LANDUSEF.description = "LANDUSE FRACTION BY CATEGORY"
    LANDUSEF.units = ""
    LANDUSEF.stagger = "Z"
    LANDUSEF.coordinates = "XLONG XLAT"
    
    ### Soil cat fraction (top)
    SOILCTOP = ncfile.createVariable('SOILCTOP','f4',('Time','soil_cat_stag','south_north','west_east'))
    SOILCTOP.FieldType = 104
    SOILCTOP.MemoryOrder = "XYZ"
    SOILCTOP.description = "SOIL CAT FRACTION (TOP)"
    SOILCTOP.units = ""
    SOILCTOP.stagger = "Z"
    SOILCTOP.coordinates = "XLONG XLAT"
    
    ### Soil cat fraction (bottom)
    SOILCBOT = ncfile.createVariable('SOILCBOT','f4',('Time','soil_cat_stag','south_north','west_east'))
    SOILCBOT.FieldType = 104
    SOILCBOT.MemoryOrder = "XYZ"
    SOILCBOT.description = "SOIL CAT FRACTION (BOTOM)"
    SOILCBOT.units = ""
    SOILCBOT.stagger = "Z"
    SOILCBOT.coordinates = "XLONG XLAT"
    
    ### Soil temperature
    TSLB = ncfile.createVariable('TSLB','f4',('Time','soil_layers_stag','south_north','west_east'))
    TSLB.FieldType = 104
    TSLB.MemoryOrder = "XYZ"
    TSLB.description = "SOIL TEMPERATURE"
    TSLB.units = "K"
    TSLB.stagger = "Z"
    TSLB.coordinates = "XLONG XLAT"
    
    ### Soil moisture
    SMOIS = ncfile.createVariable('SMOIS','f4',('Time','soil_layers_stag','south_north','west_east'))
    SMOIS.FieldType = 104
    SMOIS.MemoryOrder = "XYZ"
    SMOIS.description = "SOIL MOISTURE"
    SMOIS.units = "m3 m-3"
    SMOIS.stagger = "Z"
    SMOIS.coordinates = "XLONG XLAT"
    
    ### Soil liquid water
    SH2O = ncfile.createVariable('SH2O','f4',('Time','soil_layers_stag','south_north','west_east'))
    SH2O.FieldType = 104
    SH2O.MemoryOrder = "XYZ"
    SH2O.description = "SOIL LIQUID WATER"
    SH2O.units = "m3 m-3"
    SH2O.stagger = "Z"
    SH2O.coordinates = "XLONG XLAT"
    
    ### Relative soil moisture
    SMCREL = ncfile.createVariable('SMCREL','f4',('Time','soil_layers_stag','south_north','west_east'))
    SMCREL.FieldType = 104
    SMCREL.MemoryOrder = "XYZ"
    SMCREL.description = "RELATIVE SOIL MOISTURE"
    SMCREL.units = ""
    SMCREL.stagger = "Z"
    SMCREL.coordinates = "XLONG XLAT"
    
    ### Sea ice flag
    SEAICE = ncfile.createVariable('SEAICE','f4',('Time','south_north','west_east'))
    SEAICE.FieldType = 104
    SEAICE.MemoryOrder = "XY "
    SEAICE.description = "SEA ICE FLAG"
    SEAICE.units = ""
    SEAICE.stagger= ""
    SEAICE.coordinates = "XLONG XLAT"
    
    ### Dominant vegetation category
    IVGTYP = ncfile.createVariable('IVGTYP','i4',('Time','south_north','west_east'))
    IVGTYP.FieldType = 106
    IVGTYP.MemoryOrder = "XY "
    IVGTYP.description = "DOMINANT VEGETATION CATEGORY"
    IVGTYP.units = ""
    IVGTYP.stagger = ""
    IVGTYP.coordinates = "XLONG XLAT"
    
    ### Dominant soil category
    ISLTYP = ncfile.createVariable('ISLTYP','i4',('Time','south_north','west_east'))
    ISLTYP.FieldType = 106
    ISLTYP.MemoryOrder = "XY "
    ISLTYP.description = "DOMINANT SOIL CATEGORY"
    ISLTYP.units = ""
    ISLTYP.stagger = ""
    ISLTYP.coordinates = "XLONG XLAT"
    
    ### Vegetation fraction
    VEGFRA = ncfile.createVariable('VEGFRA','f4',('Time','south_north','west_east'))
    VEGFRA.FieldType = 104
    VEGFRA.MemoryOrder = "XY "
    VEGFRA.description = "VEGETATION FRACTION"
    VEGFRA.units = ""
    VEGFRA.stagger = ""
    VEGFRA.coordinates = "XLONG XLAT"
    
    ### Snow water equivalent
    SNOW = ncfile.createVariable('SNOW','f4',('Time','south_north','west_east'))
    SNOW.FieldType = 104
    SNOW.MemoryOrder = "XY "
    SNOW.description = "SNOW WATER EQUIVALENT"
    SNOW.units = "kg m-2"
    SNOW.stagger = ""
    SNOW.coordinates = "XLONG XLAT"
    
    ### Physical snow depth
    SNOWH = ncfile.createVariable('SNOWH','f4',('Time','south_north','west_east'))
    SNOWH.FieldType = 104
    SNOWH.MemoryOrder = "XY "
    SNOWH.description = "PHYSICAL SNOW DEPTH"
    SNOWH.units = "m"
    SNOWH.stagger = ""
    SNOWH.coordinates = "XLONG XLAT"
    
    ### Canopy water
    CANWAT = ncfile.createVariable('CANWAT','f4',('Time','south_north','west_east'))
    CANWAT.FieldType = 104
    CANWAT.MemoryOrder = "XY "
    CANWAT.description = "CANOPY WATER"
    CANWAT.units = "kg m-2"
    CANWAT.stagger = ""
    CANWAT.coordinates = "XLONG XLAT"
    
    ### Logical physical snow depth
    FNDSNOWH = ncfile.createVariable('FNDSNOWH','i4','Time')
    FNDSNOWH.FieldType = 106
    FNDSNOWH.MemoryOrder = "0  "
    FNDSNOWH.description = "SNOWH_LOGICAL"
    FNDSNOWH.units = "-"
    FNDSNOWH.stagger = ""
    
    ### Logical soil water
    FNDSOILW = ncfile.createVariable('FNDSOILW','i4','Time')
    FNDSOILW.FieldType = 106
    FNDSOILW.MemoryOrder = "0  "
    FNDSOILW.description = "SOILW_LOGICAL"
    FNDSOILW.units = "-"
    FNDSOILW.stagger = ""
    
    ### Logical albedo
    FNDALBSI = ncfile.createVariable('FNDALBSI','i4','Time')
    FNDALBSI.FieldType = 106
    FNDALBSI.MemoryOrder = "0  "
    FNDALBSI.description = "ALBSI_LOGICAL"
    FNDALBSI.units = "-"
    FNDALBSI.stagger = ""
    
    ### Logical snow
    FNDSNOWSI = ncfile.createVariable('FNDSNOWSI','i4','Time')
    FNDSNOWSI.FieldType = 106
    FNDSNOWSI.MemoryOrder = "0  "
    FNDSNOWSI.description = "SNOWSI_LOGICAL"
    FNDSNOWSI.units = "-"
    FNDSNOWSI.stagger = ""
    
    ### Logical ice depth
    FNDICEDEPTH = ncfile.createVariable('FNDICEDEPTH','i4','Time')
    FNDICEDEPTH.FieldType = 106
    FNDICEDEPTH.MemoryOrder = "0  "
    FNDICEDEPTH.description = "ICEDEPTH_LOGICAL"
    FNDICEDEPTH.units = "-"
    FNDICEDEPTH.stagger = ""
    
    ### Lake depth
    LAKE_DEPTH = ncfile.createVariable('LAKE_DEPTH','f4',('Time','south_north','west_east'))
    LAKE_DEPTH.FieldType = 104
    LAKE_DEPTH.MemoryOrder = "XY "
    LAKE_DEPTH.description = "lake depth"
    LAKE_DEPTH.units = "m"
    LAKE_DEPTH.stagger = ""
    LAKE_DEPTH.coordinates = "XLONG XLAT"
    
    ### Sea surface zonal currents
    UOCE = ncfile.createVariable('UOCE','f4',('Time','south_north','west_east'))
    UOCE.FieldType = 104
    UOCE.MemoryOrder = "XY "
    UOCE.description = "SEA SURFACE ZONAL CURRENTS"
    UOCE.units = "m s-1"
    UOCE.stagger = ""
    UOCE.coordinates = "XLONG XLAT"
    
    ### Sea surface meridional currents
    VOCE = ncfile.createVariable('VOCE','f4',('Time','south_north','west_east'))
    VOCE.FieldType = 104
    VOCE.MemoryOrder = "XY "
    VOCE.description = "SEA SURFACE MERIDIONAL CURRENTS"
    VOCE.units = "m s-1"
    VOCE.stagger = ""
    VOCE.coordinates = "XLONG XLAT"
    
    ### Urban fraction
    FRC_URB2D = ncfile.createVariable('FRC_URB2D','f4',('Time','south_north','west_east'))
    FRC_URB2D.FieldType = 104
    FRC_URB2D.MemoryOrder = "XY "
    FRC_URB2D.description = "URBAN FRACTION"
    FRC_URB2D.units = "dimensionless"
    FRC_URB2D.stagger = ""
    FRC_URB2D.coordinates = "XLONG XLAT"
    
    ### Leaf area index
    LAI = ncfile.createVariable('LAI','f4',('Time','south_north','west_east'))
    LAI.FieldType = 104
    LAI.MemoryOrder = "XY "
    LAI.description = "LEAF AREA INDEX"
    LAI.units = "m-2/m-2"
    LAI.stagger = ""
    LAI.coordinates = "XLONG XLAT"
    
    ### Orographic variance
    VAR = ncfile.createVariable('VAR','f4',('Time','south_north','west_east'))
    VAR.FieldType = 104
    VAR.MemoryOrder = "XY "
    VAR.description = "OROGRAPHIC VARIANCE"
    VAR.units = ""
    VAR.stagger = ""
    VAR.coordinates = "XLONG XLAT"
    
    ### Map scale factor on mass grid
    MAPFAC_M = ncfile.createVariable('MAPFAC_M','f4',('Time','south_north','west_east'))
    MAPFAC_M.FieldType = 104
    MAPFAC_M.MemoryOrder = "XY "
    MAPFAC_M.description = "Map scale factor on mass grid"
    MAPFAC_M.units = ""
    MAPFAC_M.stagger = ""
    MAPFAC_M.coordinates = "XLONG XLAT"
    
    ### Map scale factor on u-grid
    MAPFAC_U = ncfile.createVariable('MAPFAC_U','f4',('Time','south_north','west_east_stag'))
    MAPFAC_U.FieldType = 104
    MAPFAC_U.MemoryOrder = "XY "
    MAPFAC_U.description = "Map scale factor on u-grid"
    MAPFAC_U.units = ""
    MAPFAC_U.stagger = "X"
    MAPFAC_U.coordinates = "XLONG_U XLAT_U"
    
    ### Map scale factor on v-grid
    MAPFAC_V = ncfile.createVariable('MAPFAC_V','f4',('Time','south_north_stag','west_east'))
    MAPFAC_V.FieldType = 104
    MAPFAC_V.MemoryOrder = "XY "
    MAPFAC_V.description = "Map scale factor on v-grid"
    MAPFAC_V.units = ""
    MAPFAC_V.stagger = "Y"
    MAPFAC_V.coordinates = "XLONG_V XLAT_V"
    
    ### Map scale factor on mass grid, x direction
    MAPFAC_MX = ncfile.createVariable('MAPFAC_MX','f4',('Time','south_north','west_east'))
    MAPFAC_MX.FieldType = 104
    MAPFAC_MX.MemoryOrder = "XY "
    MAPFAC_MX.description = "Map scale factor on mass grid, x direction"
    MAPFAC_MX.units = ""
    MAPFAC_MX.stagger = ""
    MAPFAC_MX.coordinates = "XLONG XLAT"
    
    ### Map scale factor on mass grid, y direction
    MAPFAC_MY = ncfile.createVariable('MAPFAC_MY','f4',('Time','south_north','west_east'))
    MAPFAC_MY.FieldType = 104
    MAPFAC_MY.MemoryOrder = "XY "
    MAPFAC_MY.description = "Map scale factor on mass grid, y direction"
    MAPFAC_MY.units = ""
    MAPFAC_MY.stagger = ""
    MAPFAC_MY.coordinates = "XLONG XLAT"
    
    ### Map scale factor on u-grid, x direction
    MAPFAC_UX = ncfile.createVariable('MAPFAC_UX','f4',('Time','south_north','west_east_stag'))
    MAPFAC_UX.FieldType = 104
    MAPFAC_UX.MemoryOrder = "XY "
    MAPFAC_UX.description = "Map scale factor on u-grid, x direction"
    MAPFAC_UX.units = ""
    MAPFAC_UX.stagger = "X"
    MAPFAC_UX.coordinates = "XLONG_U XLAT_U"
    
    ### Map scale factor on u-grid, y direction
    MAPFAC_UY = ncfile.createVariable('MAPFAC_UY','f4',('Time','south_north','west_east_stag'))
    MAPFAC_UY.FieldType = 104
    MAPFAC_UY.MemoryOrder = "XY "
    MAPFAC_UY.description = "Map scale factor on u-grid, x direction"
    MAPFAC_UY.units = ""
    MAPFAC_UY.stagger = "X"
    MAPFAC_UY.coordinates = "XLONG_U XLAT_U"
    
    ### Map scale factor on v-grid, x direction
    MAPFAC_VX = ncfile.createVariable('MAPFAC_VX','f4',('Time','south_north_stag','west_east'))
    MAPFAC_VX.FieldType = 104
    MAPFAC_VX.MemoryOrder = "XY "
    MAPFAC_VX.description = "Map scale factor on v-grid, x direction"
    MAPFAC_VX.units = ""
    MAPFAC_VX.stagger = "Y"
    MAPFAC_VX.coordinates = "XLONG_V XLAT_V"
    
    ### Inverse map scale factor on v-grid, x direction
    MF_VX_INV = ncfile.createVariable('MF_VX_INV','f4',('Time','south_north_stag','west_east'))
    MF_VX_INV.FieldType = 104
    MF_VX_INV.MemoryOrder = "XY "
    MF_VX_INV.description = "Inverse map scale factor on v-grid, x direction"
    MF_VX_INV.units = ""
    MF_VX_INV.stagger = "Y"
    MF_VX_INV.coordinates = "XLONG_V XLAT_V"
    
    ### Map scale factor on v-grid, x direction
    MAPFAC_VY = ncfile.createVariable('MAPFAC_VY','f4',('Time','south_north_stag','west_east'))
    MAPFAC_VY.FieldType = 104
    MAPFAC_VY.MemoryOrder = "XY "
    MAPFAC_VY.description = "Map scale factor on v-grid, x direction"
    MAPFAC_VY.units = ""
    MAPFAC_VY.stagger = "Y"
    MAPFAC_VY.coordinates = "XLONG_V XLAT_V"
    
    ### Coriolis sine latitude term
    F = ncfile.createVariable('F','f4',('Time','south_north','west_east'))
    F.FieldType = 104
    F.MemoryOrder = "XY "
    F.description = "Coriolis sine latitude term"
    F.units = "s-1"
    F.stagger = ""
    F.coordinates = "XLONG XLAT"
    
    ### Coriolis cosine latitude term
    E = ncfile.createVariable('E','f4',('Time','south_north','west_east'))
    E.FieldType = 104
    E.MemoryOrder = "XY "
    E.description = "Coriolis cosine latitude term"
    E.units = "s-1"
    E.stagger = ""
    E.coordinates = "XLONG XLAT"
    
    ### Local sine of map rotation
    SINALPHA = ncfile.createVariable('SINALPHA','f4',('Time','south_north','west_east'))
    SINALPHA.FieldType = 104
    SINALPHA.MemoryOrder = "XY "
    SINALPHA.description = "Local sine of map rotation"
    SINALPHA.units = ""
    SINALPHA.stagger = ""
    SINALPHA.coordinates = "XLONG XLAT"
    
    ### Local cosine of map rotation
    COSALPHA = ncfile.createVariable('COSALPHA','f4',('Time','south_north','west_east'))
    COSALPHA.FieldType = 104
    COSALPHA.MemoryOrder = "XY "
    COSALPHA.description = "Local cosine of map rotation"
    COSALPHA.units = ""
    COSALPHA.stagger = ""
    COSALPHA.coordinates = "XLONG XLAT"
    
    ### Terrain Height
    HGT = ncfile.createVariable('HGT','f4',('Time','south_north','west_east'))
    HGT.FieldType = 104
    HGT.MemoryOrder = "XY "
    HGT.description = "Terrain Height"
    HGT.units = "m"
    HGT.stagger = ""
    HGT.coordinates = "XLONG XLAT"
    
    ### Surface skin temperature
    TSK = ncfile.createVariable('TSK','f4',('Time','south_north','west_east'))
    TSK.FieldType = 104
    TSK.MemoryOrder = "XY "
    TSK.description = "SURFACE SKIN TEMPERATURE"
    TSK.units = "K"
    TSK.stagger = ""
    TSK.coordinates = "XLONG XLAT"
    
    ### Base state x wind in idealized cases
    U_BASE = ncfile.createVariable('U_BASE','f4',('Time','bottom_top'))
    U_BASE.FieldType = 104
    U_BASE.MemoryOrder = "Z  "
    U_BASE.description = "BASE STATE X WIND IN IDEALIZED CASES"
    U_BASE.units = ""
    U_BASE.stagger = ""
    
    ### Base state y wind in idealized cases
    V_BASE = ncfile.createVariable('V_BASE','f4',('Time','bottom_top'))
    V_BASE.FieldType = 104
    V_BASE.MemoryOrder = "Z  "
    V_BASE.description = "BASE STATE Y WIND IN IDEALIZED CASES"
    V_BASE.units = ""
    V_BASE.stagger = ""
    
    ### Base state qv in idealized cases
    QV_BASE = ncfile.createVariable('QV_BASE','f4',('Time','bottom_top'))
    QV_BASE.FieldType = 104
    QV_BASE.MemoryOrder = "Z  "
    QV_BASE.description = "BASE STATE QV IN IDEALIZED CASES"
    QV_BASE.units = ""
    QV_BASE.stagger = ""
    
    ### Base state height in idealized cases
    Z_BASE = ncfile.createVariable('Z_BASE','f4',('Time','bottom_top'))
    Z_BASE.FieldType = 104
    Z_BASE.MemoryOrder = "Z  "
    Z_BASE.description = "BASE STATE HEIGHT IN IDEALIZED CASES"
    Z_BASE.units = ""
    Z_BASE.stagger = ""
    
    ### Frame x wind
    U_FRAME = ncfile.createVariable('U_FRAME','f4','Time')
    U_FRAME.FieldType = 104
    U_FRAME.MemoryOrder = "0  "
    U_FRAME.description = "FRAME X WIND"
    U_FRAME.units = "m s-1"
    U_FRAME.stagger = ""
    
    ### Frame y wind
    V_FRAME = ncfile.createVariable('V_FRAME','f4','Time')
    V_FRAME.FieldType = 104
    V_FRAME.MemoryOrder = "0  "
    V_FRAME.description = "FRAME Y WIND"
    V_FRAME.units = "m s-1"
    V_FRAME.stagger = ""
    
    ### Pressure top of the model
    P_TOP = ncfile.createVariable('P_TOP','f4','Time')
    P_TOP.FieldType = 104
    P_TOP.MemoryOrder = "0  "
    P_TOP.description = "PRESSURE TOP OF THE MODEL"
    P_TOP.units = "Pa"
    P_TOP.stagger = ""
    
    ### Base state temperature
    T00 = ncfile.createVariable('T00','f4','Time')
    T00.FieldType = 104
    T00.MemoryOrder = "0  "
    T00.description = "BASE STATE TEMPERATURE"
    T00.units = "K"
    T00.stagger = ""
    
    ### Base state pressure
    P00 = ncfile.createVariable('P00','f4','Time')
    P00.FieldType = 104
    P00.MemoryOrder = "0  "
    P00.description = "BASE STATE PRESSURE"
    P00.units = "Pa"
    P00.stagger = ""
    
    ### Base state lapse rate
    TLP = ncfile.createVariable('TLP','f4','Time')
    TLP.FieldType = 104
    TLP.MemoryOrder = "0  "
    TLP.description = "BASE STATE LAPSE RATE"
    TLP.units = ""
    TLP.stagger = ""
    
    ### Temperature at which the base theta turns constant
    TISO = ncfile.createVariable('TISO','f4','Time')
    TISO.FieldType = 104
    TISO.MemoryOrder = "0  "
    TISO.description = "TEMP AT WHICH THE BASE T TURNS CONST"
    TISO.units = "K"
    TISO.stagger = ""
    
    ### Base state lapse rate in stratosphere
    TLP_STRAT = ncfile.createVariable('TLP_STRAT','f4','Time')
    TLP_STRAT.FieldType = 104
    TLP_STRAT.MemoryOrder = "0  "
    TLP_STRAT.description = "BASE STATE LAPSE RATE (DT/D(LN(P))) IN STRATOSPHERE"
    TLP_STRAT.units = "K"
    TLP_STRAT.stagger = ""
    
    ### Base state pressure at bottom of stratosphere
    P_STRAT = ncfile.createVariable('P_STRAT','f4','Time')
    P_STRAT.FieldType = 104
    P_STRAT.MemoryOrder = "0  "
    P_STRAT.description = "BASE STATE PRESSURE AT BOTTOM OF STRATOSPHERE"
    P_STRAT.units = "Pa"
    P_STRAT.stagger = ""
    
    ### Latitude, x-staggered
    XLAT_U = ncfile.createVariable('XLAT_U','f4',('Time','south_north','west_east_stag'))
    XLAT_U.FieldType = 104
    XLAT_U.MemoryOrder = "XY "
    XLAT_U.description = "LATITUDE, SOUTH IS NEGATIVE"
    XLAT_U.units = "degree_north"
    XLAT_U.stagger = "X"
    XLAT_U.coordinates = "XLONG_U XLAT_U"
    
    ### Longitude, x-staggered
    XLONG_U = ncfile.createVariable('XLONG_U','f4',('Time','south_north','west_east_stag'))
    XLONG_U.FieldType = 104
    XLONG_U.MemoryOrder = "XY "
    XLONG_U.description = "LONGITUDE, WEST IS NEGATIVE"
    XLONG_U.units = "degree_east"
    XLONG_U.stagger = "X"
    XLONG_U.coordinates = "XLONG_U XLAT_U"
    
    ### Latitude, y-staggered
    XLAT_V = ncfile.createVariable('XLAT_V','f4',('Time','south_north_stag','west_east'))
    XLAT_V.FieldType = 104
    XLAT_V.MemoryOrder = "XY "
    XLAT_V.description = "LATITUDE, SOUTH IS NEGATIVE"
    XLAT_V.units = "degree_north"
    XLAT_V.stagger = "Y"
    XLAT_V.coordinates = "XLONG_V XLAT_V"
    
    ### Longitude, y-staggered
    XLONG_V = ncfile.createVariable('XLONG_V','f4',('Time','south_north_stag','west_east'))
    XLONG_V.FieldType = 104
    XLONG_V.MemoryOrder = "XY "
    XLONG_V.description = "LONGITUDE, WEST IS NEGATIVE"
    XLONG_V.units = "degree_east"
    XLONG_V.stagger = "Y"
    XLONG_V.coordinates = "XLONG_V XLAT_V"
    
    ### Computational grid latitude
    CLAT = ncfile.createVariable('CLAT','f4',('Time','south_north','west_east'))
    CLAT.FieldType = 104
    CLAT.MemoryOrder = "XY "
    CLAT.description = "COMPUTATIONAL GRID LATITUDE, SOUTH IS NEGATIVE"
    CLAT.units = "degree_north"
    CLAT.stagger = ""
    CLAT.coordinates = "XLONG XLAT"
    
    ### Background albedo
    ALBBCK = ncfile.createVariable('ALBBCK','f4',('Time','south_north','west_east'))
    ALBBCK.FieldType = 104
    ALBBCK.MemoryOrder = "XY "
    ALBBCK.description = "BACKGROUND ALBEDO"
    ALBBCK.units = ""
    ALBBCK.stagger = ""
    ALBBCK.coordinates = "XLONG XLAT"
    
    ### Soil temperature at lower boundary
    TMN = ncfile.createVariable('TMN','f4',('Time','south_north','west_east'))
    TMN.FieldType = 104
    TMN.MemoryOrder = "XY "
    TMN.description = "SOIL TEMPERATURE AT LOWER BOUNDARY"
    TMN.units = "K"
    TMN.stagger = ""
    TMN.coordinates = "XLONG XLAT"
    
    ### Land mask, x direction
    XLAND = ncfile.createVariable('XLAND','f4',('Time','south_north','west_east'))
    XLAND.FieldType = 104
    XLAND.MemoryOrder = "XY "
    XLAND.description = "LAND MASK (1 FOR LAND, 2 FOR WATER)"
    XLAND.units = ""
    XLAND.stagger = ""
    XLAND.coordinates = "XLONG XLAT"
    
    ### Coupling mask
    CPLMASK = ncfile.createVariable('CPLMASK','f4',('Time','num_ext_model_couple_dom_stag','south_north','west_east'))
    CPLMASK.FieldType = 104
    CPLMASK.MemoryOrder = "XYZ "
    CPLMASK.description = "COUPLING MASK (0:VALUE FROM SST UPDATE; 1:VALUE FROM COUPLED OCEAN), vertical di"
    CPLMASK.units = ""
    CPLMASK.stagger= "Z"
    CPLMASK.coordinates = "XLONG XLAT"
    
    ### Flag indicating snow coverage
    SNOWC = ncfile.createVariable('SNOWC','f4',('Time','south_north','west_east'))
    SNOWC.FieldType = 104
    SNOWC.MemoryOrder = "XY "
    SNOWC.description = "FLAG INDICATING SNOW COVERAGE (1 FOR SNOW COVER)"
    SNOWC.units = ""
    SNOWC.stagger= ""
    SNOWC.coordinates = "XLONG XLAT"
    
    ### Fraction of frozen precipitation
    SR = ncfile.createVariable('SR','f4',('Time','south_north','west_east'))
    SR.FieldType = 104
    SR.MemoryOrder = "XY "
    SR.description = "fraction of frozen precipitation"
    SR.units = "-"
    SR.stagger= ""
    SR.coordinates = "XLONG XLAT"
    
    ### Save topography
    SAVE_TOPO_FROM_REAL = ncfile.createVariable('SAVE_TOPO_FROM_REAL','i4','Time')
    SAVE_TOPO_FROM_REAL.FieldType = 106
    SAVE_TOPO_FROM_REAL.MemoryOrder = "0  "
    SAVE_TOPO_FROM_REAL.description = "1=original topo from real/0=topo modified by WRF"
    SAVE_TOPO_FROM_REAL.units = "flag"
    SAVE_TOPO_FROM_REAL.stagger = ""
    
    ### Lake flag
    LAKEFLAG = ncfile.createVariable('LAKEFLAG','i4','Time')
    LAKEFLAG.FieldType = 106
    LAKEFLAG.MemoryOrder = "0  "
    LAKEFLAG.description = "Flag for lake in the global attributes for metgrid data"
    LAKEFLAG.units = "-"
    LAKEFLAG.stagger = ""
    
    ### Lake depth flag
    LAKE_DEPTH_FLAG = ncfile.createVariable('LAKE_DEPTH_FLAG','i4','Time')
    LAKE_DEPTH_FLAG.FieldType = 106
    LAKE_DEPTH_FLAG.MemoryOrder = "0  "
    LAKE_DEPTH_FLAG.description = "Flag for lake depth in the global attributes for metgrid data"
    LAKE_DEPTH_FLAG.units = "-"
    LAKE_DEPTH_FLAG.stagger = ""
    
    ### Land mask
    LANDMASK = ncfile.createVariable('LANDMASK','f4',('Time','south_north','west_east'))
    LANDMASK.FieldType = 104
    LANDMASK.MemoryOrder = "XY "
    LANDMASK.description = "LAND MASK (1 FOR LAND, 0 FOR WATER)"
    LANDMASK.units = ""
    LANDMASK.stagger = ""
    LANDMASK.coordinates = "XLONG XLAT"
    
    ### Lake mask
    LAKEMASK = ncfile.createVariable('LAKEMASK','f4',('Time','south_north','west_east'))
    LAKEMASK.FieldType = 104
    LAKEMASK.MemoryOrder = "XY "
    LAKEMASK.description = "LAKE MASK (1 FOR LAKE, 0 FOR NON-LAKE)"
    LAKEMASK.units = ""
    LAKEMASK.stagger = ""
    LAKEMASK.coordinates = "XLONG XLAT"
    
    ### Sea surface temperature
    SST = ncfile.createVariable('SST','f4',('Time','south_north','west_east'))
    SST.FieldType = 104
    SST.MemoryOrder = "XY "
    SST.description = "SEA SURFACE TEMPERATURE"
    SST.units = "K"
    SST.stagger = ""
    SST.coordinates = "XLONG XLAT"
    
    time_char = list(time_str)
    Times[0,:] = time_char[:]
    XLAT[0,:,:] = np.zeros((ny,nx))
    XLONG[0,:,:] = np.zeros((ny,nx))
    LU_INDEX[0,:,:] = np.ones((ny,nx))*num_soil_cat
    ZNU[0,:] = znu[:] 
    ZNW[0,:] = znw[:] 
    ZS[0,:] = np.zeros(num_soil_layers)
    DZS[0,:] = np.zeros(num_soil_layers)
    VAR_SSO[0,:,:] = np.zeros((ny,nx))
    U[0,:,:,:] = u[:,:,:]
    V[0,:,:,:] = v[:,:,:]
    W[0,:,:,:] = np.zeros((nz+1,ny,nx))
    PH[0,:,:,:] = ph[:,:,:]
    PHB[0,:,:,:] = phb[:,:,:]
    T[0,:,:,:] = t[:,:,:]
    T_INIT[0,:,:,:] = t_init[:,:,:]
    MU[0,:,:] = mu[:,:]
    MUB[0,:,:] = mub[:,:]
    P[0,:,:,:] = p[:,:,:]
    PB[0,:,:,:] = pb[:,:,:]
    FNM[0,:] = fnm[:]
    FNP[0,:] = fnp[:] 
    RDNW[0,:] = rdnw[:]
    RDN[0,:] = rdn[:] 
    DNW[0,:] = dnw[:] 
    DN[0,:] = dn[:]
    T_BASE[0,:] = t_base[:]
    CFN[0] = cfn
    CFN1[0] = cfn1
    STEP_NUMBER[0] = 0
    P_HYD[0,:,:,:] = pb[:,:,:]
    Q2[0,:,:] = np.zeros((ny,nx))
    T2[0,:,:] = np.zeros((ny,nx))
    TH2[0,:,:] = np.zeros((ny,nx))
    PSFC[0,:,:] = np.zeros((ny,nx))
    U10[0,:,:] = np.zeros((ny,nx))
    V10[0,:,:] = np.zeros((ny,nx))
    RDX[0] = rdx
    RDY[0] = rdy
    DTS[0] = 0
    DTSEPS[0] = 0
    RESM[0] = 0
    ZETATOP[0] = 0
    CF1[0] = cf1
    CF2[0] = cf2
    CF3[0] = cf3
    QVAPOR[0,:,:,:] = moist[:,:,:]
    QCLOUD[0,:,:,:] = np.zeros((nz,ny,nx))
    QRAIN[0,:,:,:] = np.zeros((nz,ny,nx))
    QICE[0,:,:,:] = np.zeros((nz,ny,nx))
    QSNOW[0,:,:,:] = np.zeros((nz,ny,nx))
    QGRAUP[0,:,:,:] = np.zeros((nz,ny,nx))
    QHAIL[0,:,:,:] = np.zeros((nz,ny,nx))
    QNCLOUD[0,:,:,:] = np.zeros((nz,ny,nx))
    QNRAIN[0,:,:,:] = np.zeros((nz,ny,nx))
    QNICE[0,:,:,:] = np.zeros((nz,ny,nx))
    QNSNOW[0,:,:,:] = np.zeros((nz,ny,nx))
    QNGRAUPEL[0,:,:,:] = np.zeros((nz,ny,nx))
    QNHAIL[0,:,:,:] = np.zeros((nz,ny,nx))
    QVGRAUPEL[0,:,:,:] = np.zeros((nz,ny,nx))
    FCX[0,:] = np.zeros(spec_bdy_width)
    GCX[0,:] = np.zeros(spec_bdy_width)
    DTBC[0] = 0
    TOPOSLPX[0,:,:] = np.zeros((ny,nx))
    TOPOSLPY[0,:,:] = np.zeros((ny,nx))
    SHDMAX[0,:,:] = np.zeros((ny,nx))
    SHDMIN[0,:,:] = np.zeros((ny,nx))
    SNOALB[0,:,:] = np.zeros((ny,nx))
    LANDUSEF[0,:,:,:] = np.zeros((num_land_cat,ny,nx))
    SOILCTOP[0,:,:,:] = np.zeros((num_soil_cat,ny,nx))
    SOILCBOT[0,:,:,:] = np.zeros((num_soil_cat,ny,nx))
    TSLB[0,:,:,:] = np.zeros((num_soil_layers,ny,nx))
    SMOIS[0,:,:,:] = np.zeros((num_soil_layers,ny,nx))
    SH2O[0,:,:,:] = np.zeros((num_soil_layers,ny,nx))
    SMCREL[0,:,:,:] = np.zeros((num_soil_layers,ny,nx))
    SEAICE[0,:,:] = np.zeros((ny,nx))
    IVGTYP[0,:,:] = np.zeros((ny,nx))
    ISLTYP[0,:,:] = np.zeros((ny,nx))
    VEGFRA[0,:,:] = np.zeros((ny,nx))
    SNOW[0,:,:] = np.zeros((ny,nx))
    SNOWH[0,:,:] = np.zeros((ny,nx))
    CANWAT[0,:,:] = np.zeros((ny,nx))
    FNDSNOWH[0] = 0
    FNDSOILW[0] = 0
    FNDALBSI[0] = 0
    FNDSNOWSI[0] = 0
    FNDICEDEPTH[0] = 0
    LAKE_DEPTH[0,:,:] = np.zeros((ny,nx))
    UOCE[0,:,:] = np.zeros((ny,nx))
    VOCE[0,:,:] = np.zeros((ny,nx))
    FRC_URB2D[0,:,:] = np.zeros((ny,nx))
    LAI[0,:,:] = np.zeros((ny,nx))
    VAR[0,:,:] = np.zeros((ny,nx))
    MAPFAC_M[0,:,:] = np.ones((ny,nx))
    MAPFAC_U[0,:,:] = np.ones((ny,nx+1))
    MAPFAC_V[0,:,:] = np.ones((ny+1,nx))
    MAPFAC_MX[0,:,:] = np.ones((ny,nx))
    MAPFAC_MY[0,:,:] = np.ones((ny,nx))
    MAPFAC_UX[0,:,:] = np.ones((ny,nx+1))
    MAPFAC_UY[0,:,:] = np.ones((ny,nx+1))
    MAPFAC_VX[0,:,:] = np.ones((ny+1,nx))
    MF_VX_INV[0,:,:] = np.ones((ny+1,nx))
    MAPFAC_VY[0,:,:] = np.ones((ny+1,nx))
    F[0,:,:] = np.ones((ny,nx))*FF
    E[0,:,:] = np.zeros((ny,nx))
    SINALPHA[0,:,:] = np.zeros((ny,nx))
    COSALPHA[0,:,:] = np.ones((ny,nx))
    HGT[0,:,:] = np.zeros((ny,nx))
    TSK[0,:,:] = tsk[:,:]
    U_BASE[0,:] = u_base[:]
    V_BASE[0,:] = v_base[:]
    QV_BASE[0,:] = qv_base[:]
    Z_BASE[0,:] = np.zeros(nz)
    U_FRAME[0] = 0
    V_FRAME[0] = 0
    P_TOP[0] = int(ptop)
    T00[0] = 0
    P00[0] = 0
    TLP[0] = 0
    TISO[0] = 0
    TLP_STRAT[0] = 0
    P_STRAT[0] = 0
    XLAT_U[0,:,:] = np.zeros((ny,nx+1))
    XLONG_U[0,:,:] = np.zeros((ny,nx+1))
    XLAT_V[0,:,:] = np.zeros((ny+1,nx))
    XLONG_V[0,:,:] = np.zeros((ny+1,nx))
    CLAT[0,:,:] = np.zeros((ny,nx))
    ALBBCK[0,:,:] = np.zeros((ny,nx))
    TMN[0,:,:] = tmn[:,:]
    XLAND[0,:,:] = np.ones((ny,nx))*2
    CPLMASK[0,:,:,:] = np.zeros((1,ny,nx))
    SNOWC[0,:,:] = np.zeros((ny,nx))
    SR[0,:,:] = np.zeros((ny,nx))
    SAVE_TOPO_FROM_REAL[0] = 0
    LAKEFLAG[0] = 0
    LAKE_DEPTH_FLAG[0] = 0
    LANDMASK[0,:,:] = np.zeros((ny,nx))
    LAKEMASK[0,:,:] = np.zeros((ny,nx))
    SST[0,:,:] = np.ones((ny,nx))*RD

    return ncfile
