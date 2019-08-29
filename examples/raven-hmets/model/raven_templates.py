#!/usr/bin/env python

# Copyright 2019 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of the EEE code library for "Computationally inexpensive identification
# of noninformative model parameters by sequential screening: Efficient Elementary Effects (EEE)".
#
# The EEE code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The MVA code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with The EEE code library.
# If not, see <https://github.com/julemai/EEE/blob/master/LICENSE>.
#
# If you use this method in a publication please cite:
#
#    M Cuntz & J Mai et al. (2015).
#    Computationally inexpensive identification of noninformative model parameters by sequential screening.
#    Water Resources Research, 51, 6417-6441.
#    https://doi.org/10.1002/2015WR016907.
#

from __future__ import print_function

"""
Template files for Efficient Elementary Effects sensitivity analysis of RAVEN

History
-------
Written,  JM, Jun 2019
"""

RVI = """
#########################################################################                                  
:FileType          rvi ASCII Raven rev217 (v2.9)                                                                             
:WrittenBy         James Craig & Juliane Mai                                                                       
:CreationDate      June 2019
#
# RAVEN run of Salmon River near Prince George using HMETS model setup                                                        
#------------------------------------------------------------------------
#
:StartDate               1989-01-01 00:00:00 # 1954-01-01 00:00:00 
:EndDate                 2010-12-31 00:00:00  
:EvaluationTime          1991-01-01 00:00:00                                                        
# :Duration              20819                                                                                   
:TimeStep                1.0                                                                                   
:Method                  ORDERED_SERIES 

:PotentialMeltMethod     POTMELT_HMETS
:RainSnowFraction        RAINSNOW_DATA
:SWRadiationMethod       SW_RAD_NONE         # no radiation is faster
:Evaporation             PET_DATA
:CatchmentRoute          ROUTE_DUMP
:Routing                 ROUTE_NONE 
:SoilModel               SOIL_TWO_LAYER

:Alias DELAYED_RUNOFF CONVOLUTION[1] 

:HydrologicProcesses
  :Precipitation   RAVEN_DEFAULT           ATMOS_PRECIP   MULTIPLE 
  :Infiltration    INF_HMETS               PONDED_WATER   MULTIPLE 
    :Overflow      OVERFLOW_RAVEN          SOIL[0]        DELAYED_RUNOFF
  :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[0]        SURFACE_WATER   # interflow, really
  :Percolation     PERC_LINEAR             SOIL[0]        SOIL[1]         # recharge
    :Overflow      OVERFLOW_RAVEN          SOIL[1]        DELAYED_RUNOFF
  :SoilEvaporation SOILEVAP_ALL            SOIL[0]        ATMOSPHERE      # AET
  :Convolve        CONVOL_GAMMA            CONVOLUTION[0] SURFACE_WATER   # 'surface runoff'
  :Convolve        CONVOL_GAMMA_2          DELAYED_RUNOFF SURFACE_WATER   # 'delayed runoff'
  :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[1]        SURFACE_WATER
  :SnowBalance     SNOBAL_HMETS            MULTIPLE       MULTIPLE
:EndHydrologicProcesses

#:CreateRVPTemplate

#---------------------------------------------------------
# Output Options
#
# :WriteForcingFunctions
# :WriteNetcdfFormat

# Accumulated Infiltration volume
:CustomOutput DAILY AVERAGE Between:PONDED_WATER.And.SOIL[0] BY_BASIN

:EvaluationMetrics NASH_SUTCLIFFE RMSE
:SilentMode
:DontWriteWatershedStorage
#

"""

RVP = """
#########################################################################                                  
:FileType          rvp ASCII Raven rev217 (v2.9)                                                                             
:WrittenBy         James Craig & Juliane Mai                                                                       
:CreationDate      June 2019
#
# RAVEN run of Salmon River near Prince George using HMETS model setup                                                               
#------------------------------------------------------------------------                                 
#

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x06" and "par_x10" and "par_x10" wouldn't be detectable)
#    para_sum_x05_x06 = {dpar[sum_x05_x06]} =  par_x05 + par_x06  = {par[x05]} + {par[x06]}
#    para_sum_x09_x10 = {dpar[sum_x09_x10]} =  par_x09 + par_x10  = {par[x09]} + {par[x10]}

#-----------------------------------------------------------------
# Soil Classes
#-----------------------------------------------------------------
:SoilClasses
  :Attributes,
  :Units,
  TOPSOIL,
  PHREATIC,    
:EndSoilClasses

#-----------------------------------------------------------------
# Land Use Classes
#-----------------------------------------------------------------
:LandUseClasses, 
  :Attributes,        IMPERM,    FOREST_COV, 
       :Units,          frac,          frac, 
       FOREST,           0.0,           1.0,    
:EndLandUseClasses

#-----------------------------------------------------------------
# Vegetation Classes
#-----------------------------------------------------------------
:VegetationClasses, 
  :Attributes,        MAX_HT,       MAX_LAI, MAX_LEAF_COND, 
       :Units,             m,          none,      mm_per_s, 
       FOREST,             4,             5,             5,     
:EndVegetationClasses

#-----------------------------------------------------------------
# Soil Profiles
#-----------------------------------------------------------------
:SoilProfiles
         LAKE, 0
         ROCK, 0
  DEFAULT_P, 2, TOPSOIL, {par[x20]}, PHREATIC, {par[x21]}, 
# DEFAULT_P, 2, TOPSOIL,      x(20), PHREATIC,      x(21), 
:EndSoilProfiles

#-----------------------------------------------------------------
# Global Parameters
#-----------------------------------------------------------------
:GlobalParameter         SNOW_SWI_MIN {par[x09]}           # x(9)    
:GlobalParameter         SNOW_SWI_MAX {dpar[sum_x09_x10]}  # x(9)+x(10) 
:GlobalParameter     SWI_REDUCT_COEFF {par[x11]}           # x(11)
:GlobalParameter             SNOW_SWI 0.05   #not sure why/if needed... 

#-----------------------------------------------------------------
# Soil Parameters
#-----------------------------------------------------------------
:SoilParameterList
  :Parameters,        POROSITY,      PERC_COEFF,  PET_CORRECTION, BASEFLOW_COEFF
       :Units,               -,             1/d,               -,            1/d
      TOPSOIL,             1.0,      {par[x17]},      {par[x15]},      {par[x18]} 
     PHREATIC,             1.0,             0.0,             0.0,      {par[x19]}    
 #    TOPSOIL,             1.0,           x(17),           x(15),          x(18)
 #   PHREATIC,             1.0,             0.0,             0.0,          x(19)
:EndSoilParameterList

#-----------------------------------------------------------------
# Land Use Parameters
#-----------------------------------------------------------------
:LandUseParameterList
  :Parameters, MIN_MELT_FACTOR,     MAX_MELT_FACTOR,    DD_MELT_TEMP,  DD_AGGRADATION, REFREEZE_FACTOR,    REFREEZE_EXP, DD_REFREEZE_TEMP, HMETS_RUNOFF_COEFF,
       :Units,          mm/d/C,              mm/d/C,               C,            1/mm,          mm/d/C,               -,                C,                  -,
    [DEFAULT],      {par[x05]}, {dpar[sum_x05_x06]},      {par[x07]},      {par[x08]},      {par[x13]},      {par[x14]},       {par[x12]},         {par[x16]},
#                         x(5),           x(5)+x(6),            x(7),            x(8),           x(13),           x(14),            x(12),              x(16),        
:EndLandUseParameterList
:LandUseParameterList
  :Parameters,   GAMMA_SHAPE,     GAMMA_SCALE,    GAMMA_SHAPE2,    GAMMA_SCALE2, 
       :Units,             -,               -,               -,               -, 
    [DEFAULT],    {par[x01]},      {par[x02]},      {par[x03]},      {par[x04]}, 
    #                   x(1),            x(2),            x(3),            x(4), 
:EndLandUseParameterList

#-----------------------------------------------------------------
# Vegetation Parameters
#-----------------------------------------------------------------
:VegetationParameterList
  :Parameters,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT, 
       :Units,               -,               -,
    [DEFAULT],             0.0,             0.0,
:EndVegetationParameterList

"""

RVC = """
#########################################################################                                  
:FileType          rvc ASCII Raven rev217 (v2.9)                                                                             
:WrittenBy         James Craig & Juliane Mai                                                                       
:CreationDate      June 2019
#
# RAVEN run of Salmon River near Prince George using HMETS model setup                                                           
#------------------------------------------------------------------------                                 
#

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x20" and "par_x21" wouldn't be detectable)
#    para_half_x20 = para_x20 * 1000. / 2. = {par[x20]} / 2. [m] = {dpar[half_x20]} [mm]
#    para_half_x21 = para_x21 * 1000. / 2. = {par[x21]} / 2. [m] = {dpar[half_x21]} [mm]

# initialize to 1/2 full
#:UniformInitialConditions SOIL[0] {dpar[half_x20]} # x(20)*1000/2 [mm]         
#:UniformInitialConditions SOIL[1] {dpar[half_x21]} # x(21)*1000/2 [mm]         

:HRUStateVariableTable (formerly :IntialConditionsTable)
   :Attributes SOIL[0] SOIL[1]
   :Units mm mm
   1 {dpar[half_x20]} {dpar[half_x21]}
:EndHRUStateVariableTable

"""

RVT = """
#########################################################################                                  
:FileType          rvt ASCII Raven rev217 (v2.9)                                                                             
:WrittenBy         James Craig & Juliane Mai                                                                       
:CreationDate      June 2019
#
# RAVEN run of Salmon River near Prince George using HMETS model setup                                                          
#------------------------------------------------------------------------

# meteorological forcings
:Gauge
  :Latitude    54.09639
  :Longitude -122.67972
  :Elevation  606.0
  :RedirectToFile data_obs/Salmon-River-Near-Prince-George_meteo_daily.rvt
:EndGauge

# observed streamflow
:RedirectToFile data_obs/Salmon-River-Near-Prince-George_Qobs_daily.rvt

"""

RVH = """
#########################################################################                                  
:FileType          rvh ASCII Raven rev217 (v2.9)                                                                             
:WrittenBy         James Craig & Juliane Mai                                                                       
:CreationDate      June 2019
#
# RAVEN run of Salmon River near Prince George using HMETS model setup                                                            
#------------------------------------------------------------------------                            
#                                                                                                           
#                                                                                                               
:SubBasins                                                                                                              
        :Attributes     NAME    DOWNSTREAM_ID   PROFILE   REACH_LENGTH    GAUGED                                                          
        :Units          none    none            none      km              none                                                                                                    
        1,             hmets,   -1,             NONE,     _AUTO,          1
:EndSubBasins                                                                                                                           
                                                                                                                
:HRUs                                                                                                           
        :Attributes     AREA    ELEVATION  LATITUDE    LONGITUDE  BASIN_ID  LAND_USE_CLASS  VEG_CLASS SOIL_PROFILE AQUIFER_PROFILE TERRAIN_CLASS    SLOPE   ASPECT  
        :Units           km2            m       deg          deg      none            none       none         none            none          none      deg      deg     
                 1    4230.0,       606.0, 54.09639,  -122.67972,         1          FOREST     FOREST    DEFAULT_P         [NONE]        [NONE]      0.0        0
:EndHRUs

"""

