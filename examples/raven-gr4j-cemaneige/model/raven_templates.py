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
:FileType          rvi ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Emulation of GR4J simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#
#:RunName              run1
:StartDate             1989-01-01 00:00:00 # 1954-01-01 00:00:00
:EndDate               2010-12-31 00:00:00
:EvaluationTime        1991-01-01 00:00:00
# :Duration            20819
:TimeStep              1.0
:Method                ORDERED_SERIES

:SoilModel             SOIL_MULTILAYER  4
:Routing               ROUTE_NONE
:CatchmentRoute        ROUTE_DUMP
:Evaporation           PET_DATA
:RainSnowFraction      RAINSNOW_DINGMAN
:PotentialMeltMethod   POTMELT_DEGREE_DAY
:OroTempCorrect        OROCORR_SIMPLELAPSE
:OroPrecipCorrect      OROCORR_SIMPLELAPSE

#------------------------------------------------------------------------
# Soil Layer Alias Definitions
#
:Alias PRODUCT_STORE      SOIL[0]
:Alias ROUTING_STORE      SOIL[1]
:Alias TEMP_STORE         SOIL[2]
:Alias GW_STORE           SOIL[3]

#------------------------------------------------------------------------
# Hydrologic process order for GR4J Emulation
#
:HydrologicProcesses
 :Precipitation            PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE
 :SnowTempEvolve           SNOTEMP_NEWTONS    SNOW_TEMP
 :SnowBalance              SNOBAL_CEMA_NIEGE  SNOW            PONDED_WATER
 :OpenWaterEvaporation     OPEN_WATER_EVAP    PONDED_WATER    ATMOSPHERE     			 # Pn
 :Infiltration             INF_GR4J           PONDED_WATER    MULTIPLE       			 # Ps-
 :SoilEvaporation          SOILEVAP_GR4J      PRODUCT_STORE   ATMOSPHERE     			 # Es
 :Percolation              PERC_GR4J          PRODUCT_STORE   TEMP_STORE     			 # Perc
 :Flush                    RAVEN_DEFAULT      SURFACE_WATER   TEMP_STORE     			 # Pn-Ps
 :Split                    RAVEN_DEFAULT      TEMP_STORE      CONVOLUTION[0] CONVOLUTION[1] 0.9  # Split Pr
 :Convolve                 CONVOL_GR4J_1      CONVOLUTION[0]  ROUTING_STORE  			 # Q9
 :Convolve                 CONVOL_GR4J_2      CONVOLUTION[1]  TEMP_STORE     			 # Q1
 :Percolation              PERC_GR4JEXCH      ROUTING_STORE   GW_STORE       			 # F(x1)
 :Percolation              PERC_GR4JEXCH2     TEMP_STORE      GW_STORE       			 # F(x1)
 :Flush                    RAVEN_DEFAULT      TEMP_STORE      SURFACE_WATER  			 # Qd
 :Baseflow                 BASE_GR4J          ROUTING_STORE   SURFACE_WATER  			 # Qr
:EndHydrologicProcesses
#------------------------------------------------------------------------

#---------------------------------------------------------
# Output Options
:EvaluationMetrics NASH_SUTCLIFFE LOG_NASH KLING_GUPTA RMSE PCT_BIAS ABSERR ABSMAX PDIFF TMVOL RCOEF NSC RSR R2 PERSINDEX
#:SuppressOutput
:SilentMode
:DontWriteWatershedStorage

"""

RVP = """
#########################################################################
:FileType          rvp ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Emulation of GR4J simulation of Salmon River near Prince George
#------------------------------------------------------------------------

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x06" wouldn't be detectable)
#    para_1_minus_x06 = {dpar[par_1_minus_x06]} =  1.0 - par_x06 = 1.0 - {par[x06]}


# -Global snow parameters-------------------------------------
:RainSnowTransition {par[x08]} {par[x09]}             # 0.0 1.0 :: base temp, range around base temp where rain<-->snow, [-3,3], [0.5,4]
:AirSnowCoeff       {dpar[par_1_minus_x06]}  	      # [1/d] = 1.0 - CEMANEIGE_X2 = 1.0 - x6
:AvgAnnualSnow      {par[x05]}          	          # [mm]  =       CEMANEIGE_X1 =       x5

# -Orographic Corrections-------------------------------------
:PrecipitationLapseRate 0.0004
:AdiabaticLapseRate 0.0065

# - Soil classes ---------------------------------------------
:SoilClasses
  :Attributes
  :Units
   SOIL_PROD
   SOIL_ROUT
   SOIL_TEMP
   SOIL_GW
:EndSoilClasses
:SoilParameterList
 :Parameters, POROSITY ,      GR4J_X3,     GR4J_X2
 :Units     ,     none ,           mm,        mm/d
   [DEFAULT],      1.0 ,   {par[x03]},  {par[x02]}
:EndSoilParameterList

# ----Soil Profiles--------------------------------------------
#     name,#horizons,[soiltype,thickness]x[#horizons]
#     GR4J_X1 is thickness of first layer (SOIL_PROD), here 0.696
:SoilProfiles
  DEFAULT_P,      4, SOIL_PROD   , {par[x01]}, SOIL_ROUT  ,   0.300, SOIL_TEMP  ,   1.000, SOIL_GW  ,   1.000,
:EndSoilProfiles

# ----Vegetation Classes---------------------------------------
:VegetationClasses
   :Attributes,       MAX_HT,       MAX_LAI,      MAX_LEAF_COND
        :Units,            m,          none,           mm_per_s
       VEG_ALL,           0.0,          0.0,                0.0
:EndVegetationClasses

# --Land Use Classes------------------------------------------
:LandUseClasses
  :Attributes, IMPERM, FOREST_COV
  :Units     ,   frac,       frac
       LU_ALL,    0.0,        0.0
:EndLandUseClasses
:LandUseParameterList
 :Parameters,     GR4J_X4,         MELT_FACTOR
 :Units     ,           d,              mm/d/C
   [DEFAULT],  {par[x04]},          {par[x07]}              # p4=GR4J_X4 and p7=CEMANEIGE_X3 = [2,8]
:EndLandUseParameterList

"""

RVC = """
#########################################################################
:FileType          rvc ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Emulation of GR4J simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#
# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x01" wouldn't be detectable)
#    para_half_x01 = para_x01 * 1000. / 2. = {par[x01]} / 2. [m] = {dpar[half_x01]} [mm]

# initialize to 1/2 full
# GR4J_X1 * 1000. / 2.0
#:UniformInitialConditions SOIL[0] {dpar[half_x01]} # x(1)*1000/2 [mm]

:HRUStateVariableTable (formerly :IntialConditionsTable)
   :Attributes SOIL[0]          SOIL[1]
   :Units      mm               mm
   1           {dpar[half_x01]} 15.0
:EndHRUStateVariableTable

"""

RVT = """
#########################################################################
:FileType          rvt ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Emulation of GR4J simulation of Salmon River near Prince George
#------------------------------------------------------------------------

# meteorological forcings
:Gauge 08KC001
  :Latitude    54.4848
  :Longitude -123.3659
  :Elevation  843.0
  :RedirectToFile data_obs/Salmon-River-Near-Prince-George_meteo_daily.rvt
:EndGauge

# observed streamflow
:RedirectToFile data_obs/Salmon-River-Near-Prince-George_Qobs_daily.rvt

"""

RVH = """
#########################################################################
:FileType          rvh ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Emulation of GR4J simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#
:SubBasins
  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,GAUGED
  :Units     ,          none,          none,   none,          km,  none
            1,   gr4j-salmon,            -1,   NONE,       _AUTO,     1
:EndSubBasins
:HRUs
  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS,SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT
  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,        none,            none,          none,   deg,    deg
            1,4250.6,     843.0,  54.4848, -123.3659,        1,        LU_ALL,   VEG_ALL,   DEFAULT_P,          [NONE],        [NONE],     0,      0
:EndHRUs

"""
