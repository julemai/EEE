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
Template files for Efficient Elementary Effects sensitivity analysis of ROBIN (coupled with Raven)

History
-------
Adapted from JM (Jun 2019) by SL (Mar 2022)
"""

PAR = """##### crop paramter
59
Veg_Name                Forest                  Forest2
extcoef                 {par[x01]}              {par[x02]}              # 5.750153E-01  5.750153E-01
rue                     {par[x03]}              {par[x04]}              # 3.423747      4.712127
alpha_r		        0.4	        	0.4
base_temp		{par[x05]}		{par[x06]}		# 5.008859E+00	7.940524E+00
opt_temp		{par[x07]}		{par[x08]}		# 15		15
max_temp		{par[x09]}		{par[x10]}		# 40		40
declr_rue	        0.5	        	9.106525E-03
phi_r_min		{par[x11]}		{par[x12]}		# 0.65	        6.996634E-01
phi_r_max		{par[x13]}		{par[x14]}		# 0.65		6.612798E-01
alpha_leaf_stems2	0.1			{par[x15]}		# 0.1           3.054366E-01
alpha_leaf_stems20	0.2			{par[x16]}		# 0.2           2.096685E-01
lambda_leaf_min		{par[x17]}		{par[x18]}		# 8.455921E-02	8.455921E-02 
lambda_leaf_max		{par[x19]}		{par[x20]}		# 8.455921E-02	8.455921E-02
yr_lambdaleaf		{par[x21]}		{par[x22]}		# 35		30
lambda_r		{par[x23]}		{par[x24]}		# 1.020223E-06	1.020223E-06
lambda_nv_min		{par[x25]}		{par[x26]}		# 4.286613E-05 	4.715102E-04 
lambda_nv_max		{par[x27]}		{par[x28]}		# 4.286613E-05	8.967972E-04
yr_lambdanv		{par[x29]}		{par[x30]}		# 35           	1
aplha_nv_power		{par[x31]}		{par[x32]}		# 9.277507E-01	0.5
lambda_nv_leaf		{par[x33]}		{par[x34]}		# 0.0	        0.0
lambda_nv_root		{par[x35]}		{par[x36]}		# 1	      	1
lambda_nv_stems		{par[x37]}		{par[x38]}		# 1	        1 
w_stems_max_1000nv	{par[x39]}		{par[x40]}		# 230	        100
alpha_n_max		{par[x41]}		{par[x42]}		# 1.5         	2
kc_min			{par[x43]}		{par[x44]}		# 0.4	        0.4
kc_max	             	{dpar[sum_x43_x45]}	{dpar[sum_x44_x46]}	# 0.9686279	0.9686279 para_sum_x43_x45 = {dpar[sum_x43_x45]} = par_x43 + par_x45 = {par[x43]} + {par[x45]}; para_sum_x44_x46 = {dpar[sum_x44_x46]} = par_x44 + par_x46 = {par[x44]} + {par[x46]}
laimin			{par[x47]}		{par[x48]}		# 0.26943    	0.26943
laimax	             	{dpar[sum_x47_x49]}	{dpar[sum_x48_x50]} 	# 4.5	        4.5 para_sum_x47_x49 = {dpar[sum_x47_x49]} = par_x47 + par_x49 = {par[x47]} + {par[x49]}; para_sum_x48_x50 = {dpar[sum_x48_x50]} = par_x48 + par_x50 = {par[x48]} + {par[x50]}
virtual_lai          	5.821848E-01   		1.100414E-01
sla_min			{par[x51]}		{par[x52]}		# 10     	1.962361E+01
sla_max			{par[x53]}		{par[x54]}		# 24.72         5.466540E+00 
yr_mid_sla		{par[x55]}		{par[x56]}		# 35           	2
alpha_dbh_a		{par[x57]}		{par[x58]}		# 0.145         0.45
alpha_dbh_n		{par[x59]}		{par[x60]}		# 2.45	        2
alpha_ht_a		{par[x61]}		{par[x62]}		# 4.63        	4.63
alpha_ht_n		{par[x63]}		{par[x64]}		# 0.44       	0.44
alpha_vol_n1		{par[x65]}		{par[x66]}		# 1.555103E+00	1.555103E+00
alpha_vol_n2		{par[x67]}		{par[x68]}		# 2.787484E-01	2.787484E-01
alpha_vol_n3		{par[x69]}		{par[x70]}		# 1.795241E-01	1.795241E-01
alpha_vol_a		{par[x71]}		{par[x72]}		# 2.411041E-04	2.411041E-04
veg_in_delay	     	60	            	30
veg_de_delay	     	60	            	60
r_wleaf_max	        0.7453703      		0.6
min_stem_tree        	20             		10
root_depth_max       	2              		2
identinum	        2	            	2
cum_day_t_lg_tbase   	2              		2
cum_day_t_lt_tcold	3	            	3
max_wleaf_a          	0.1122544      		7.539280E-02
max_wleaf_n	        1.2603	        	1.5203
cold_temp	        1.375298E+01		1.375298E+01
lambda_leaf_max_w    	0.00161878	    	0.00161878
lambda_leaf_max_w_n  	4.710727       		4.710727
lambda_leaf_max_t    	0.002266099    		0.002266099
lambda_leaf_max_t_n  	1.610429       		1.610429    
day_length_leaf_off  	11	            	11
reproduct_rate       	6.675921e-05   		8.639454E-04
pec_full_occupied    	0.5            		0.7
alpha_et             	0.7233215      		0.7233215
############		
0		

"""

HRUCROP = """subid hruid rotationname temp_ave
2	11	R1	18
2	33	R1	18
2	34	R1	18
2	73	R1	18
2	74	R1	18
11	12	R1	18
11	35	R1	18
11	36	R1	18
11	75	R1	18
11	76	R1	18
16	13	R1	18
16	77	R1	18
16	78	R1	18
16	79	R1	18
29	1	R1	18
29	14	R1	18
29	37	R1	18
29	38	R1	18
29	39	R1	18
29	80	R1	18
36	2	R1	18
36	15	R1	18
36	40	R1	18
36	41	R1	18
36	81	R1	18
36	82	R1	18
38	3	R1	18
38	83	R1	18
38	84	R1	18
38	85	R1	18
38	86	R1	18
40	87	R1	18
40	88	R1	18
40	89	R1	18
40	90	R1	18
40	91	R1	18
41	92	R1	18
41	93	R1	18
41	94	R1	18
41	95	R1	18
50	16	R1	18
50	42	R1	18
50	43	R1	18
50	96	R1	18
53	10	R1	18
53	17	R1	18
53	44	R1	18
53	45	R1	18
59	4	R1	18
59	97	R1	18
59	98	R1	18
59	99	R1	18
59	100	R1	18
73	18	R1	18
73	46	R1	18
73	47	R1	18
73	101	R1	18
73	102	R1	18
93	19	R1	18
93	48	R1	18
93	49	R1	18
93	103	R1	18
103	5	R1	18
103	20	R1	18
103	50	R1	18
103	51	R1	18
106	104	R1	18
106	105	R1	18
106	106	R1	18
106	107	R1	18
106	108	R1	18
115	21	R1	18
115	52	R1	18
115	53	R1	18
116	22	R2	18
116	54	R2	18
116	109	R2	18
116	110	R2	18
119	23	R1	18
119	55	R1	18
119	111	R1	18
133	112	R1	18
133	113	R1	18
133	114	R1	18
133	115	R1	18
134	24	R1	18
134	56	R1	18
134	57	R1	18
134	116	R1	18
142	117	R1	18
142	118	R1	18
142	119	R1	18
142	120	R1	18
150	121	R1	18
150	122	R1	18
150	123	R1	18
150	124	R1	18
150	125	R1	18
163	6	R1	18
163	25	R1	18
163	58	R1	18
163	59	R1	18
163	60	R1	18
164	7	R1	18
164	126	R1	18
164	127	R1	18
164	128	R1	18
164	129	R1	18
166	130	R1	18
166	131	R1	18
166	132	R1	18
166	133	R1	18
167	134	R1	18
167	135	R1	18
167	136	R1	18
167	137	R1	18
170	138	R1	18
170	139	R1	18
170	140	R1	18
170	141	R1	18
170	142	R1	18
173	26	R1	18
173	61	R1	18
173	62	R1	18
173	143	R1	18
178	27	R1	18
178	63	R1	18
178	144	R1	18
181	145	R1	18
181	146	R1	18
181	147	R1	18
181	148	R1	18
195	28	R1	18
195	149	R1	18
195	150	R1	18
205	151	R1	18
205	152	R1	18
205	153	R1	18
221	8	R1	18
221	29	R1	18
221	64	R1	18
221	65	R1	18
221	154	R1	18
222	155	R1	18
222	156	R1	18
222	157	R1	18
222	158	R1	18
263	9	R1	18
263	30	R1	18
263	66	R1	18
263	67	R1	18
263	68	R1	18
263	159	R1	18
275	31	R1	18
275	69	R1	18
275	160	R1	18
275	161	R1	18
275	162	R1	18
276	163	R1	18
276	164	R1	18
276	165	R1	18
276	166	R1	18
276	167	R1	18
288	32	R1	18
288	70	R1	18
288	71	R1	18
288	72	R1	18
288	168	R1	18
288	169	R1	18
rotationname	cropname1	cropname2	cropname3	cropname4	cropname5	cropname6	cropname7	cropname8	cropname9	cropname10	cropname11	cropname12	cropname13	cropname14	cropname15	cropname16	cropname17	cropname18	cropname19	cropname20	cropname21	cropname22	cropname23	cropname24	cropname25	cropname26	cropname27	cropname28	cropname29	cropname30	cropname31	cropname32	cropname33	cropname34	cropname35	cropname36	cropname37	cropname38	cropname39	cropname40	cropname41	cropname42	cropname43	cropname44	cropname45	cropname46	cropname47	cropname48	cropname49	cropname50	cropname51	cropname52	cropname53	cropname54	cropname55	cropname56	cropname57	cropname58	cropname59	cropname60	cropname61	cropname62	cropname63	cropname64	cropname65	cropname66	cropname67	cropname68	cropname69	cropname70	cropname71	cropname72	cropname73	cropname74	cropname75	cropname76	cropname77	cropname78	cropname79	cropname80	cropname81	cropname82	cropname83	cropname84	cropname85	cropname86	cropname87	cropname88
R2	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2	Forest2
R1	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest	Forest
"""

INFO = """1981      # read (ifile,*) Begin_Year
1       #  read (ifile,*) Begin_juliday
12050       #read (ifile,*) ndays
38         #read (ifile,*) nsub
169         #read (ifile,*) nhru
2         #read (ifile,*) ncrop
88         # nseeding event
88         # nharvest event
2         # nrotation
88         # number of growing season
999         # max hru id
999         # max sub id
1  #  NoNutrient
1981 01 01 ## YEAR MON DAY
1  # 1 suppresses screen output (0 turns it back on)
###
./cropout/
0  # is the number of customized output options. 0 means no customized outputs.
TIMESERIES		# Output type
-1				# OutputVars -1 means all variable
27 63 144 		# HRU IDs
-1				# sub-basin IDs -1 means no sub-basin
-1   			# timestep -1 means all model time steps
TIMESERIES
-1
-1
178
-1 
TIMESERIES
-1
22 110 109
-1
-1   
TIMESERIES
-1
-1
116
-1 
"""

INIC = """#####	crop	Initial	condition
13			
CropName	Forest	        Forest2	
w_total	    305725.7438	    4	
w_leaf	    5399.823839	    1056
w_stems	    1.358499E+05    36114.64968
w_seed	    0	            1
w_root	    33000	        500	
ht	        15	            5	
rootdep	    5	            5	
dbh	        17.5	        20	
basa_area	22	            6.3	
volume	    176	            1	
yr_c	    182	            1	
vg_number	1.095562E+03	201	
stage    	1	            1	
"""

MGT = """####
88
mgtype	index	Forest	Forest2
seeding	1	31	31
seeding	2	31	31
seeding	3	30	30
seeding	4	30	30
seeding	5	30	30
seeding	6	31	31
seeding	7	30	30
seeding	8	30	30
seeding	9	30	30
seeding	10	31	31
seeding	11	30	30
seeding	12	30	30
seeding	13	30	30
seeding	14	31	31
seeding	15	30	30
seeding	16	30	30
seeding	17	30	30
seeding	18	31	31
seeding	19	30	30
seeding	20	30	30
seeding	21	30	30
seeding	22	31	31
seeding	23	30	30
seeding	24	30	30
seeding	25	30	30
seeding	26	31	31
seeding	27	30	30
seeding	28	30	30
seeding	29	30	30
seeding	30	31	31
seeding	31	30	30
seeding	32	30	30
seeding	33	30	30
seeding	34	31	31
seeding	35	30	30
seeding	36	30	30
seeding	37	30	30
seeding	38	31	31
seeding	39	30	30
seeding	40	30	30
seeding	41	30	30
seeding	42	31	31
seeding	43	30	30
seeding	44	30	30
harvest	1	335	335
harvest	2	336	336
harvest	3	335	335
harvest	4	335	335
harvest	5	335	335
harvest	6	336	336
harvest	7	335	335
harvest	8	335	335
harvest	9	335	335
harvest	10	336	336
harvest	11	335	335
harvest	12	335	335
harvest	13	335	335
harvest	14	336	336
harvest	15	335	335
harvest	16	335	335
harvest	17	335	335
harvest	18	336	336
harvest	19	335	335
harvest	20	335	335
harvest	21	335	335
harvest	22	336	336
harvest	23	335	335
harvest	24	335	335
harvest	25	335	335
harvest	26	336	336
harvest	27	335	335
harvest	28	335	335
harvest	29	335	335
harvest	30	336	336
harvest	31	335	335
harvest	32	335	335
harvest	33	335	335
harvest	34	336	336
harvest	35	335	335
harvest	36	335	335
harvest	37	335	335
harvest	38	336	336
harvest	39	335	335
harvest	40	335	335
harvest	41	335	335
harvest	42	336	336
harvest	43	335	335
harvest	44	335	335
"""

MODEL = """21		
VEG_NAME	                  Forest	         Forest2
LAI_MODEL	                  3PG_LAI	         3PG_LAI
PHENOLOGY_MODEL	              3PG_PNLG	         3PG_PNLG
HT_MODEL	                  Ht_3PG	         Ht_3PG
WATER_STRESS_MODEL	          AETPET	         AETPET
TEMP_STRESS_MODEL	          Temp_3PG	         Temp_3PG
VPD_STRESS_MODEL	          VPD_3PG	         VPD_3PG
TOTAL_STRESS_MASS_ACC_MODEL	  MA_3PG	         MA_3PG
OPTIMAL_GPP_MODEL	          Beer_LAW	         Beer_LAW
RESPIRATION_MODEL	          R_Const_C	         R_Const_C
MASS_DIS_ROOT_MODEL	          BioA_Rt_3PG	     BioA_Rt_3PG
MASS_DIS_LEAF_MODEL	          BioA_LF_3PG	     BioA_LF_3PG
MASS_DIS_SEED_MODEL	          BioA_SD_3PG	     BioA_SD_3PG
DBH_MODEL	                  DBH_3PG	         DBH_3PG
BASAL_AREA_MODEL	          BasA_3PG	         BasA_3PG
STAND_VOlUME_MODEL	          StVol_3PG	         StVol_3PG
MASS_LOSS_ROOT_MODEL	      WL_Rt_3PG	         WL_Rt_3PG
MASS_LOSS_LEAF_MODEL	      WL_LF_3PG	         WL_LF_3PG
MASS_LOSS_VEG_NUM_MODEL	      WL_NV_3PG	         WL_NV_3PG
ROOT_DISTRIBUTION_MODEL	      ROOT_SWAT	         ROOT_SWAT
REPRODUCTION_MODEL            ReProd_3PG         ReProd_3PG
"""

OBS = """#add observation file definiation in here
nummber_of_observation_files 7
# obsfile 1
./cropmodel/obs/LAI.txt
obs_state_variables LAI
level_of_observation Subbasin
number_of_obs_time_step 166
number_of_obs_levels 2
# obsfile 2
./cropmodel/obs/Stocking.txt
obs_state_variables Stocking
level_of_observation Subbasin
number_of_obs_time_step 10
number_of_obs_levels 2
# obsfile 3
./cropmodel/obs/stembio.txt
obs_state_variables W_Stem
level_of_observation Subbasin
number_of_obs_time_step 10
number_of_obs_levels 2
# obsfile 4
./cropmodel/obs/BasalArea.txt
obs_state_variables BasalArea
level_of_observation Subbasin
number_of_obs_time_step 10
number_of_obs_levels 2
# obsfile 5
./cropmodel/obs/Volume.txt
obs_state_variables Volume
level_of_observation Subbasin
number_of_obs_time_step 10
number_of_obs_levels 2
# obsfile 6
./cropmodel/obs/stembiotree.txt
obs_state_variables W_Stem_Tree
level_of_observation Subbasin
number_of_obs_time_step 10
number_of_obs_levels 2
# obsfile 7
./cropmodel/obs/leafbio.txt
obs_state_variables W_Leaf
level_of_observation Subbasin
number_of_obs_time_step 10
number_of_obs_levels 2
"""
