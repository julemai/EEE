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
Template files for Efficient Elementary Effects sensitivity analysis of CEQUEAU

History
-------
Written,  JM, Jun 2019
"""

EXECUTION_XML = """<?xml version="1.0"?>
<!-- ==========================================                                                             -->
<!-- Written using Efficient Elementary Effects (EEE) Toolbox                                               -->
<!-- https://github.com/julemai/EEE                                                                         -->
<!-- ==========================================                                                             -->
<!-- M Cuntz & J Mai et al. (2015).                                                                         -->
<!-- Computationally inexpensive identification of noninformative model parameters by sequential screening. -->
<!-- Water Resources Research, 51, 6417-6441.                                                               -->
<!-- https://doi.org/10.1002/2015WR016907.                                                                  -->
<!-- ==========================================                                                             -->
<root>
  <!-- simulation start date -->
  <dateDebut>{setup[start_day]}</dateDebut>
  <!-- simulation end date -->
  <dateFin>{setup[end_day]}</dateFin>
  <fichierParamSimulation>{setup[tmp_folder]}/parametres.xml</fichierParamSimulation>
  <fichierBassinVersant>{setup[tmp_folder]}/bassinVersant.xml</fichierBassinVersant>
  <fichierMeteo>{setup[tmp_folder]}/data_obs/meteoCorrected.nc</fichierMeteo>
  <!-- initial states initialization -->
  <fichierEtats>{setup[tmp_folder]}/data_obs/resultats_1953-01-01_2019-09-20_meteoCorrected.nc</fichierEtats>
  <!-- results of this run -->
  <fichierResultats>{setup[tmp_folder]}/output/resultats.nc</fichierResultats>
</root>

"""

PARAMETRES_XML = """<?xml version="1.0"?>
<!-- ==========================================                                                             -->
<!-- Written using Efficient Elementary Effects (EEE) Toolbox                                               -->
<!-- https://github.com/julemai/EEE                                                                         -->
<!-- ==========================================                                                             -->
<!-- M Cuntz & J Mai et al. (2015).                                                                         -->
<!-- Computationally inexpensive identification of noninformative model parameters by sequential screening. -->
<!-- Water Resources Research, 51, 6417-6441.                                                               -->
<!-- https://doi.org/10.1002/2015WR016907.                                                                  -->
<!-- ==========================================                                                             -->
<root>
  <option>
    <ipassim>24</ipassim>
    <moduleFonte>2</moduleFonte>
    <moduleEvapo>5</moduleEvapo>
    <calculQualite>0</calculQualite>
    <logNeigeAjustee>1</logNeigeAjustee>
    <sortiesMatricielles>1</sortiesMatricielles>
  </option>
  <sol>
    <cin_s>0.02750	0.02750	0.02750	0.02750	0.02750	0.04000	0.04000	0.06700	0.06700	0.06700	0.005000</cin_s>
    <cvnb_s>0.005000	0.005000	0.005000	0.005000	0.005000	0.006840	0.006840	0.01880	0.01880	0.01880	0.005000</cvnb_s>
    <cvnh_s>0	0	0	0	0	0	0	0	0	0	0</cvnh_s>
    <cvsi_s>5	5	5	5	5	1.0430	1.0430	3.5430	3.5430	3.5430	3</cvsi_s>
    <hinf_s>82	82	82	82	82	105	105	100	100	100	10</hinf_s>
    <hint_s>72	72	72	72	72	100.0	100.0	24.010	24.010	24.010	72.030</hint_s>
    <hnap_s>0 0 0 0 0 0 0 0 0 0 0</hnap_s>
    <hpot_s>0 0 0 0 0 0 0 0 0 0 0</hpot_s>
    <hsol_s>250	250	250	250	250	396	396	337.60	337.60	337.60	420</hsol_s>
    <hrimp_s>0 0 0 0 0 0 0 0 0 0 0</hrimp_s>
    <cvmar_s>0 0 0 0 0 0 0 0 0 0 0</cvmar_s>
    <cvsb_s>0.005000	0.005000	0.005000	0.005000	0.005000	0.008680	0.008680	0.012840	0.012840	0.012840	0.01670</cvsb_s>
    <xinfma_s>48	48	48	48	48	5	5	5	5	5	36.0200000000000</xinfma_s>
    <hmar_s>0 0 0 0 0 0 0 0 0 0 0</hmar_s>
    <xla_s>4915	4915	4915	4915	4915	4915	4915	4915	4915	4915	4915</xla_s>
    <trifonte_s>0	0	0	0	0	0	0	0	0	0	0</trifonte_s>
    <vidintmax_s>15 15 15 15 15 15 15 15 15 15 15</vidintmax_s>
  </sol>
  <solInitial>
    <hsini>80</hsini>
    <hnini>40</hnini>
    <hmini>0</hmini>
    <q0>326</q0>
  </solInitial>
  <transfert>
    <exxkt>0</exxkt>
    <zn>0</zn>
    <HUDebit>10  7	10	10	90	7	90	90	0	7	0	0	0	7	0	0	0	7	0	0	0	7	0	0	0	7	0	0</HUDebit>
    <HUProd>30	52	10	0	60	180	30	10	20	265	60	30	0	262	20	60	0	237	0	20	0	205	0	0	0	0	0	0</HUProd>
  </transfert>
  <lac>0</lac>
  <surface>0</surface>
  <interpolation>
    <type>3</type>
    <coep>0.75</coep>
    <coet>-2.5</coet>
  </interpolation>
  <fonte>
    <cemaneige>
      <Kf>5</Kf>
      <Tf>-0.5000   -0.5000   -0.5000   -0.5000   -0.5000   -0.5000   -0.5000   -0.5000   -0.5000   -0.5000   -0.5000</Tf>
      <CTg>0.43</CTg>
      <Gseuil>300</Gseuil>
      <Vmin>0.05</Vmin>
      <strne_s>-0.882</strne_s>
      <PourPluieNeige>1</PourPluieNeige>
    </cemaneige>
  </fonte>
  <evapo>
    <oudin>
      <evnap_s>0 0 0 0 0 0 0 0 0 0 0</evnap_s>
      <K1>1.3730</K1>
      <K2>6.4800    6.4800    6.4800    6.4800    6.4800   -2.0000   -2.0000    2.0000    6.4800    2.0000    3.0000</K2>
    </oudin>
  </evapo>
  <relevesNeige/> 
  <ajustements>
    <!-- par_v1 -->
	<niveauEauSol isProp = "false">{par[v1]}</niveauEauSol>
	<!-- par_v2 -->
	<niveauEauNappe isProp = "false">{par[v2]}</niveauEauNappe>
	<!-- par_v3 -->	
	<G isProp = "false">{par[v3]}</G>
	<!-- par_v4 -->
	<eTg isProp = "false">{par[v4]}</eTg>
    <!-- par_v5 NOT USED -->
    <indexTempNeige isProp = "false">{par[v5]}</indexTempNeige>
	<!-- par_v6 -->
	<precipitation_a>{par[v6]}</precipitation_a>
	<!-- par_v7 -->
	<precipitation_b>{par[v7]}</precipitation_b>
	<!-- par_v8 -->
	<temperature_a>{par[v8]}</temperature_a>
	<!-- par_v9 -->
	<temperature_b>{par[v9]}</temperature_b>
  </ajustements>
</root>

"""

BASSINVERSANT_XML = """<?xml version="1.0"?>
<!-- ==========================================                                                             -->
<!-- Written using Efficient Elementary Effects (EEE) Toolbox                                               -->
<!-- https://github.com/julemai/EEE                                                                         -->
<!-- ==========================================                                                             -->
<!-- M Cuntz & J Mai et al. (2015).                                                                         -->
<!-- Computationally inexpensive identification of noninformative model parameters by sequential screening. -->
<!-- Water Resources Research, 51, 6417-6441.                                                               -->
<!-- https://doi.org/10.1002/2015WR016907.                                                                  -->
<!-- ==========================================                                                             -->
<root>
  <nomBassinVersant>SystemeLSJ</nomBassinVersant>
  <carreauxEntiers>	
	<carreauEntier>
		<id>PD01</id>
		<ij>2841 2842 2843 2844 2941 2942 2943</ij>
		<latitude>51.09601695399794</latitude>
		<longitude>-71.22386756047717</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>100</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>0</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>PD02</id>
		<ij>2536 2537 2538 2539 2636 2637 2638 2639 2736 2737 2738 2839 2840</ij>
		<latitude>50.67689204899216</latitude>
		<longitude>-71.51494581329033</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>100</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>0</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>PD03</id>
		<ij>2735 2833 2834 2835 2836 2837 2838 2933 2934 2935 2936 2937 2938 2939 2940 3034 3035</ij>
		<latitude>50.51586633532165</latitude>
		<longitude>-71.16564011123055</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>100</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>0</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>PD05</id>
		<ij>2433 2434 2532 2533 2534 2535 2632 2633 2634 2635 2733 2734</ij>
		<latitude>50.2979928133249</latitude>
		<longitude>-71.59763828790162</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>100</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>0</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>PD15</id>
		<ij>2730 2731 2732 2830</ij>
		<latitude>50.05616102800558</latitude>
		<longitude>-71.34009179492774</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>0</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>100</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>PERIB01</id>
		<ij>2444 2445 2543 2544 2545 2546 2643 2644 2645 2646 2647 2648 2649 2650 2742 2743 2744 2745 2746 2747 2748 2749 2750 2751 2841 2842 2843 2844 2845 2846 2847 2848 2849 2850 2851 2852 2853 2942 2943 2944 2945 2946 2947 2948 2949 2950 2951 2952 2953 2954 3041 3042 3043 3044 3045 3046 3047 3048 3049 3050 3051 3052 3053 3054 3148 3149 3150 3151 3152 3153 3154 3155 3249 3250 3251 3252 3253 3254 3255 3353 3355</ij>
		<latitude>51.60984551420005</latitude>
		<longitude>-71.2080078005961</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>100</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>0</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>PERIB02</id>
		<ij>2242 2243 2341 2342 2343 2344 2440 2441 2442 2443 2540 2541 2542 2640 2641 2642 2739 2740 2741 2839 2840</ij>
		<latitude>50.99082207133508</latitude>
		<longitude>-71.72824337761911</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>100</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>0</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>LM01</id>
		<ij>2941 3041 3042 3142 3242 3243 3244 3245 3343 3344 3345 3443 3444</ij>
		<latitude>51.16848166358147</latitude>
		<longitude>-70.72704231258305</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>100</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>0</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>LM03</id>
		<ij>2936 2937 2938 2939 3034 3035 3036 3037 3039 3040 3135 3136 3138 3139 3140 3141 3235 3236 3239 3240 3241 3336 3337 3338 3339 3340 3341 3342 3440 3441 3442</ij>
		<latitude>50.74413438660478</latitude>
		<longitude>-70.77818833804743</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>100</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>0</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>LM14</id>
		<ij>3038 3137 3237 3238</ij>
		<latitude>50.67280033655355</latitude>
		<longitude>-70.80445572778304</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>0</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>100</pctLacExutoire>
	</carreauEntier>	
	<carreauEntier>
		<id>MBLANC01</id>
		<ij>3043 3044 3045 3046 3047 3143 3144 3145 3146 3147 3148 3149 3242 3243 3244 3245 3246 3247 3248 3249 3345 3346 3347 3348 3349 3447 3448</ij>
		<latitude>51.4343835904832</latitude>
		<longitude>-70.76881034842944</longitude>
		<pctLacRiviere>0</pctLacRiviere>
		<pctForet>100</pctForet>
		<pctMarais>0</pctMarais>
		<pctLacExutoire>0</pctLacExutoire>
	</carreauEntier>	
  </carreauxEntiers>
  <carreauxPartiels>
	    <carreauPartiel>
		<id>PD</id>
		<nom>Passes Dangereuses</nom>
		<superficie>4026</superficie>
		<superficieReservoir>345</superficieReservoir>
		<idCPAval>0</idCPAval>
		<idCPsAmont>2 0 0 0 0 0 0</idCPsAmont>
		<idCE>1 2 3 4 5</idCE>
		<superficieCE>383 989 1407 902 345</superficieCE>
	</carreauPartiel>
    <carreauPartiel>
		<id>PERIB</id>
		<nom>Peribonka</nom>
		<superficie>7505</superficie>
		<superficieReservoir>0</superficieReservoir>
		<idCPAval>1</idCPAval>
		<idCPsAmont>0 0 0 0 0 0 0</idCPsAmont>
		<idCE>6 7 0 0 0</idCE>
		<superficieCE>6194 1311 0 0 0</superficieCE>
	</carreauPartiel>
    <carreauPartiel>
		<id>LM</id>
		<nom>Lac Manouane</nom>
		<superficie>2998</superficie>
		<superficieReservoir>432</superficieReservoir>
		<idCPAval>0</idCPAval>
		<idCPsAmont>4 0 0 0 0 0 0</idCPsAmont>
		<idCE>8 9 10 0 0</idCE>
		<superficieCE>631 1935 432 0 0</superficieCE>
	</carreauPartiel>
    <carreauPartiel>
		<id>MBLANC</id>
		<nom>Montagnes Blanches</nom>
		<superficie>1900</superficie>
		<superficieReservoir>0</superficieReservoir>
		<idCPAval>3</idCPAval>
		<idCPsAmont>0 0 0 0 0 0 0</idCPsAmont>
		<idCE>11 0 0 0 0</idCE>
		<superficieCE>1900 0 0 0 0</superficieCE>
	</carreauPartiel>
  </carreauxPartiels>
  <barrage/> 
</root>
"""

