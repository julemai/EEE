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

RVC = """
:UniformInitialConditions SOIL[2] 300.0
"""

RVH = """
#----------------------------------------------
# This is a Raven HRU rvh input file generated
# by Routing toolbox
#----------------------------------------------
:SubBasins
  :Attributes   NAME  DOWNSTREAM_ID       PROFILE REACH_LENGTH  GAUGED
  :Units        none           none          none           km    none
  29     sub29     36     Chn_29     ZERO-     0
  36     sub36     40     Chn_36     ZERO-     0
  38     sub38     41     Chn_38     ZERO-     0
  59     sub59     73     Chn_59     ZERO-     0
  103     sub103     134     Chn_103     ZERO-     0
  163     sub163     167     Chn_163     ZERO-     0
  164     sub164     205     Chn_164     ZERO-     0
  221     sub221     222     Chn_221     ZERO-     0
  263     sub263     275     Chn_263     ZERO-     1
  53     sub53     93     Chn_53     ZERO-     0
  2     sub2     36     Chn_2     ZERO-     1
  11     sub11     36     Chn_11     ZERO-     1
  16     sub16     36     Chn_16     ZERO-     1
  50     sub50     59     Chn_50     ZERO-     1
  73     sub73     142     Chn_73     0.394658849     1
  93     sub93     163     Chn_93     ZERO-     1
  115     sub115     288     Chn_115     ZERO-     1
  116     sub116     115     Chn_116     ZERO-     1
  119     sub119     276     Chn_119     ZERO-     1
  134     sub134     133     Chn_134     0.382246247     1
  173     sub173     164     Chn_173     ZERO-     1
  178     sub178     288     Chn_178     ZERO-     1
  195     sub195     276     Chn_195     0.40112006     1
  275     sub275     288     Chn_275     1.2977472680000002     0
  288     sub288     -1     Chn_288     1.027288638     0
  40     sub40     38     Chn_40     0.312945921     0
  41     sub41     59     Chn_41     0.256715584     0
  106     sub106     150     Chn_106     ZERO-     1
  133     sub133     263     Chn_133     0.100518492     0
  142     sub142     163     Chn_142     0.8169321169999999     0
  150     sub150     263     Chn_150     ZERO-     1
  166     sub166     164     Chn_166     0.07929807700000001     0
  167     sub167     166     Chn_167     0.297416061     1
  170     sub170     164     Chn_170     ZERO-     1
  181     sub181     275     Chn_181     ZERO-     1
  205     sub205     221     Chn_205     0.112994205     0
  222     sub222     263     Chn_222     0.057152743     0
  276     sub276     288     Chn_276     0.5627813899999999     0
:EndSubBasins


:HRUs
  :Attributes AREA ELEVATION  LATITUDE  LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT
  :Units       km2         m       deg        deg       none            none       none           none             none            none     deg      deg
  1     0.06612592123449999     495.881755315     47.0665044989     -84.3934565628     29     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.279407910319     178.153687421     
  2     0.0632384416956     495.378877878     47.063162549     -84.3909881881     36     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.318732815545     155.052415911     
  3     0.00540873851459     493.441697593     47.0628314235     -84.3953242709     38     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.236296604018     137.727288722     
  4     0.0165607112489     478.127251493     47.0616177717     -84.3989762408     59     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.727084545903     218.629370071     
  5     0.0500344610312     409.011980149     47.0359913005     -84.4188262365     103     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.25749956271     160.818488243     
  6     0.19423060140499998     386.488645132     47.049435491     -84.3992257326     163     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.185689223303     161.688632409     
  7     0.0383486959141     373.586192548     47.0450548842     -84.4006012604     164     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.120158052676     140.508239934     
  8     0.21290737624199998     370.678483747     47.0423478928     -84.4081798677     221     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.217765933542     178.038377193     
  9     0.558974609372     368.297732233     47.0484414262     -84.4219227181     263     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.161634411713     210.327781114     
  10     0.00262962940252     480.95328144     47.0537329715     -84.3927402774     53     LU_Lake     Veg_Lake     LAKE      [NONE]     [NONE]     0.26953736582     171.331810868     
  11     0.0440990842447     563.841072022     47.0678287688     -84.386932814     2     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.6645066425     138.399886474     
  12     0.06803248293929999     558.86107469     47.0661808819     -84.3846828603     11     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.65756422261     150.054015949     
  13     0.0123965221921     564.913362929     47.0629423444     -84.3845280828     16     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.67920066101     77.39995271499998     
  14     0.0920174465682     532.276901112     47.0682430165     -84.3916991692     29     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.57715328414     170.776099723     
  15     0.098178464527     511.770511224     47.0617706772     -84.3892340042     36     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.26191424293     197.214372578     
  16     0.14863323322299998     553.040032358     47.0658378389     -84.4046060973     50     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.99344627408     206.587232407     
  17     0.0318920089734     492.418421083     47.0549619256     -84.393292573     53     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.89878575814     210.031888169     
  18     0.151938685949     486.13047495     47.0598786132     -84.3972402893     73     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.38864765181     178.351049957     
  19     0.0490025240737     463.506165205     47.0535576707     -84.3952159024     93     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.9556496592     164.218420456     
  20     0.06765721598819999     428.292678793     47.035478285     -84.4210000094     103     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.08666658141     249.633467194     
  21     0.008524027659969999     370.009877613     47.0632918347     -84.4290538598     115     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.96150615363     129.44927267     
  22     0.0266093481899     404.745195061     47.0643083658     -84.4265003892     116     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.80957482627     173.130632137     
  23     0.191687818096     470.087101365     47.0626245685     -84.4170230697     119     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.21522656331     154.086629966     
  24     0.15467958139000001     413.409153833     47.0389052333     -84.418572712     134     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.3143420759     199.741422181     
  25     0.370199720245     426.920655729     47.0511305141     -84.397245413     163     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.86469571199     156.150958087     
  26     0.103682790668     416.410802768     47.0418736689     -84.3948458643     173     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.34976773614     149.18415272     
  27     0.0356215042121     407.901996621     47.0630177944     -84.4241538557     178     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.88178101602     160.697397533     
  28     0.352809367591     480.164959818     47.06059469     -84.4101122756     195     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.02881245259     164.165926538     
  29     0.361257711371     401.077130329     47.0393229504     -84.4072955613     221     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.44863785145     139.094771869     
  30     0.659721787625     397.210130468     47.0471089748     -84.4249068382     263     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.6305142181     206.359425718     
  31     0.223027700449     386.955931606     47.0546447402     -84.4204010749     275     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.76832653486     223.362009063     
  32     0.0962351250438     376.291021859     47.0611578346     -84.4256574396     288     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.1152803759     149.454077272     
  33     0.0129212712892     547.797961778     47.0672064003     -84.3874489457     2     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.53807991124     160.511425035     
  34     0.0075004158269500005     553.014927455     47.0676992178     -84.3874474139     2     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.756315440428     159.297900133     
  35     0.025891732023199998     549.421098341     47.0649468634     -84.3847612011     11     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.69237581863     89.511509536     
  36     0.0127798236854     550.779592276     47.0660957937     -84.3851310883     11     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.845436811727     141.809645168     
  37     0.09398453497959999     521.127828922     47.0672333794     -84.3951907833     29     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.15404300139     248.461278536     
  38     0.0263602386248     506.466317026     47.0678255418     -84.3935025157     29     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.669628872839     172.37212743     
  39     0.02965106826     529.265931561     47.0683343641     -84.3912824094     29     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     9.15130297436     136.211176877     
  40     0.0593395119746     511.988942558     47.0634874046     -84.3884229753     36     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.55792748306     147.342408559     
  41     0.0215685463248     499.43889211     47.0629396131     -84.3902180224     36     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.703581099632     187.611184149     
  42     0.141787723099     535.297049185     47.0640714755     -84.4047372168     50     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.78392670061     282.4758062119     
  43     0.10957036382     546.943577266     47.0661568953     -84.4027638669     50     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.70852192183     156.616751382     
  44     0.00880040208372     491.825125327     47.0547697833     -84.3935275404     53     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.09422124159     239.49550770899998     
  45     0.0108441051371     486.345199025     47.0544649027     -84.3928897144     53     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.660024296315     182.724112999     
  46     0.0824313569693     494.125806623     47.058287435     -84.3968199948     73     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.63491672183     201.078733309     
  47     0.0657343320028     479.23661777     47.0593575608     -84.3958509635     73     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.728811924551     170.378196049     
  48     0.0404802339833     486.31960894     47.0523930337     -84.3917824554     93     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.66354301501     151.670257112     
  49     0.0102098906798     479.498417752     47.0533427674     -84.3927612044     93     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.760402325246     141.726328539     
  50     0.0532206800243     433.190187912     47.0363287665     -84.4217651166     103     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.04264149898     224.387089068     
  51     0.022184803156     416.065400744     47.0358718052     -84.4186551009     103     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.889868476383     184.759603355     
  52     0.00284081978016     346.805089227     47.0629029823     -84.4299368313     115     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.1902676915     68.49412966199998     
  53     0.00248326942042     381.750816268     47.0636794311     -84.4287621339     115     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.02250644294     176.705234701     
  54     0.0104461820371     405.971938227     47.064129975     -84.4260595983     116     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     5.39252229799     79.63822633299998     
  55     0.0484331030214     470.728312937     47.0625067983     -84.4169864858     119     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.04413470791     139.237886006     
  56     0.129219808707     417.203215337     47.0393960559     -84.419106808     134     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.06863202702     185.094797354     
  57     0.037345780462     406.493306308     47.0383626533     -84.417578742     134     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.834515920339     180.414799044     
  58     0.182684841554     415.480567883     47.0489848954     -84.3975728095     163     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.35824193742     217.810947862     
  59     0.07091762197550001     408.987033071     47.0505251004     -84.3992419482     163     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.927985799452     184.750241689     
  60     0.13411448748700003     445.723515919     47.051342293     -84.3953120308     163     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.31617911426     139.075249777     
  61     0.018962482669299997     407.727695647     47.0434549523     -84.3928333544     173     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.2754737001     135.589216263     
  62     0.0144375166291     395.622391586     47.0429549398     -84.3949167852     173     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.891788742908     156.886594655     
  63     0.0128449774457     421.819875557     47.0630265273     -84.4234122781     178     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     5.59404944528     91.83938479099999     
  64     0.250827661965     402.70397995     47.0405903436     -84.408890007     221     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.88068375812     182.303418297     
  65     0.10676543058     398.861679373     47.0375238179     -84.4110044592     221     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.839318208581     180.968959292     
  66     0.45794361657399996     399.39783116     47.0503691305     -84.4234695415     263     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.80433503642     180.832283745     
  67     0.215521115747     384.343120205     47.0491405313     -84.4261542925     263     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.772875015059     185.009024539     
  68     0.17634116993100002     402.457873963     47.0473339874     -84.4231118057     263     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.51787473301     197.972194278     
  69     0.22380640934999999     389.880578267     47.0547527615     -84.4161868983     275     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.27469126547     158.424273605     
  70     0.0837350497009     374.446121029     47.0595904522     -84.428297166     288     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.81671425354     255.505953051     
  71     0.047399757802799994     345.128922986     47.0600775622     -84.425448652     288     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.655553848427     165.217669103     
  72     0.06670414607970002     383.84470344     47.0608943613     -84.4238262609     288     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.80964130543     152.11888638     
  73     0.00665796501057     544.097626527     47.0671708661     -84.3876001466     2     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.73118002992     166.965487872     
  74     0.020780870348     542.103681837     47.0670232098     -84.387888398     2     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.2507729186     141.268798314     
  75     0.019749447310399997     553.803902662     47.0642547163     -84.3846375727     11     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     9.22316892517     114.03548559999999     
  76     0.022171358451700002     544.831094274     47.0658984963     -84.385296233     11     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.84729069801     144.28987844     
  77     0.00278698129305     554.697559022     47.0626825328     -84.3848677876     16     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.01010793552     74.80636328999998     
  78     0.0103783062234     539.618229594     47.0631032826     -84.3855412966     16     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.99680067707     134.921591478     
  79     0.0104570963033     545.568761532     47.0628919991     -84.3853300628     16     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     10.9929187447     112.75028335100001     
  80     0.0408540509488     542.307077684     47.0670996666     -84.3978907956     29     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     9.1557441005     301.8069616411     
  81     0.0273609407563     497.3420386     47.0642727929     -84.3902834003     36     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.690157693787     174.865440842     
  82     0.0226905688113     524.385219313     47.0631095034     -84.3869907965     36     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.72038411726     123.38578729     
  83     0.00533807893731     497.45591794     47.0619752658     -84.3944073661     38     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.640878902442     156.991073133     
  84     0.008373074457070001     499.504320585     47.0618599257     -84.3945660473     38     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.25509789665     227.249361006     
  85     0.00184229479976     497.318296358     47.0621898018     -84.394928576     38     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.696158509387     182.709917688     
  86     0.00275896632002     502.260563521     47.0614739323     -84.3939966634     38     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.28641409224     130.928114804     
  87     0.0122844771265     496.830499318     47.0633641007     -84.394256441     40     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.747488545829     128.306275879     
  88     0.0103256961679     519.959838139     47.0639662698     -84.3951167547     40     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.29645321534     185.813868332     
  89     0.0245066744937     496.75239305     47.0637082411     -84.3946420108     40     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.725068498149     184.213822895     
  90     0.0388192558581     516.602124323     47.0646785651     -84.3960279972     40     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.44961863959     211.929557509     
  91     0.010973242861999999     530.341240276     47.0651248328     -84.3968091881     40     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.80314393152     237.678271692     
  92     0.00514571325012     487.954668017     47.0621120416     -84.3967210837     41     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.747196791494     175.613641903     
  93     0.0042942587493299995     494.750496212     47.0621236166     -84.3960695052     41     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.75512656048     133.391227418     
  94     0.00675136488191     487.918176896     47.0623891132     -84.3969081371     41     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.694156563172     178.425428307     
  95     0.00825473488448     492.297450045     47.0627052468     -84.3967470122     41     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.07189493546     171.339414532     
  96     0.0318917210193     533.279884798     47.0645015043     -84.4041226998     50     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.61660007428     273.6951111778     
  97     0.018682076739500002     487.050354126     47.06183643     -84.4002278277     59     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.80483184333     279.9085150426     
  98     0.0150430015325     480.70572347     47.0622567747     -84.398919357     59     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.83142193056     174.125180916     
  99     0.0482757365793     510.338195801     47.0635231197     -84.39863406     59     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.17122965684     191.949287197     
  100     0.0178699762252     521.831139567     47.0640489281     -84.3990287461     59     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.46850537058     192.648981737     
  101     0.058484581652199996     479.982786735     47.0585723674     -84.3953735049     73     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.69597612394     159.061423243     
  102     0.038697159616700005     493.099892117     47.058901326     -84.3967962291     73     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.64666544133     191.558721621     
  103     0.0307860803033     464.109801699     47.0533625159     -84.3950296213     93     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.7746049554     136.583570747     
  104     0.00920990187752     407.896172914     47.0473482551     -84.4094687539     106     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.448715433093     189.25444859     
  105     0.0174801852496     418.389601393     47.0469233014     -84.409233256     106     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.68814101497     203.990752464     
  106     0.0157553882656     421.479462204     47.0469825845     -84.4098172476     106     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     9.1955974591     223.378859217     
  107     0.00802005051489     408.879251803     47.0476948071     -84.4096157933     106     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.641775838083     184.656509611     
  108     0.0143746241057     418.272729205     47.0471872012     -84.4090220663     106     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.95935304533     169.583173436     
  109     0.0037967248635999997     400.549537854     47.0638818083     -84.4260062207     116     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.825990917209     135.466438189     
  110     0.00731605137782     403.763731447     47.0640126473     -84.4264769627     116     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.60595547826     153.565236131     
  111     0.0006939779651830001     482.952127391     47.0634403804     -84.4187273385     119     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.871227335313     109.44146490099999     
  112     0.00418910262297     369.730183613     47.0419888778     -84.4165742843     133     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.65995258636     202.574485416     
  113     0.00326629370282     389.881978826     47.0414916815     -84.416765912     133     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.21050129196     206.831495869     
  114     0.00421284221398     398.814280674     47.0418823091     -84.4156373736     133     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     9.20386556761     126.705382985     
  115     0.006604074220230001     399.614510525     47.0419794208     -84.4153275378     133     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.57401428659     107.00640389099999     
  116     0.0354279979543     420.634174563     47.039936473     -84.4189550369     134     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.91553011422     188.67517829     
  117     0.09043711486279998     439.663084712     47.0548341102     -84.4047066733     142     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.9584923565     227.407598172     
  118     0.037692874507     445.974558135     47.0544562734     -84.4061680418     142     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.43449851805     248.140673194     
  119     0.17872811198799998     442.637804205     47.056223513     -84.403900566     142     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.37419759619     182.329951699     
  120     0.04003806530660001     448.233898129     47.0557082256     -84.4039292298     142     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.61065654503     198.526053195     
  121     0.00828309045523     398.021061831     47.0480795783     -84.4067728306     150     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.661746535313     188.291479891     
  122     0.0288469486859     398.970229772     47.0481659885     -84.4067842873     150     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.71509775939     219.299458934     
  123     0.00879362801215     411.278176568     47.0475530943     -84.4068747981     150     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     9.03234998611     293.4015254253     
  124     0.00895462628665     390.93939833     47.0485297778     -84.4064692466     150     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.643904059006     170.60836889     
  125     0.0331331694503     401.568047755     47.0493233095     -84.4064694336     150     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.94980016736     174.958147161     
  126     0.0201138524736     378.481788902     47.0440992421     -84.3978715394     164     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.791322955349     141.940493529     
  127     0.08579361399789999     396.582307134     47.0435675121     -84.398844254     164     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.57419941681     242.746185374     
  128     0.08765903792610001     390.856278159     47.0454863188     -84.3982518547     164     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.55885975971     158.608625732     
  129     0.0196695877824     403.659679119     47.0451035266     -84.3964828716     164     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.6369128613     154.647347876     
  130     0.00656550952256     399.380850277     47.045577408     -84.4045441094     166     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.74602649139     296.9875612263     
  131     0.00149043906856     382.380780672     47.0455444546     -84.403617248     166     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.742815664628     179.266982898     
  132     0.0038723966872999996     404.522594034     47.045644294     -84.4047347356     166     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.52910539719     236.237217319     
  133     0.0015945522317299997     394.490939617     47.0455080109     -84.4042158038     166     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.67369093746     250.943524957     
  134     0.021871574620700003     393.850160682     47.0468269946     -84.4033167852     167     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.0068402176     206.378223546     
  135     0.0114509792874     400.919751658     47.0465257051     -84.4048722262     167     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.26520635367     300.815997418     
  136     0.00804217828426     385.737598278     47.0471839763     -84.4028454108     167     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.834166868887     181.085553299     
  137     0.0304990553407     392.092582052     47.0469369264     -84.4036649914     167     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.7029923535     180.919455001     
  138     0.0022106106957000002     387.550058764     47.0470652715     -84.3990606389     170     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.535172805136     202.569428142     
  139     0.008109557174300001     398.639678308     47.0471120032     -84.398547036     170     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.40384549137     147.780076119     
  140     0.00287207987106     403.605394398     47.0471047021     -84.3979249073     170     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.12608095082     69.45367709099997     
  141     0.00382793625908     388.595835168     47.0472803572     -84.3990826105     170     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.556580752723     177.6203901     
  142     0.013574930366600001     395.518056663     47.0469509045     -84.3986393722     170     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.44565961286     153.517122307     
  143     0.018516890584100002     428.79883532     47.0414163553     -84.3964215571     173     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.49007714767     324.0415547902     
  144     0.018424113115799997     415.225715492     47.0630267809     -84.4238860999     178     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.50721512125     140.944742446     
  145     0.0110614352942     444.800712794     47.0545482427     -84.4105726085     181     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.58931226785     57.22105688400001     
  146     0.0118707445427     434.891689903     47.0546558705     -84.4108326479     181     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.16040568427     72.76895327900002     
  147     0.00374521686889     454.684303083     47.0545348647     -84.4101990225     181     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.63348725987     106.385456788     
  148     0.0049678342111399995     437.884019861     47.0551156847     -84.4106239667     181     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.43135959839     97.54594965799998     
  149     0.13637797000100002     436.325345015     47.0580396153     -84.4094586914     195     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.67000867402     102.23347626399999     
  150     0.19960858080399999     488.205934765     47.0605667969     -84.40925433     195     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.41305020314     168.760522726     
  151     0.00685656431456     378.8669674     47.0441795454     -84.4027394298     205     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.93459382772     86.35691509100002     
  152     0.00199330817155     373.554113215     47.0443026299     -84.4032369414     205     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.645654160578     179.325179602     
  153     0.00569784867104     378.089039     47.0444505457     -84.4030754413     205     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.12974740099     158.370437588     
  154     0.10643305089     405.837234843     47.0420215617     -84.4070421868     221     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.17724593072     156.387626116     
  155     0.000544314445159     371.314527429     47.0442625941     -84.4118937868     222     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.895057636759     186.116371061     
  156     0.00118129007873     373.79030463     47.0441234128     -84.4117537465     222     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.22957915702     103.651208543     
  157     0.000230317842524     370.552285331     47.0443491249     -84.4118383427     222     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.700074521559     186.178913661     
  158     0.000318736167674     371.771006993     47.0444458065     -84.4119605965     222     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     2.43243925486     179.661129543     
  159     0.199077061961     396.681525427     47.050315896     -84.4214125754     263     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.70702416269     177.133335849     
  160     0.0651887500128     370.770626087     47.0547492836     -84.4195650467     275     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.86884553571     179.452106464     
  161     0.0881568722817     389.077670169     47.0537834638     -84.4168053013     275     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.30394886229     244.04126512300002     
  162     0.06554849125750001     414.177923142     47.0538003874     -84.4142793435     275     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     8.0753434835     133.324417424     
  163     0.0140847548348     346.226569547     47.0571523296     -84.4205679033     276     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.798505575499     101.76397093000003     
  164     0.0221915282596     365.373764244     47.0570559398     -84.4181114708     276     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     4.03504252132     79.54974892199999     
  165     0.019874306810699998     347.63837431     47.0577663068     -84.4205613782     276     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.884027102635     134.852015867     
  166     0.103541113942     383.168017136     47.0588951002     -84.4183654709     276     LU_Forest     Veg_Forest     SoilPf_Slp_2     [NONE]     [NONE]     3.78713141198     139.474322459     
  167     0.015207994895899999     425.033193135     47.0599086623     -84.4168196949     276     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     7.87443832854     146.605292379     
  168     0.0612158817001     344.326774816     47.0596361223     -84.4259377139     288     LU_Forest     Veg_Forest     SoilPf_Slp_1     [NONE]     [NONE]     0.639532703554     202.760933699     
  169     0.0697509421142     382.088895073     47.0587849657     -84.4274946666     288     LU_Forest     Veg_Forest     SoilPf_Slp_3     [NONE]     [NONE]     9.20082310429     294.4880303834     
:EndHRUs

:HRUGroup C31_C32_HRUS
22,54,109,110,27,63,144
:EndHRUGroup

:PopulateHRUGroup Lake_HRUs With LANDUSE EQUALS LU_Lake 
:PopulateHRUGroup Non_C31_C32_HRUS With HRUS NOTWITHIN C31_C32_HRUS
    
#:RedirectToFile Turkey_Lake_Lake.rvh
:SubBasinGroup   Allsubbasins
       29   36   38   59   103   163   164   221   263   53
       2   11   16   50   73   93   115   116   119   134
       173   178   195   275   288   40   41   106   133   142
       150   166   167   170   181   205   222   276
:EndSubBasinGroup   
# :SBGroupPropertyOverride Allsubbasins   MANNINGS_N 0.001
# :SBGroupPropertyMultiplier Allsubbasins  MANNINGS_N 1.0
:SubBasinGroup   AllLakesubbasins
       29   36   38   59   103   163   164   221   263   53
:EndSubBasinGroup   
# :SBGroupPropertyOverride   AllLakesubbasins   RESERVOIR_CREST_WIDTH 12.0
# :SBGroupPropertyMultiplier  AllLakesubbasins   RESERVOIR_CREST_WIDTH 1.0
 
"""

RVH_LAKE = """
#----------------------------------------------
# This is a Raven lake rvh file generated
# by Routing toolbox
#----------------------------------------------
:Reservoir   Lake_10   ######## 
  :SubBasinID  29
  :HRUID   1
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    66125.9212345
:EndReservoir   
#############################################
###New Lake starts
:Reservoir   Lake_14   ######## 
  :SubBasinID  36
  :HRUID   2
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    63238.4416956
:EndReservoir   
#############################################
###New Lake starts
:Reservoir   Lake_11   ######## 
  :SubBasinID  38
  :HRUID   3
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    5408.73851459
:EndReservoir   
#############################################
###New Lake starts
:Reservoir   Lake_15   ######## 
  :SubBasinID  59
  :HRUID   4
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    16560.7112489
:EndReservoir   
#############################################
###New Lake starts
:Reservoir   Lake_20   ######## 
  :SubBasinID  103
  :HRUID   5
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    50034.4610312
:EndReservoir   
#############################################
###New Lake starts
:Reservoir   Lake_12   ######## 
  :SubBasinID  163
  :HRUID   6
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    194230.601405
:EndReservoir   
#############################################
###New Lake starts
:Reservoir   Lake_18   ######## 
  :SubBasinID  164
  :HRUID   7
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    38348.6959141
:EndReservoir   
#############################################
###New Lake starts
:Reservoir   Lake_19   ######## 
  :SubBasinID  221
  :HRUID   8
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    212907.376242
:EndReservoir   
#############################################
###New Lake starts
:Reservoir   Lake_17   ######## 
  :SubBasinID  263
  :HRUID   9
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    558974.609372
:EndReservoir   
#############################################
###New Lake starts
:Reservoir   Lake_16   ######## 
  :SubBasinID  53
  :HRUID   10
  :Type RESROUTE_STANDARD   
  :WeirCoefficient  0.6
  :CrestWidth 1.2345
  :MaxDepth 3.0
  :LakeArea    2629.62940252
:EndReservoir   
#############################################
###New Lake starts
"""

RVI = """
# ----------------------------------------------
# Raven Input file
# HBV-EC Nith River emulation test case
# ----------------------------------------------
# --Simulation Details -------------------------
:StartDate       1981-01-01 00:00:00
:Duration        12050 #12050 #5844  #12050   # 5844 to 1996               # 1 1 1981 to 12 31 2016
:Method          ORDERED_SERIES
:TimeStep        1.0
#:RunName         calibration
:OutputDirectory 	 output/ 
#
# --Model Details -------------------------------
:Method                 ORDERED_SERIES
:Interpolation          INTERP_NEAREST_NEIGHBOR
:SoilModel              SOIL_MULTILAYER 3

:VegetationProcesses Robin

:Routing                ROUTE_DIFFUSIVE_WAVE #ROUTE_HYDROLOGIC 
:CatchmentRoute         ROUTE_TRI_CONVOLUTION
:OW_Evaporation         PET_PRIESTLEY_TAYLOR
:RainSnowFraction       RAINSNOW_DATA

##Vegetation impacted processes 
:SWCanopyCorrect        SW_CANOPY_CORR_NONE    #SW_CANOPY_CORR_STATIC  #SW_CANOPY_CORR_NONE #SW_CANOPY_CORR_DYNAMIC
:PrecipIceptFract       PRECIP_ICEPT_LAI
:Evaporation            PET_PRIESTLEY_TAYLOR   #PET_PENMAN_MONTEITH  #
:PotentialMeltMethod    POTMELT_HBV            #POTMELT_EB #POTMELT_HBV #POTMELT_EB

###about forcing data
:RelativeHumidityMethod RELHUM_DATA
:SWRadiationMethod      SW_RAD_DATA
:LWRadiationMethod      LW_RAD_DEFAULT
:WindspeedMethod        WINDVEL_DATA


#
# --Hydrologic Processes-------------------------
:Alias       Forest_Floor   SOIL[0]
:Alias       Ablation_Till  SOIL[1]
:Alias       Basal_Till     SOIL[2]
#

:DefineHRUGroups  C31_C32_HRUS  Non_C31_C32_HRUS
:DisableHRUGroup Non_C31_C32_HRUS
:EvaluationPeriod CALIBRATION 1984-01-01 1996-12-01
:EvaluationPeriod VALIDATION 1997-09-01 2012-12-01 
:EvaluationPeriod VALIDATION_H 1998-01-01 2001-12-01 


:HydrologicProcesses

  :SnowRefreeze      FREEZE_DEGREE_DAY  SNOW_LIQ         SNOW
  :Precipitation     PRECIP_RAVEN       ATMOS_PRECIP     MULTIPLE
  
  :CanopyDrip        CANDRIP_RUTTER      CANOPY           PONDED_WATER
  :CanopyEvaporation CANEVP_MAXIMUM      CANOPY           ATMOSPHERE
  :CanopySnowEvap    CANEVP_MAXIMUM      CANOPY_SNOW      ATMOSPHERE
  :SnowBalance       SNOBAL_SIMPLE_MELT SNOW             SNOW_LIQ                    #:SnowBalance       SNOBAL_TWO_LAYER   MULTIPLE         MULTIPLE     # :SnowBalance       SNOBAL_SIMPLE_MELT SNOW             PONDED_WATER
  :-->Overflow       RAVEN_DEFAULT      SNOW_LIQ         PONDED_WATER                                                                                                                                                            #:-->Overflow     RAVEN_DEFAULT      SNOW_LIQ         PONDED_WATER
  :Infiltration      INF_GREEN_AMPT     PONDED_WATER     MULTIPLE
  ### Forest_Floor
  :SoilEvaporation   SOILEVAP_ROOT_VEG        Forest_Floor     ATMOSPHERE
  :Percolation       PERC_GAWSER          Forest_Floor     Ablation_Till   
  :Baseflow          BASE_THRESH_POWER    Forest_Floor     SURFACE_WATER
  ### Ablation_Till
  :Baseflow          BASE_THRESH_POWER     Ablation_Till    SURFACE_WATER
  :Percolation       PERC_GAWSER           Ablation_Till    Basal_Till
  ### Basal_Till
  :Baseflow          BASE_THRESH_POWER  Basal_Till       SURFACE_WATER
  :CapillaryRise     CRISE_HBV          Basal_Till       Ablation_Till
  
:EndHydrologicProcesses



#------------------------------------------------------------------------# Transport calculations for hydrograph separation
# FixedConcentrations 2,3,4 are boundary conditions specifying that water passing 
# through to underlying compartments get a concentration of 0 mg/L (avoid double-counting)

#:Transport BASAL_TILL
#:FixedConcentration BASAL_TILL Basal_Till 1.0

#:Transport FOREST_FLOOR
#:FixedConcentration FOREST_FLOOR Forest_Floor 1.0
#:FixedConcentration FOREST_FLOOR Ablation_Till 0.0

#:Transport ABLATION_TILL
#:FixedConcentration ABLATION_TILL Ablation_Till 1.0
#:FixedConcentration ABLATION_TILL Basal_Till 0.0

# compare results from this tracer with Q_total - Q_baseflow - they should be identical.
#:Transport OVERLAND_FLOW
#:FixedConcentration OVERLAND_FLOW PONDED_WATER 1.0
#:FixedConcentration OVERLAND_FLOW Forest_Floor 0.0 
#:FixedConcentration OVERLAND_FLOW Ablation_Till 0.0
#:FixedConcentration OVERLAND_FLOW Basal_Till 0.0



#
# --Output Options-------------------------------
#
#:CreateRVPTemplate
:SilentMode
#:SuppressOutput 
#:WriteForcingFunctions 
:EvaluationMetrics NASH_SUTCLIFFE  PCT_BIAS KLING_GUPTA RMSE KLING_GUPTA_SUMMER KLING_GUPTA_NOT_SUMMER
#:WriteMassBalanceFile 

"""

RVP = """
#-----------------------------------------------------------------
# Raven Properties file Template. Created by Raven v2.9 w/ netCDF
#-----------------------------------------------------------------
# all expressions of format *xxx* need to be specified by the user 
# all parameter values of format ** need to be specified by the user 
# soil, land use, and vegetation classes should be made consistent with user-generated .rvh file 
#-----------------------------------------------------------------
:AvgAnnualRunoff  100

#-----------------------------------------------------------------
# Soil Classes
#-----------------------------------------------------------------
:SoilClasses          
:Attributes		      SAND	  %CLAY	 %SILT	 %ORGANIC
	:Units			  none	  none	 none	 none
	Forest_Floor,     0.66,   0.1,  0.24,    0.05,    
	Ablation_Till,    0.66,   0.1,  0.24,    0.03,    
	Basal_Till,       0.22,   0.13, 0.65,    0.02,    
:EndSoilClasses

#-----------------------------------------------------------------
# Land Use Classes
#-----------------------------------------------------------------
:LandUseClasses, 
  :Attributes,        IMPERM,          FOREST_COV, 
       :Units,          frac,          frac, 
  LU_Forest,            0,              1, 
  LU_Lake,              0,              0,         
:EndLandUseClasses

#-----------------------------------------------------------------
# Vegetation Classes
#-----------------------------------------------------------------

:VegetationClasses, 
  :Attributes,        MAX_HT,       MAX_LAI,   MAX_LEAF_COND, 
       :Units,             m,          none,        mm_per_s, 
    Veg_Forest,           15,           7.5,             5.3,     
    Veg_Lake,              0,             0,               0,          
:EndVegetationClasses  

#-----------------------------------------------------------------
# Soil Profiles
#-----------------------------------------------------------------
:SoilProfiles
         LAKE, 0
         ROCK, 0
   SoilPf_Slp_1,      3, Forest_Floor,  0.27445793,  	Ablation_Till,  0.44388355,   Basal_Till,        3, 
   SoilPf_Slp_2,      3, Forest_Floor,  0.17445793,    Ablation_Till,  0.3998468,   Basal_Till,        3, 
   SoilPf_Slp_3,      3, Forest_Floor,  0.1,  	Ablation_Till,  0.3898468,   Basal_Till,        3,
   
:EndSoilProfiles

#-----------------------------------------------------------------
# Global Parameters
#-----------------------------------------------------------------
:GlobalParameter        RAINSNOW_TEMP 0.07195228
:GlobalParameter       RAINSNOW_DELTA 2.855478
:GlobalParameter     SNOW_TEMPERATURE -2.9
:GlobalParameter             SNOW_SWI 0.05 #(Raven manual 0.04 - 0.07) 

#-----------------------------------------------------------------
# Soil Parameters
#-----------------------------------------------------------------
###
#POROSITY and FIELD_CAPACITY HYDRAUL_COND for ablation till obtained from (Figure2 Murray2005)
#ALBEDO_WET ALBEDO_DRY obtained from Figure 2 of https://en.wikipedia.org/wiki/Albedo ; https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GL036377   
#WETTING_FRONT_PSI table 2 of  https://ascelibrary.org/doi/pdf/10.1061/%28ASCE%290733-9429%281983%29109%3A1%2862%29?casa_token=Db6omitRaNAAAAAA:NJEfRhzP2102wLfyfQ6-niyrgad8YenXl3aAK5eU_FnPN9Aorx6EV31zTiGxUW07XFYGQV_t-ODb
####
:SoilParameterList
  :Parameters,        POROSITY,      ALBEDO_WET,      ALBEDO_DRY,    HYDRAUL_COND,      WETTING_FRONT_PSI,  PET_CORRECTION,   FIELD_CAPACITY,        SAT_WILT,       MAX_PERC_RATE,                 MAX_BASEFLOW_RATE,  BASEFLOW_COEFF,       BASEFLOW_N,    MAX_CAP_RISE_RATE, BASEFLOW_THRESH 
       :Units,            none,            none,            none,            mm/d,                     mm,                -,               -,               -,                mm/d,                               mm/d,            none,             none,                mm/d,            none 
Forest_Floor,     6.543765E-01,            0.15,            0.15,    977.6504,                 14.94443,                1,         0.2482568,    1.007189E-02,  65.30797,             100.0,             0.0,  0.4242617,                 0.0,         0.2482568,  
Ablation_Till,    6.882158E-01,            0.00,            0.00,    5.283446E+02,                 14.94443,                1,         0.4098241,    5.475422E-02,  11.58644,             115.4073,    5.168251E-02,  2.287313,                 0.0,         0.4098241,
Basal_Till                0.2,             0.00,            0.00,            9.52,                    0.0,                0,         0.5,    5.764653E-03,        7.501124E-01,             221.6116,       0.9351502,  5.136656, 13.12126,         0.5,
:EndSoilParameterList


#-----------------------------------------------------------------
# Land Use Parameters
#-----------------------------------------------------------------
:LandUseParameterList
  :Parameters, FOREST_SPARSENESS,       ROUGHNESS,    REFREEZE_FACTOR,    MELT_FACTOR,      MIN_MELT_FACTOR , HBV_MELT_FOR_CORR, HBV_MELT_ASP_CORR,    CC_DECAY_COEFF
       :Units,               -,                 m,             mm/d/C,        mm/d/K ,             mm/d/K,                 none,              none,             1/day
   [DEFAULT] ,                0,               1.9,         0.0,        2.89891738,        2.887848,                    1,     0.3772205,                 0  
   LU_Forest,                 0,               1.9,          _DEFAULT,       _DEFAULT,           _DEFAULT,             _DEFAULT,          _DEFAULT,          _DEFAULT 
   LU_Lake,                  1,                 0,           _DEFAULT,       _DEFAULT,           _DEFAULT,             _DEFAULT,          _DEFAULT,          _DEFAULT 
:EndLandUseParameterList

#-----------------------------------------------------------------
# Vegetation Parameters
#-----------------------------------------------------------------
## SAI_HT_RATIO did not find any information. use zero to make ignore the SAI, the canopy interception only 
## determined by LAI
:VegetationParameterList
  :Parameters,      MAX_HEIGHT,   MAX_LEAF_COND,       ALBEDO,  SVF_EXTINCTION,    SAI_HT_RATIO,    RAIN_ICEPT_FACT,    SNOW_ICEPT_FACT,    MAX_CAPACITY,    MAX_SNOW_CAPACITY, 
       :Units,               m,        mm_per_s,            -,               -,               -,                  -,                  -,              mm,                   mm, 
    Veg_Forest              15,             5.3,         0.15,              0.5,            0.0,  0.3318293,               0.00,  16.62024,                 0.00, 
    Veg_Lake,                0,             0.0,         0.06,              0.5,            0.0,                0.0,                0.0,               0,                    0,  
:EndVegetationParameterList

:SeasonalRelativeLAI
  Veg_Forest, 0.2,0.2,0.2,0.4,0.8,1,1,0.7,0.6,0.0.3,0.0.2,0.0.2
  Veg_Lake, 0.5,0.6,0.7,1,1,1,1,0.7,0.6,0.6,0.6,0.5
:EndSeasonalRelativeLAI

:SeasonalRelativeHeight
  Veg_Forest, 1,1,1,1,1,1,1,1,1,1,1,1
  Veg_Lake, 0.5,0.6,0.7,1,1,1,1,0.7,0.6,0.6,0.6,0.5
:EndSeasonalRelativeHeight

:RedirectToFile Turkey_Lake_channel.rvp
"""

RVP_CHANNEL = """
#----------------------------------------------
# This is a Raven channel properties file generated
# by Routing toolbox
#----------------------------------------------
:ChannelProfile          Chn_29          
  :Bedslope          0.0100908642754
  :SurveyPoints
    0          525.343033987
    16.0          521.343033987
    18.469          521.343033987
    18.777625          520.108533987
    19.394875          520.108533987
    19.7035          521.343033987
    22.172499999999996          521.343033987
    38.1725          525.343033987
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0391944629923
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_36          
  :Bedslope          0.0100908642754
  :SurveyPoints
    0          510.982653339
    16.0          506.982653339
    18.469          506.982653339
    18.777625          505.748153339
    19.394875          505.748153339
    19.7035          506.982653339
    22.172499999999996          506.982653339
    38.1725          510.982653339
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0391944629923
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_38          
  :Bedslope          0.318345272341
  :SurveyPoints
    0          501.816108712
    16.0          497.816108712
    18.469          497.816108712
    18.777625          496.581608712
    19.394875          496.581608712
    19.7035          497.816108712
    22.172499999999996          497.816108712
    38.1725          501.816108712
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0982156219482
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_59          
  :Bedslope          0.107388234734
  :SurveyPoints
    0          503.965449486
    16.0          499.965449486
    18.469          499.965449486
    18.777625          498.730949486
    19.394875          498.730949486
    19.7035          499.965449486
    22.172499999999996          499.965449486
    38.1725          503.965449486
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0718840076687
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_103          
  :Bedslope          0.0644325210483
  :SurveyPoints
    0          427.21447922
    16.0          423.21447922
    18.469          423.21447922
    18.777625          421.97997921999996
    19.394875          421.97997921999996
    19.7035          423.21447922
    22.172499999999996          423.21447922
    38.1725          427.21447922
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.119253449079
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_163          
  :Bedslope          0.0257457088668
  :SurveyPoints
    0          421.773949103
    16.0          417.773949103
    18.469          417.773949103
    18.777625          416.539449103
    19.394875          416.539449103
    19.7035          417.773949103
    22.172499999999996          417.773949103
    38.1725          421.773949103
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0382734077682
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_164          
  :Bedslope          1.8736411823e-05
  :SurveyPoints
    0          394.191028831
    16.0          390.191028831
    18.469          390.191028831
    18.777625          388.95652883099996
    19.394875          388.95652883099996
    19.7035          390.191028831
    22.172499999999996          390.191028831
    38.1725          394.191028831
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.01
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_221          
  :Bedslope          0.00523219383959
  :SurveyPoints
    0          399.501884619
    16.0          395.501884619
    18.469          395.501884619
    18.777625          394.267384619
    19.394875          394.267384619
    19.7035          395.501884619
    22.172499999999996          395.501884619
    38.1725          399.501884619
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0249872664141
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_263          
  :Bedslope          0.0119156500919
  :SurveyPoints
    0          393.660666209
    16.0          389.660666209
    18.469          389.660666209
    18.777625          388.426166209
    19.394875          388.426166209
    19.7035          389.660666209
    22.172499999999996          389.660666209
    38.1725          393.660666209
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0170151145712
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_53          
  :Bedslope          0.179015843091
  :SurveyPoints
    0          494.550656902
    16.0          490.550656902
    18.469          490.550656902
    18.777625          489.316156902
    19.394875          489.316156902
    19.7035          490.550656902
    22.172499999999996          490.550656902
    38.1725          494.550656902
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.121361546403
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_2          
  :Bedslope          0.230660950822
  :SurveyPoints
    0          558.349418734
    16.0          554.349418734
    18.469          554.349418734
    18.777625          553.114918734
    19.394875          553.114918734
    19.7035          554.349418734
    22.172499999999996          554.349418734
    38.1725          558.349418734
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.15
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_11          
  :Bedslope          0.104755347593
  :SurveyPoints
    0          557.750606883
    16.0          553.750606883
    18.469          553.750606883
    18.777625          552.516106883
    19.394875          552.516106883
    19.7035          553.750606883
    22.172499999999996          553.750606883
    38.1725          557.750606883
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.139480768379
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_16          
  :Bedslope          0.524660628693
  :SurveyPoints
    0          555.288330866
    16.0          551.288330866
    18.469          551.288330866
    18.777625          550.053830866
    19.394875          550.053830866
    19.7035          551.288330866
    22.172499999999996          551.288330866
    38.1725          555.288330866
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0929872352239
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_50          
  :Bedslope          0.125478272608
  :SurveyPoints
    0          548.207855579
    16.0          544.207855579
    18.469          544.207855579
    18.777625          542.973355579
    19.394875          542.973355579
    19.7035          544.207855579
    22.172499999999996          544.207855579
    38.1725          548.207855579
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.124431025288
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_73          
  :Bedslope          0.0953753264322
  :SurveyPoints
    0          490.418913852
    16.0          486.418913852
    18.469          486.418913852
    18.777625          485.184413852
    19.394875          485.184413852
    19.7035          486.418913852
    22.172499999999996          486.418913852
    38.1725          490.418913852
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.114069683381
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_93          
  :Bedslope          0.179015843091
  :SurveyPoints
    0          476.027638818
    16.0          472.027638818
    18.469          472.027638818
    18.777625          470.793138818
    19.394875          470.793138818
    19.7035          472.027638818
    22.172499999999996          472.027638818
    38.1725          476.027638818
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.121361546403
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_115          
  :Bedslope          0.381299119196
  :SurveyPoints
    0          371.29358581
    16.0          367.29358581
    18.469          367.29358581
    18.777625          366.05908581
    19.394875          366.05908581
    19.7035          367.29358581
    22.172499999999996          367.29358581
    38.1725          371.29358581
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.15
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_116          
  :Bedslope          0.251605310115
  :SurveyPoints
    0          408.517704785
    16.0          404.517704785
    18.469          404.517704785
    18.777625          403.28320478499995
    19.394875          403.28320478499995
    19.7035          404.517704785
    22.172499999999996          404.517704785
    38.1725          408.517704785
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.15
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_119          
  :Bedslope          0.218572964499
  :SurveyPoints
    0          474.237166886
    16.0          470.237166886
    18.469          470.237166886
    18.777625          469.002666886
    19.394875          469.002666886
    19.7035          470.237166886
    22.172499999999996          470.237166886
    38.1725          474.237166886
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.15
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_134          
  :Bedslope          0.0644325210483
  :SurveyPoints
    0          418.789688416
    16.0          414.789688416
    18.469          414.789688416
    18.777625          413.55518841599996
    19.394875          413.55518841599996
    19.7035          414.789688416
    22.172499999999996          414.789688416
    38.1725          418.789688416
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.119253449079
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_173          
  :Bedslope          0.07745972134
  :SurveyPoints
    0          418.89245378
    16.0          414.89245378
    18.469          414.89245378
    18.777625          413.65795377999996
    19.394875          413.65795377999996
    19.7035          414.89245378
    22.172499999999996          414.89245378
    38.1725          418.89245378
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.115194742223
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_178          
  :Bedslope          0.271268490128
  :SurveyPoints
    0          416.58781382
    16.0          412.58781382
    18.469          412.58781382
    18.777625          411.35331382
    19.394875          411.35331382
    19.7035          412.58781382
    22.172499999999996          412.58781382
    38.1725          416.58781382
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.15
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_195          
  :Bedslope          0.104496289369
  :SurveyPoints
    0          477.840204763
    16.0          473.840204763
    18.469          473.840204763
    18.777625          472.60570476299995
    19.394875          472.60570476299995
    19.7035          473.840204763
    22.172499999999996          473.840204763
    38.1725          477.840204763
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.106343384711
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_275          
  :Bedslope          0.0221767491895
  :SurveyPoints
    0          393.320472141
    16.0          389.320472141
    18.469          389.320472141
    18.777625          388.08597214099996
    19.394875          388.08597214099996
    19.7035          389.320472141
    22.172499999999996          389.320472141
    38.1725          393.320472141
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0687385987483
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_288          
  :Bedslope          0.0064282797319
  :SurveyPoints
    0          373.966661209
    16.0          369.966661209
    18.469          369.966661209
    18.777625          368.73216120899997
    19.394875          368.73216120899997
    19.7035          369.966661209
    22.172499999999996          369.966661209
    38.1725          373.966661209
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0351680880579
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_40          
  :Bedslope          0.0586318627713
  :SurveyPoints
    0          514.9887169680001
    16.0          510.988716968
    18.469          510.988716968
    18.777625          509.754216968
    19.394875          509.754216968
    19.7035          510.988716968
    22.172499999999996          510.988716968
    38.1725          514.9887169680001
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0454876694712
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_41          
  :Bedslope          0.0624422833113
  :SurveyPoints
    0          494.59626621
    16.0          490.59626621
    18.469          490.59626621
    18.777625          489.36176621
    19.394875          489.36176621
    19.7035          490.59626621
    22.172499999999996          490.59626621
    38.1725          494.59626621
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.117397209569
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_106          
  :Bedslope          0.103506433286
  :SurveyPoints
    0          420.429215123
    16.0          416.429215123
    18.469          416.429215123
    18.777625          415.194715123
    19.394875          415.194715123
    19.7035          416.429215123
    22.172499999999996          416.429215123
    38.1725          420.429215123
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.124829254432
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_133          
  :Bedslope          0.245020481957
  :SurveyPoints
    0          394.840479308
    16.0          390.840479308
    18.469          390.840479308
    18.777625          389.605979308
    19.394875          389.605979308
    19.7035          390.840479308
    22.172499999999996          390.840479308
    38.1725          394.840479308
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.100298459328
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_142          
  :Bedslope          0.0765436924053
  :SurveyPoints
    0          446.882977289
    16.0          442.882977289
    18.469          442.882977289
    18.777625          441.64847728899997
    19.394875          441.64847728899997
    19.7035          442.882977289
    22.172499999999996          442.882977289
    38.1725          446.882977289
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.128344177038
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_150          
  :Bedslope          0.099467565886
  :SurveyPoints
    0          404.286507463
    16.0          400.286507463
    18.469          400.286507463
    18.777625          399.052007463
    19.394875          399.052007463
    19.7035          400.286507463
    22.172499999999996          400.286507463
    38.1725          404.286507463
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0876988594996
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_166          
  :Bedslope          0.163465425479
  :SurveyPoints
    0          402.523708272
    16.0          398.523708272
    18.469          398.523708272
    18.777625          397.289208272
    19.394875          397.289208272
    19.7035          398.523708272
    22.172499999999996          398.523708272
    38.1725          402.523708272
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0814093887996
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_167          
  :Bedslope          0.0435837051062
  :SurveyPoints
    0          397.336627528
    16.0          393.336627528
    18.469          393.336627528
    18.777625          392.102127528
    19.394875          392.102127528
    19.7035          393.336627528
    22.172499999999996          393.336627528
    38.1725          397.336627528
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0980800006884
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_170          
  :Bedslope          1.17206255463
  :SurveyPoints
    0          399.655767072
    16.0          395.655767072
    18.469          395.655767072
    18.777625          394.421267072
    19.394875          394.421267072
    19.7035          395.655767072
    22.172499999999996          395.655767072
    38.1725          399.655767072
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.15
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_181          
  :Bedslope          3.80516598977
  :SurveyPoints
    0          445.293222942
    16.0          441.293222942
    18.469          441.293222942
    18.777625          440.058722942
    19.394875          440.058722942
    19.7035          441.293222942
    22.172499999999996          441.293222942
    38.1725          445.293222942
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.115556912026
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_205          
  :Bedslope          0.0261419418762
  :SurveyPoints
    0          381.865484233
    16.0          377.865484233
    18.469          377.865484233
    18.777625          376.630984233
    19.394875          376.630984233
    19.7035          377.865484233
    22.172499999999996          377.865484233
    38.1725          381.865484233
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.0735809432661
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_222          
  :Bedslope          0.0461639648126
  :SurveyPoints
    0          376.658721588
    16.0          372.658721588
    18.469          372.658721588
    18.777625          371.42422158799997
    19.394875          371.42422158799997
    19.7035          372.658721588
    22.172499999999996          372.658721588
    38.1725          376.658721588
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.100941542553
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
:ChannelProfile          Chn_276          
  :Bedslope          0.0643051252693
  :SurveyPoints
    0          381.569968521
    16.0          377.569968521
    18.469          377.569968521
    18.777625          376.335468521
    19.394875          376.335468521
    19.7035          377.569968521
    22.172499999999996          377.569968521
    38.1725          381.569968521
  :EndSurveyPoints
  :RoughnessZones
    0          0.10000000149
    18.469          0.110045125251
    19.7035          0.10000000149
  :EndRoughnessZones
:EndChannelProfile


##############new channel ##############################
 
"""

RVT = """
# --------------------------------------------
# Raven Time Series Input file
# --------------------------------------------
#
# --------------------------------------------
#   Meteorological Gauges
# --------------------------------------------



:Gauge Mid Station
  :RainCorrection  1
  :SnowCorrection  0.8
  :RedirectToFile ./Forcing/Mid_2021-02-22.rvt
:EndGauge



# --------------------------------------------
#  Stream Gauges
# --------------------------------------------

:RedirectToFile ./obs/C31_Sub116.rvt
:RedirectToFile ./obs/C32_Sub178.rvt
:RedirectToFile ./obs/C33_Sub119.rvt
:RedirectToFile ./obs/C34_Sub195.rvt
:RedirectToFile ./obs/C37_Sub150.rvt
:RedirectToFile ./obs/C38_Sub106.rvt
:RedirectToFile ./obs/C39_Sub173.rvt
:RedirectToFile ./obs/C47_Sub16.rvt
:RedirectToFile ./obs/C49_Sub11.rvt
:RedirectToFile ./obs/C50_Sub2.rvt


:RedirectToFile ./obs/C46_Sub50.rvt
:RedirectToFile ./obs/C35_Sub181.rvt
:RedirectToFile ./obs/C42_Sub93.rvt
:RedirectToFile ./obs/S6_Sub134.rvt
:RedirectToFile ./obs/S4_Sub263.rvt
:RedirectToFile ./obs/S3_Sub221.rvt
:RedirectToFile ./obs/S2_Sub167.rvt
:RedirectToFile ./obs/S0_Sub36.rvt
:RedirectToFile ./obs/S1_Sub73.rvt
#:RedirectToFile ./obs/S5_Sub317.rvt

:RedirectToFile ./obs/C31_Sub116_ObsWts.rvt
:RedirectToFile ./obs/C32_Sub178_ObsWts.rvt
:RedirectToFile ./obs/C33_Sub119_ObsWts.rvt
:RedirectToFile ./obs/C34_Sub195_ObsWts.rvt
:RedirectToFile ./obs/C35_Sub181_ObsWts.rvt
:RedirectToFile ./obs/C37_Sub150_ObsWts.rvt
:RedirectToFile ./obs/C38_Sub106_ObsWts.rvt
:RedirectToFile ./obs/C39_Sub173_ObsWts.rvt
:RedirectToFile ./obs/C42_Sub93_ObsWts.rvt
:RedirectToFile ./obs/C46_Sub50_ObsWts.rvt
:RedirectToFile ./obs/C47_Sub16_ObsWts.rvt
:RedirectToFile ./obs/C49_Sub11_ObsWts.rvt
:RedirectToFile ./obs/C50_Sub2_ObsWts.rvt
:RedirectToFile ./obs/S6_Sub134_ObsWts.rvt
:RedirectToFile ./obs/S4_Sub263_ObsWts.rvt
:RedirectToFile ./obs/S3_Sub221_ObsWts.rvt
:RedirectToFile ./obs/S2_Sub167_ObsWts.rvt
:RedirectToFile ./obs/S0_Sub36_ObsWts.rvt
:RedirectToFile ./obs/S1_Sub73_ObsWts.rvt
#:RedirectToFile ./obs/S5_Sub317_ObsWts.rvt
"""
