#!/usr/bin/env python
from __future__ import print_function

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
#
#
# python 3_derive_elementary_effects.py \
#                       -i example_ishigami-homma/model_output.pkl \
#                       -d example_ishigami-homma/parameters.dat \
#                       -m example_ishigami-homma/parameter_sets_1_para3_M.dat \
#                       -v example_ishigami-homma/parameter_sets_1_para3_v.dat  \
#                       -o example_ishigami-homma/eee_results.dat

"""
Derives the Elementary Effects based on model outputs stored as dictionary in a pickle file (option -i)
using specified model parameters (option -d). The model parameters were sampled beforehand as Morris
trajectories. The Morris trajectory information is stored in two files (option -m and option -v). The
Elementary Effects are stored in a file (option -o).

History
-------
Written,  JM, Mar 2019
"""



# -------------------------------------------------------------------------
# Command line arguments
#
modeloutputs   = 'example_ishigami-homma/model_output.pkl'
modeloutputkey = 'All'
maskfile       = 'example_ishigami-homma/parameters.dat'
morris_M       = 'example_ishigami-homma/parameter_sets_1_para3_M.dat'
morris_v       = 'example_ishigami-homma/parameter_sets_1_para3_v.dat'
outfile        = 'example_ishigami-homma/eee_results.dat'
skip           = None                                                              # number of lines to skip in Morris files

import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="Derives the Elementary Effects based on model outputs stored as dictionary in a pickle file (option -i) using specified model parameters (option -d). The model parameters were sampled beforehand as Morris trajectories. The Morris trajectory information is stored in two files (option -m and option -v). The Elementary Effects are stored in a file (option -o).")

parser.add_option('-i', '--modeloutputs', action='store',
                    default=modeloutputs, dest='modeloutputs', metavar='modeloutputs',
                    help="Name of file used to save (scalar) model outputs in a pickle file (default: 'model_output.pkl').")
parser.add_option('-k', '--modeloutputkey', action='store',
                    default=modeloutputkey, dest='modeloutputkey', metavar='modeloutputkey',
                    help="Key of model output dictionary stored in pickle output file. If 'All', all model outputs are taken into account and multi-objective EEE is applied. (default: 'All').")
parser.add_option('-d', '--maskfile', action='store', dest='maskfile', type='string',
                  default=maskfile, metavar='File',
                  help='Name of file where all model parameters are specified including their distribution, distribution parameters, default value and if included in analysis or not. (default: maskfile=parameters.dat).')
parser.add_option('-m', '--morris_M', action='store', dest='morris_M', type='string',
                  default=morris_M, metavar='morris_M',
                  help="Morris trajectory information: The UNSCALED parameter sets. (default: 'parameter_sets_1_para3_M.dat').")
parser.add_option('-v', '--morris_v', action='store', dest='morris_v', type='string',
                  default=morris_v, metavar='morris_v',
                  help="Morris trajectory information: The indicator which parameter changed between subsequent sets in a trajectory. (default: 'parameter_sets_1_para3_v.dat').")
parser.add_option('-s', '--skip', action='store',
                    default=skip, dest='skip', metavar='skip',
                    help="Number of lines to skip in Morris output files (default: None).")
parser.add_option('-o', '--outfile', action='store', dest='outfile', type='string',
                  default=outfile, metavar='File',
                  help='File containing Elementary Effect estimates of all model parameters listed in parameter information file. (default: eee_results.dat).')
(opts, args) = parser.parse_args()

modeloutputs   = opts.modeloutputs
modeloutputkey = opts.modeloutputkey
maskfile       = opts.maskfile
morris_M       = opts.morris_M
morris_v       = opts.morris_v
outfile        = opts.outfile
skip           = opts.skip

del parser, opts, args


# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/lib')

import numpy       as np
import pickle
from fsread          import fsread              # in lib/
from autostring      import astr                # in lib/


# -------------------------
# read parameter info file
# -------------------------
# parameter info file has following header:
#       # para   dist       lower     upper     default   informative(0)_or_noninformative(1)
#       #                   mean      stddev
nc,snc = fsread(maskfile, comment="#",cskip=1,snc=[0,1],nc=[2,3,4,5])
snc = np.array(snc)
para_name   = snc[:,0]
para_dist   = snc[:,1]
lower_bound = nc[:,0]
upper_bound = nc[:,1]
initial     = nc[:,2]
# if informative(0)    -> maskpara=False
# if noninformative(1) -> maskpara=True
mask_para = np.where((nc[:,3].flatten())==1.,True,False)

dims_all  = np.shape(mask_para)[0]
idx_para  = np.arange(dims_all)[mask_para]  # indexes of parameters which will be changed [0,npara-1]
dims      = np.sum(mask_para)

# pick only non-masked bounds
lower_bound_mask = lower_bound[np.where(mask_para)]
upper_bound_mask = upper_bound[np.where(mask_para)]
para_dist_mask   = para_dist[np.where(mask_para)]
para_name_mask   = para_name[np.where(mask_para)]

# -------------------------
# read model outputs
# -------------------------
model_output = pickle.load( open( modeloutputs, "rb" ) )

if modeloutputkey == 'All':
    keys = list(model_output.keys())
else:
    keys = [ modeloutputkey ]

model_output = [ np.array(model_output[ikey]) for ikey in keys ]
nkeys = len(model_output)


# -------------------------
# read Morris M
# -------------------------
ff = open(morris_M, "r")
parasets = ff.readlines()
ff.close()
if skip is None:
    skip = int(parasets[0].strip().split(':')[1])
else:
    skip = int(skip)
parasets = parasets[skip:]
for iparaset,paraset in enumerate(parasets):
    parasets[iparaset] = list(map(float,paraset.strip().split()))
parasets = np.array(parasets)

# -------------------------
# read Morris v
# -------------------------
ff = open(morris_v, "r")
parachanged = ff.readlines()
ff.close()
if skip is None:
    skip = int(parachanged[0].strip().split(':')[1])
else:
    skip = int(skip)
parachanged = parachanged[skip:]
for iparachanged,parachan in enumerate(parachanged):
    parachanged[iparachanged] = int(parachan.strip())
parachanged = np.array(parachanged)

# -------------------------
# calculate Elementary Effects
# -------------------------
ee         = np.zeros([dims_all,nkeys],dtype=float)
ee_counter = np.zeros([dims_all,nkeys],dtype=int)
ntraj      = int( np.shape(parasets)[0] / (dims+1) )
nsets      = np.shape(parasets)[0]

for ikey in range(nkeys):
    for iset in range(nsets):

        ipara_changed = parachanged[iset]
        if ipara_changed != -1:
            ee_counter[ipara_changed,ikey] += 1
            if ( len(np.shape(model_output[ikey])) == 1):
                # scalar model output
                ee[ipara_changed,ikey] += np.abs(model_output[ikey][iset]-model_output[ikey][iset+1]) / np.abs(parasets[iset,ipara_changed] - parasets[iset+1,ipara_changed])
            elif ( len(np.shape(model_output[ikey])) == 2):
                # 1D model output
                ee[ipara_changed,ikey] += np.mean(np.abs(model_output[ikey][iset,:]-model_output[ikey][iset+1,:]) / np.abs(parasets[iset,ipara_changed] - parasets[iset+1,ipara_changed]))
            else:
                raise ValueError('Only scalar and 1D model outputs are supported!')

for ikey in range(nkeys):
    for ipara in range(dims_all):
        if ee_counter[ipara,ikey] > 0:
            ee[ipara,ikey] /= ee_counter[ipara,ikey]

# -------------------------
# write final file
# -------------------------
#     format:
#     # model output #1: 'out1'
#     # model output #2: 'out2'
#     # ii     para_name    elemeffect(ii),ii=1:3,jj=1:1      counter(ii),ii=1:3,jj=1:1
#       1      'x_1'        0.53458196335158181               5
#       2      'x_2'        7.0822368906630215                5
#       3      'x_3'        3.5460086652980554                5
f = open(outfile, 'w')
for ikey in range(nkeys):
    f.write('# model output #'+str(ikey+1)+': '+keys[ikey]+'\n')
f.write('# ii     para_name    elemeffect(ii),ii=1:'+str(dims_all)+',jj=1:'+str(nkeys)+'      counter(ii),ii=1:'+str(dims_all)+',jj=1:'+str(nkeys)+' \n')
for ipara in range(dims_all):
    f.write(str(ipara)+'   '+para_name[ipara]+'   '+' '.join(astr(ee[ipara,:],prec=8))+'   '+' '.join(astr(ee_counter[ipara,:]))+'\n')
f.close()
print("wrote:   '"+outfile+"'")
