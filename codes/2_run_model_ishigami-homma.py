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
# An example calling sequence to derive model outputs for previously sampled parameter sets stored
# in an ASCII file (option -i) where some lines might be skipped (option -s). The final model outputs
# are stored in a pickle file (option -o). The model outputs are stored as dictionaries. Multiple
# model outputs are possible.
#
# python 2_run_model_ishigami-homma.py \
#                       -i parameter_sets_1_scaled_para3_M.dat \
#                       -s 8
#                       -o model_output.pkl

"""
Runs a model for a bunch of parameter sets and stores model outputs in a pickle file.

History
-------
Written,  JM, Mar 2019
"""

# -------------------------------------------------------------------------
# Command line arguments - if script
#

# Comment|Uncomment - Begin
#if __name__ == '__main__':

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/lib')

import argparse
import numpy as np
import scipy.stats as stats
import copy
import pickle

infile      = 'example_ishigami-homma/parameter_sets_1_scaled_para3_M.dat'      # name of file containing sampled parameter sets to run the model
outfile     = 'example_ishigami-homma/model_output.pkl'                         # name of file used to save (scalar) model outputs
skip        = None                                                              # number of lines to skip in input file

parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''An example calling sequence to derive model outputs for previously sampled parameter sets stored in an ASCII file (option -i) where some lines might be skipped (option -s). The final model outputs are stored in a pickle file (option -o). The model outputs are stored as dictionaries. Multiple model outputs are possible..''')
parser.add_argument('-i', '--infile', action='store',
                    default=infile, dest='infile', metavar='infile',
                    help="Name of file containing sampled SCALED parameter sets to run the model (default: 'parameter_sets.out').")
parser.add_argument('-s', '--skip', action='store',
                    default=skip, dest='skip', metavar='skip',
                    help="Number of lines to skip in input file (default: None).")
parser.add_argument('-o', '--outfile', action='store',
                    default=outfile, dest='outfile', metavar='outfile',
                    help="Name of file used to save (scalar) model outputs in a pickle file (default: 'model_output.pkl').")

args     = parser.parse_args()
infile   = args.infile
outfile  = args.outfile
skip     = args.skip

del parser, args


def model_function(paraset):
    # function that takes parameter set and returns (scalar) model output
    # here: Ishigami-Homa function (Ishigami and Homma, [1990])
    #            f(x) = sin(p1) + a * sin(p2)**2 + b * p3**4 * sin(p1)
    #       with
    #            a=2.0 and b=1.0
    out = {}
    
    a = 2.0
    b = 1.0
    
    model = np.sin(paraset[0]) + a * np.sin(paraset[1])**2 + b * paraset[2]**4 * np.sin(paraset[0])

    out['out1'] = model                                            # 0D (scalar) model output
    # out['out2'] = [ model / np.float(ii+1) for ii in range(10) ]   # mimicking a 1D model output
    # ...
    # add more model outputs if necessary
    
    return out

# read parameter sets
ff = open(infile, "r")
parasets = ff.readlines()
ff.close()

if skip is None:
    skip = np.int(parasets[0].strip().split(':')[1])
else:
    skip = np.int(skip)
parasets = parasets[skip:]

model_output = {}

for iparaset,paraset in enumerate(parasets):
    paraset = list(map(float,paraset.strip().split()))
    model = model_function(paraset)

    if iparaset == 0:
        for ikey in model.keys():
            model_output[ikey] = []


    for ikey in model.keys():
            
        model_output[ikey].append(model[ikey])
            

pickle.dump( model_output, open( outfile, "wb" ) )

print("wrote:   '"+outfile+"'")
        
