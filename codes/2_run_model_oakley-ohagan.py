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
# python 2_run_model_oakley-ohagan.py \
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

infile      = 'example_oakley-ohagan/parameter_sets_1_scaled_para15_M.dat'     # name of file containing sampled parameter sets to run the model
outfile     = 'example_oakley-ohagan/model_output.pkl'                         # name of file used to save (scalar) model outputs
skip        = None                                                             # number of lines to skip in input file

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
    # here: Oakley-O'Hagan function (Oakley & O'Hagan, [2004])
    out = {}
    
    paraset = np.array(paraset)

    isone = True
    iparaset = paraset[:,np.newaxis]

    nn = 15

    assert iparaset.shape[0] == nn, 'paraset.shape[0] must '+str(nn)+'.'

    a1 = np.array([0.01, 0.05, 0.23, 0.04, 0.12, 0.39, 0.39, 0.61, 0.62, 0.40, 1.07, 1.15, 0.79, 1.12, 1.20])
    a2 = np.array([0.43, 0.09, 0.05, 0.32, 0.15, 1.04, 0.99, 0.97, 0.90, 0.81, 1.84, 2.47, 2.39, 2.00, 2.26])
    a3 = np.array([0.10, 0.21, 0.08, 0.27, 0.13, 0.75, 0.86, 1.03, 0.84, 0.80, 2.21, 2.04, 2.40, 2.05, 1.98])
    M  = np.array([
        [-0.02, -0.19,   0.13,  0.37,  0.17,  0.14, -0.44,  -0.08,  0.71, -0.44,  0.5,  -0.02, -0.05,  0.22,  0.06],
        [ 0.26,  0.05,   0.26,  0.24, -0.59, -0.08, -0.29,   0.42,  0.5,   0.08, -0.11,  0.03, -0.14, -0.03, -0.22],
        [-0.06,  0.2,    0.1,  -0.29, -0.14,  0.22,  0.15,   0.29,  0.23, -0.32, -0.29, -0.21,  0.43,  0.02,  0.04],
        [ 0.66,  0.43,   0.3,  -0.16, -0.31, -0.39,  0.18,   0.06,  0.17,  0.13, -0.35,  0.25, -0.02,  0.36, -0.33],
        [-0.12,  0.12,   0.11,  0.05, -0.22,  0.19, -0.07,   0.02, -0.1,   0.19,  0.33,  0.31, -0.08, -0.25,  0.37],
        [-0.28, -0.33,  -0.1,  -0.22, -0.14, -0.14, -0.12,   0.22, -0.03, -0.52,  0.02,  0.04,  0.36,  0.31,  0.05],
        [-0.08,  0.004,  0.89, -0.27, -0.08, -0.04, -0.19,  -0.36, -0.17,  0.09,  0.4,  -0.06,  0.14,  0.21, -0.01],
        [-0.09,  0.59,   0.03, -0.03, -0.24, -0.1,   0.03,   0.1,  -0.34,  0.01, -0.61,  0.08,  0.89,  0.14,  0.15],
        [-0.13,  0.53,   0.13,  0.05,  0.58,  0.37,  0.11,  -0.29, -0.57,  0.46, -0.09,  0.14, -0.39, -0.45, -0.15],
        [ 0.06, -0.32,   0.09,  0.07, -0.57,  0.53,  0.24,  -0.01,  0.07,  0.08, -0.13,  0.23,  0.14, -0.45, -0.56],
        [ 0.66,  0.35,   0.14,  0.52, -0.28, -0.16, -0.07,  -0.2,   0.07,  0.23, -0.04, -0.16,  0.22,  0,    -0.09],
        [ 0.32, -0.03,   0.13,  0.13,  0.05, -0.17,  0.18,   0.06, -0.18, -0.31, -0.25,  0.03, -0.43, -0.62, -0.03],
        [-0.29,  0.03,   0.03, -0.12,  0.03, -0.34, -0.41,   0.05, -0.27, -0.03,  0.41,  0.27,  0.16, -0.19,  0.02],
        [-0.24, -0.44,   0.01,  0.25,  0.07,  0.25,  0.17,   0.01,  0.25, -0.15, -0.08,  0.37, -0.3,   0.11, -0.76],
        [ 0.04, -0.26,   0.46, -0.36, -0.95, -0.17,  0.003,  0.05,  0.23,  0.38,  0.46, -0.19,  0.01,  0.17,  0.16] ])

    y = np.dot(a1,iparaset) + np.dot(a2,np.sin(iparaset)) + np.dot(a3,np.cos(iparaset))
    for i in range(iparaset.shape[1]):
        y[i] += np.dot(iparaset[:,i].T,np.dot(M,iparaset[:,i]))
    model = y[0]

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
        
