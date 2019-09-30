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


# An example calling sequence to create 10 trajectories for a model with parameters indicated in "parameters.dat"
# to perform the Efficient Elementary Effects (EEE) analysis:
#
# python 1_create_parameter_sets.py -t 10 \
#                                   -d example_ishigami-homma/parameters.dat \
#                                   -o example_ishigami-homma/parameter_sets
#

"""
Sample parameter sets using Morris trajectories. In total T trajectories will be sampled (option -t). 
The parameters can be specified to be either uniform or Gaussian distributed. In case of a uniform 
distributed parameter the lower and upper bound need to be specified. For Gaussian distributed 
parameters, the mean and standard deviation need to be given. Further, it need to be indicated if the 
parameter should be included in the analysis or not. If a parameter is not included, 
the default value will be chosen for each model run. All this need to be set in the parameters input 
file (option -d). The sampled parameter sets and trajectories will be stored in an ASCII file (option -o).

An example parameter setup file that will run the EEE analysis for three parameters, looks like:
     
     # para   dist       lower     upper     default   informative(0)_or_noninformative(1)
     #                   mean      stddev
     x_1      uniform    -3.14159  3.14159   0.00000   0
     x_2      uniform    -3.14159  3.14159   0.00000   0
     x_3      uniform    -3.14159  3.14159   0.00000   0


History
-------
Written,  JM, Mar 2019
"""

# -------------------------------------------------------------------------
# Command line arguments
#
outfile   = 'example_ishigami-homma/parameter_sets'
maskfile  = 'example_ishigami-homma/parameters.dat'
ntraj     = 5
nfiles    = 1

import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="Sample parameter sets using Morris trajectories. In total T trajectories will be sampled (option -t). The parameters can be specified to be either uniform or Gaussian distributed. In case of a uniform distributed parameter the lower and upper bound need to be specified. For Gaussian distributed parameters, the mean and standard deviation need to be given. Further, it need to be indicated if the parameter should be included in the analysis or not. If a parameter is not included, the default value will be chosen for each model run. All this need to be set in the parameters input file (option -d). The sampled parameter sets and trajectories will be stored in an ASCII file (option -o).")

parser.add_option('-o', '--outfile', action='store', dest='outfile', type='string',
                  default=outfile, metavar='File',
                  help='Name of pdf output file descriptor as [outfile]_[file_id]_para[npara].dat (default: outfile=parameter_sets).')
parser.add_option('-d', '--maskfile', action='store', dest='maskfile', type='string',
                  default=maskfile, metavar='File',
                  help='Name of file where all model parameters are specified including their distribution, distribution parameters, default value and if included in analysis or not. (default: maskfile=parameters.dat).')
parser.add_option('-n', '--nfiles', action='store', dest='nfiles', type='int',
                  default=nfiles, metavar='Amount',
                  help='Number of files to be sampled containing ntraj each. (default: nfiles=1).')
parser.add_option('-t', '--ntraj', action='store', dest='ntraj', type='int',
                  default=ntraj, metavar='Amount',
                  help='Number of trajectories sampled. (default: ntraj=5).')
(opts, args) = parser.parse_args()

outfile   = opts.outfile   # morris_database
maskfile  = opts.maskfile  # parameters.dat
ntraj     = opts.ntraj     # 100
nfiles    = opts.nfiles    # 1

# print('outfile   :: '+outfile)
# print('maskfile  :: '+maskfile)
# print('ntraj     :: '+str(ntraj))
# print('nfiles    :: '+str(nfiles))

del parser, opts, args


# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/lib')

import numpy       as np
import scipy.stats as stats
import copy
from fsread          import fsread              # in lib/
from autostring      import astr                # in lib/
from morris          import morris_sampling     # in lib/


# maskfile has following header:
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

fileID = astr(np.arange(1,nfiles+1),zero=True)

# print('dims = '+astr(dims))

for kk in range(nfiles):
    # print('Sampling #'+astr(kk+1)+' of '+astr(nfiles))

    lower_bound_01 = np.zeros(dims)
    upper_bound_01 = np.ones(dims)

    # samples N=ntraj*10 trajectories and then picks the r=ntraj ones that are most appart from each other
    [OptMatrix, OptOutVec] = morris_sampling(dims, lower_bound_01, upper_bound_01, N=ntraj*10, p=6, r=ntraj, Diagnostic=0)

    # Write M = OptMatrix (unscaled)
    outfile_name = outfile+'_'+fileID[kk]+'_para'+astr(dims)+'_M.dat'
    f = open(outfile_name, 'w')
    f.write('header lines: '+astr(dims + 5)+'\n')
    f.write('Morris sequences generated using PYTHON \n')
    f.write('first: reference set, then: '+astr(dims)+' lines OAT sample \n')
    f.write('parameter ranges are: \n')

    for ii in range(dims):
        if para_dist_mask[ii].lower() == 'uniform':
            f.write('      p_{'+astr((np.where(mask_para))[0][ii]+1)+'} = "'+para_name_mask[ii]+'" = Uniform[ '+astr(lower_bound_mask[ii],4)+' , '+astr(upper_bound_mask[ii],4)+' ]\n')
        elif para_dist_mask[ii].lower() == 'gaussian':
            f.write('      p_{'+astr((np.where(mask_para))[0][ii]+1)+'} = "'+para_name_mask[ii]+'" = Gaussian[ '+astr(lower_bound_mask[ii],4)+' , '+astr(upper_bound_mask[ii],4)+' ]\n')
        else:
            raise ValueError('This distribution is not implemented. Only ["uniform","gaussian"].')
    f.write( ' \n')

    tmpOptMatrix = copy.deepcopy(OptMatrix)
    for ii in range(dims):
        if para_dist_mask[ii].lower() == 'gaussian':
            tmpOptMatrix[:,ii]  = 0.01 + 0.98 * tmpOptMatrix[:,ii]
            tmpOptMatrix[:,ii] = stats.norm.ppf(tmpOptMatrix[:,ii], loc=0.0, scale=1.0)
        
    for ii in range(ntraj*(dims+1)):
        spread_set=np.copy(initial)
        spread_set[idx_para]=tmpOptMatrix[ii,:]
        f.write(' '.join(astr(spread_set,8)))
        f.write('\n')

    f.close()
    print("wrote:   '"+outfile_name+"'")

    # morris_sampling actually don't support 'Gaussian' distribution...
    for ii in range(dims):
        if para_dist_mask[ii].lower() == 'gaussian':
            
            # N[0,1]; cut at 1% and 99% percentiles
            OptMatrix[:,ii]  = 0.01 + 0.98 * OptMatrix[:,ii]
            OptMatrix[:,ii] = stats.norm.ppf(OptMatrix[:,ii], loc=lower_bound_mask[ii], scale=upper_bound_mask[ii])

        elif para_dist_mask[ii].lower() == 'uniform':

            # uniform scaling to [a,b]
            OptMatrix[:,ii] *= (upper_bound_mask[ii] - lower_bound_mask[ii])
            OptMatrix[:,ii] += lower_bound_mask[ii]

        else:

            raise ValueError('This distribution is not implemented. Only ["uniform","gaussian"].')

    # Write M = OptMatrix (scaled)
    outfile_name = outfile+'_'+fileID[kk]+'_scaled_para'+astr(dims)+'_M.dat'
    f = open(outfile_name, 'w')
    f.write('header lines: '+astr(dims + 5)+'\n')
    f.write('Morris sequences generated using PYTHON \n')
    f.write('first: reference set, then: '+astr(dims)+' lines OAT sample \n')
    f.write('parameter ranges are: \n')

    for ii in range(dims):
        if para_dist_mask[ii].lower() == 'uniform':
            f.write('      p_{'+astr((np.where(mask_para))[0][ii]+1)+'} = "'+para_name_mask[ii]+'" = Uniform[ '+astr(lower_bound_mask[ii],4)+' , '+astr(upper_bound_mask[ii],4)+' ]\n')
        elif para_dist_mask[ii].lower() == 'gaussian':
            f.write('      p_{'+astr((np.where(mask_para))[0][ii]+1)+'} = "'+para_name_mask[ii]+'" = Gaussian[ '+astr(lower_bound_mask[ii],4)+' , '+astr(upper_bound_mask[ii],4)+' ]\n')
        else:
            raise ValueError('This distribution is not implemented. Only ["uniform","gaussian"].')
    f.write( ' \n')

    for ii in range(ntraj*(dims+1)):
        spread_set=np.copy(initial)
        spread_set[idx_para]=OptMatrix[ii,:]
        f.write(' '.join(astr(spread_set,8)))
        f.write('\n')

    f.close()
    print("wrote:   '"+outfile_name+"'")

    # Write v = OptOutVec
    outfile_name = outfile+'_'+fileID[kk]+'_para'+astr(dims)+'_v.dat'
    f = open(outfile_name, 'w')
    f.write('header lines: '+astr(dims + 5)+'\n')
    f.write('Morris sequences generated using PYTHON \n')
    f.write('first: reference set, then: '+astr(dims)+' lines OAT sample \n')
    f.write('parameter ranges are: \n')

    for ii in range(dims):
        if para_dist_mask[ii].lower() == 'uniform':
            f.write('      p_{'+astr((np.where(mask_para))[0][ii]+1)+'} = "'+para_name_mask[ii]+'" = Uniform[ '+astr(lower_bound_mask[ii],4)+' , '+astr(upper_bound_mask[ii],4)+' ]\n')
        elif para_dist_mask[ii].lower() == 'gaussian':
            f.write('      p_{'+astr((np.where(mask_para))[0][ii]+1)+'} = "'+para_name_mask[ii]+'" = Gaussian[ '+astr(lower_bound_mask[ii],4)+' , '+astr(upper_bound_mask[ii],4)+' ]\n')
        else:
            raise ValueError('This distribution is not implemented. Only ["uniform","gaussian"].')
    f.write( ' \n')

    for ii in range(ntraj*(dims+1)):
        if (OptOutVec[ii] == -1.0):
            f.write(astr(-1.0))
        else:
            f.write(astr(1.0*idx_para[np.int(OptOutVec[ii])]))
        f.write('\n')

    f.close()
    print("wrote:   '"+outfile_name+"'")












    










