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
# python 2_run_model_cequeau.py \
#                       -i parameter_sets_1_scaled_para21_M.dat \
#                       -s XXX
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
sys.path.append(os.path.abspath(dir_path+'/lib'))
sys.path.append(os.path.abspath(dir_path+'/../examples/cequeau-nc/model'))

import argparse
import numpy as np
import scipy.stats as stats
import copy
import pickle
from   pathlib2        import Path
import subprocess
import shutil
import datetime

from   cequeau_templates import EXECUTION_XML, PARAMETRES_XML, BASSINVERSANT_XML # in examples/cequeau-nc/model/
from   cequeau_common    import writeString, makeDirectories, get_discharge      # in examples/cequeau-nc/model/
from   fread             import fread                                            # in lib/

infile      = 'example_cequeau-nc/parameter_sets_1_scaled_para9_M.dat'     # name of file containing sampled parameter sets to run the model
outfile     = 'example_cequeau-nc/model_output.pkl'                        # name of file used to save (scalar) model outputs
skip        = None                                                         # number of lines to skip in input file

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

def model_function(paras, run_id=None):
    # input:
    #     paras     ... list of model parameters scaled to their range;
    #                   values for all N model parameters have to be provided
    #                   example:
    #                        [ x1, x2, x3, x4, .... ]
    #     run_id    ... optional name of this run (to, e.g., print or store in a file)
    #                   example:
    #                        run_aset_001
    # output:
    #     model output in dictionary
    #     example:
    #          model['out'] = 7.4

    if not(run_id is None):
        print("Run ID: ",run_id)

    # ---------------
    # derive some parameters
    # ---------------
    dict_dparas = {}
    
    # dict_dparas['sum_x05_x06']  = paras[4]+paras[5]           # MAX_MELT_FACTOR > MIN_MELT_FACTOR
    # dict_dparas['sum_x09_x10']  = paras[8]+paras[9]           # SNOW_SWI_MAX > SNOW_SWI_MIN
    # dict_dparas['half_x20']     = paras[19] * 0.5 * 1000      # half the value but in [mm] not [m]
    # dict_dparas['half_x21']     = paras[20] * 0.5 * 1000      # half the value but in [mm] not [m]
    
    # ---------------
    # paste all paras into template files
    # ---------------
    #   ex.:    string = "parameter v01 = {par[v01]} and another parameter v02 = {par[v02]}"
    #           keys   = ['v01','v02']
    #           vals   = [1.0,3.0]
    #           string.format(par=dict(zip(keys,vals)))
    #
    #           --> 'parameter v01 = 1.0 and another parameter v02 = 3.0'
    #
    # to replace patterns: {par[v01]} by parameter value paras[0]
    #                      {par[v02]} by parameter value paras[1]
    #                      ...
    if len(paras) > 9 and len(paras) < 100:
        keys_paras = ["v{:02d}".format(ii)   for ii in range(1,len(paras)+1) ]
    elif len(paras) > 99 and len(paras) < 1000:
        keys_paras = ["v{:03d}"   for ii in range(1,len(paras)+1) ]
    elif len(paras) <= 9:
        keys_paras = ["v"+str(ii) for ii in range(1,len(paras)+1) ]
    else:
        raise ValueError("More than 999 parameters are not implemented yet!")
    vals_paras = paras
    dict_paras = dict(zip(keys_paras,vals_paras))

    # fill in to templates
    # templates need to have patterns:
    #         {par[v01]},  {par[v02]},                     ... for parameters
    #         {dpar[something]},  {dpar[somethingelse]},   ... for derived parameters

    # ---------------
    # create a run folder
    # ---------------
    tmp_folder = "/tmp/eee-analysis/"+str(run_id) # "/tmp/juletest" #  TODO a generic folder name in /tmp
    cequeau_exe_name    = os.path.abspath(dir_path+"/../"+"examples/cequeau-nc/model/cequeau")
    cequeau_obs_folder  = os.path.abspath(dir_path+"/../"+"examples/cequeau-nc/model/data_obs")

    if os.path.exists(tmp_folder):
        shutil.rmtree(tmp_folder)

    print("dict_paras: ",dict_paras)

    # all CEQUEAU setup files
    writeString( Path(tmp_folder,"execution.xml"),     EXECUTION_XML.format(path=tmp_folder,par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"parametres.xml"),    PARAMETRES_XML.format(path=tmp_folder,par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"bassinVersant.xml"), BASSINVERSANT_XML.format(path=tmp_folder,par=dict_paras,dpar=dict_dparas) )

    # link executable
    if not(os.path.exists(str(Path(tmp_folder,os.path.basename(cequeau_exe_name))))):
        print("from: ",os.path.realpath(cequeau_exe_name))
        print("to:   ",str(Path(tmp_folder,os.path.basename(cequeau_exe_name))))
        os.symlink(os.path.realpath(cequeau_exe_name), str(Path(tmp_folder,os.path.basename(cequeau_exe_name))))

    # link observations folder
    if not(os.path.exists(str(Path(tmp_folder,os.path.basename(cequeau_obs_folder))))):
        os.symlink(os.path.realpath(cequeau_obs_folder), str(Path(tmp_folder,os.path.basename(cequeau_obs_folder))))

    # create ouput folder
    out_folder = str(Path(tmp_folder,"output"))
    os.makedirs(out_folder)
    
    # ---------------
    # run the model with these input rv* files
    # ---------------        
    cmd = [str(Path(tmp_folder,os.path.basename(cequeau_exe_name))),str(Path(tmp_folder,"execution.xml"))]
    print("run cmd: ",' '.join(cmd))

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    print("")
    print("Cequeau standard output:")
    for line in process.stdout:
        print(">>> ",line.rstrip()) # rstrip removes trailing \n

    if not(os.path.exists(str(Path(tmp_folder,"output","resultats.nc")))):            
        print("")
        print("ERROR: No resultats.nc produced")
        print("")
        print("Cequeau error file content:")
        # ff = open(str(Path(tmp_folder,"output","Raven_errors.txt")), "r")
        # lines = ff.readlines()
        # ff.close()
        # for line in lines:
        #     print(">>> ",line.rstrip()) # rstrip removes trailing \n

        raise ValueError("ERROR: No resultats.nc produced (scroll up to see content of error file)")

    model = {}

    # # ---------------
    # # extract model output: Diagnostics: NSE
    # # ---------------
    # model['nse'] = 0.0
    # ff = open(str(Path(tmp_folder,"output","Diagnostics.csv")), "r")
    # lines = ff.readlines()
    # ff.close()

    # nse = np.float(lines[-1].strip().split(',')[2])
    # print("NSE:            ",nse)
    # model['nse'] = nse
    # print("")
    
    # ---------------
    # extract model output: Hydrographs: simulated Q
    # ---------------
    model['Q']  = 0.0
    back_days   = 0
    # 1	PD	      Passes Dangereuses
    # 2	PERIB	  Peribonka	
    # 3	LM	      Lac Manouane	
    # 4	MBLANC	  Montagnes Blanches
    basin_id    = 1    # starts with 1
    start_day   = datetime.datetime(2000,01,01,0,0) - datetime.timedelta(days=back_days) # TODO
    qsim_file   = str(Path(tmp_folder,"output","resultats.nc"))
    qsim        = get_discharge(start_day, qsim_file, 7, var="debitExutoire", ibasin=basin_id-1, group="etatsCP", ilag=0)
    model['Q']  = qsim

    print("Q:              ",model['Q'][0:4],"...",model['Q'][-4:])
    print("Q_range:         [",np.min(model['Q']),",",np.max(model['Q']),"]")
    print("shape Q:        ",np.shape(model['Q']))
    print("")

    # ---------------
    # cleanup
    # ---------------
    if os.path.exists(tmp_folder):
        shutil.rmtree(tmp_folder)

    return model

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

# this loop could be easily parallized and modified such that it
# actually submits multiple tasks to a HPC
for iparaset,paraset in enumerate(parasets):
    
    paraset = list(map(float,paraset.strip().split()))
    model = model_function(paraset,run_id='run_set_'+str(iparaset))

    if iparaset == 0:
        for ikey in model.keys():
            model_output[ikey] = []


    for ikey in model.keys():
            
        model_output[ikey].append(model[ikey])
            

pickle.dump( model_output, open( outfile, "wb" ) )

print("wrote:   '"+outfile+"'")
        
