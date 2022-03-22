#!/usr/bin/env python
from __future__ import print_function

# Copyright 2019 Juliane Mai - juliane.mai(at)uwaterloo.ca
# Adapted 2022 by Simon Lin - simon.lin(at)uwaterloo.ca
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
# python 2_run_model_raven-hmets.py \
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
sys.path.append(os.path.abspath(dir_path+'/../examples/robin/model'))

import argparse
import numpy as np
import scipy.stats as stats
import copy
import pickle
from   pathlib2        import Path
import subprocess
import shutil

from   raven_model_files import RVI, RVT, RVP, RVP_CHANNEL, RVH, RVH_LAKE, RVC        # in examples/model/robin; adapted from examples/raven-hmets/model
from   robin_model_files import PAR, HRUCROP, INFO, INIC, MODEL, MGT, OBS             # in examples/model/robin; adapted from examples/raven-hmets/model/
from   raven_common      import writeString, makeDirectories                          # for modifying model input files; copied from examples/raven-hmets/model/
from   fread             import fread                                                 # in lib/

infile      = 'examples/robin/parameter_sets_1_scaled_para15_M.dat'                   # name of file containing sampled parameter sets to run the model
outfile     = 'examples/robin/model_output.pkl'                                       # name of file used to save (scalar) model outputs
skip        = None                                                                    # number of lines to skip in input file

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
    dict_dparas['sum_x43_x45']    = paras[43]+paras[45]          # kc_max for Forest type
    dict_dparas['sum_x44_x46']    = paras[44]+paras[46]          # kc_max for Forest2 type
    dict_dparas['sum_x47_x49']    = paras[47]+paras[49]          # laimax for Forest type
    dict_dparas['sum_x48_x50']    = paras[48]+paras[50]          # laimax for Forest2 type
    # ---------------
    # paste all paras into template files
    # ---------------
    #   ex.:    string = "parameter x01 = {par[x01]} and another parameter x02 = {par[x02]}"
    #           keys   = ['x01','x02']
    #           vals   = [1.0,3.0]
    #           string.format(par=dict(zip(keys,vals)))
    #
    #           --> 'parameter x01 = 1.0 and another parameter x02 = 3.0'
    #
    # to replace patterns: {par[x01]} by parameter value paras[0]
    #                      {par[x02]} by parameter value paras[1]
    #                      ...
    if len(paras) > 9 and len(paras) < 100:
        keys_paras = ["x{:02d}".format(ii)   for ii in range(1,len(paras)+1) ]
    elif len(paras) > 99 and len(paras) < 1000:
        keys_paras = ["x{:03d}".format(ii)    for ii in range(1,len(paras)+1) ]
    elif len(paras) <= 9:
        keys_paras = ["x{:01d}".format(ii) for ii in range(1,len(paras)+1) ]
    else:
        raise ValueError("More than 999 parameters are not implemented yet!")
    vals_paras = paras
    dict_paras = dict(zip(keys_paras,vals_paras))

    # fill in to templates
    # templates need to have patterns:
    #         {par[x01]},  {par[x02]},                     ... for parameters
    #         {dpar[something]},  {dpar[somethingelse]},   ... for derived parameters

    # ---------------
    # create a run folder
    # ---------------
    tmp_folder = os.path.abspath(dir_path+"/../tmp/eee-analysis/Turkey_Lakes/"+str(run_id)) # "/tmp/juletest" #  TODO a generic folder name in /tmp
    robin_exe_name   = os.path.abspath(dir_path+"/../examples/robin/model/raven_robin")
    raven_obs_folder = os.path.abspath(dir_path+"/../examples/robin/model/obs")
    raven_forcing_folder = os.path.abspath(dir_path+"/../examples/robin/model/Forcing")
    robin_obs_folder = os.path.abspath(dir_path+"/../examples/robin/model/cropmodel/obs")
    print(robin_obs_folder)
 
    if os.path.exists(tmp_folder):
        shutil.rmtree(tmp_folder)

    # all RAVEN setup files (ASCII files can only be read once at a time, creating local copies speeds up run time)
    writeString( Path(tmp_folder,"Turkey_Lake.rvi"), RVI.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"Turkey_Lake.rvp"), RVP.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"Turkey_Lake.rvh"), RVH.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"Turkey_Lake.rvt"), RVT.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"Turkey_Lake.rvc"), RVC.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"Turkey_Lake_Lake.rvh"), RVH_LAKE.format(par=dict_paras,dpar=dict_dparas) )   
    writeString( Path(tmp_folder,"Turkey_Lake_channel.rvp"), RVP_CHANNEL.format(par=dict_paras,dpar=dict_dparas) )   

    # all ROBIN setup files (ASCII files can only be read once at a time, creating local copies speeds up run time)
    writeString( Path(tmp_folder,"cropmodel/crop.hrucrop"), HRUCROP.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"cropmodel/crop.info"), INFO.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"cropmodel/crop.inic"), INIC.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"cropmodel/crop.mgt"), MGT.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"cropmodel/crop.obs"), OBS.format(par=dict_paras,dpar=dict_dparas) )
    writeString( Path(tmp_folder,"cropmodel/crop.model"), MODEL.format(par=dict_paras,dpar=dict_dparas) )
   
   # write sampled Robin parameter set to crop.par file for use
    writeString( Path(tmp_folder,"cropmodel/crop.par"), PAR.format(par=dict_paras,dpar=dict_dparas) )

    # link executable
    if not(os.path.exists(str(Path(tmp_folder,os.path.basename(robin_exe_name))))):
        print("from: ",os.path.realpath(robin_exe_name))
        print("to:   ",str(Path(tmp_folder,os.path.basename(robin_exe_name))))
        os.symlink(os.path.realpath(robin_exe_name), str(Path(tmp_folder,os.path.basename(robin_exe_name))))

    # link observations folders (raven)
    if not(os.path.exists(str(Path(tmp_folder,os.path.basename(raven_obs_folder))))):
        os.symlink(os.path.realpath(raven_obs_folder), str(Path(tmp_folder,os.path.basename(raven_obs_folder))))

    # link forcing folders (raven)
    if not(os.path.exists(str(Path(tmp_folder,os.path.basename(raven_forcing_folder))))):
        os.symlink(os.path.realpath(raven_forcing_folder), str(Path(tmp_folder,os.path.basename(raven_forcing_folder))))

    # link observations folder in coupled model subdirectory "cropmodel" (robin)
    if not(os.path.exists(str(Path(tmp_folder,"cropmodel",os.path.basename(robin_obs_folder))))):
        os.symlink(os.path.realpath(robin_obs_folder), str(Path(tmp_folder,"cropmodel",os.path.basename(robin_obs_folder))))

    print(str(Path(tmp_folder,"cropmodel",os.path.basename(robin_obs_folder))))
 
    # create ouput folder
    out_folder = str(Path(tmp_folder,"output"))
    os.makedirs(out_folder)
    
    # ---------------
    # run the model with these input rv* files
    # ---------------        
    cmd = [str(Path(tmp_folder,os.path.basename(robin_exe_name))),str(Path(tmp_folder,"Turkey_Lake")),"-o",str(Path(tmp_folder,"output"))+'/']
    #cmd = ["./"+str(os.path.basename(robin_exe_name)),str("Turkey_Lake"),"-o",str("output")+'/']
    #cmd_up = str("cd '"+tmp_folder+"'")
    #cmd_down = str("cd '"+dir_path+"'")

    print("run cmd: ",' '.join(cmd))

    # Robin literally searches for all the input files based on the current working directory ("./cropmodel/")
    #process_up = subprocess.call(cmd_up, shell=True)
    #print(os.getcwd())
    process = subprocess.Popen(cmd, cwd = tmp_folder, stdout=subprocess.PIPE)
    #process_down = subprocess.call(cmd_down, shell=True)

    print("")
    print("Raven standard output:")
    for line in process.stdout:
        print(">>> ",line.rstrip()) # rstrip removes trailing \n

    if not(os.path.exists(str(Path(tmp_folder,"output","Diagnostics.csv")))):            
        print("")
        print("ERROR: No Diagnostics.csv produced")
        print("")
        print("Raven error file content:")
        ff = open(str(Path(tmp_folder,"output","Raven_errors.txt")), "r")
        lines = ff.readlines()
        ff.close()
        for line in lines:
            print(">>> ",line.rstrip()) # rstrip removes trailing \n

        raise ValueError("ERROR: No Diagnostics.csv produced (scroll up to see content of error file)")

    model = {}

    # ---------------
    # extract model output: Diagnostic KGE scores for all vegetation observations
    # ---------------
    ff = open(str(Path(tmp_folder,"cropout","Evaluation_Metric.csv")), "r")
    lines = ff.readlines()
    ff.close()
    
    # probably a more elegant way to do this but I'm bad at python...
    model['KGE_LAI_C31'] = np.float(lines[1].strip().split(',')[4])             # first index is row number, second index is column number
    model['KGE_LAI_C32'] = np.float(lines[2].strip().split(',')[4])
    model['KGE_stocking_C31'] = np.float(lines[3].strip().split(',')[4])
    model['KGE_stocking_C32'] = np.float(lines[4].strip().split(',')[4])
    model['KGE_W_stem_C31'] = np.float(lines[5].strip().split(',')[4])
    model['KGE_W_stem_C32'] = np.float(lines[6].strip().split(',')[4])
    model['KGE_basalarea_C31'] = np.float(lines[7].strip().split(',')[4])
    model['KGE_basalarea_C32'] = np.float(lines[8].strip().split(',')[4])
    model['KGE_volume_C31'] = np.float(lines[9].strip().split(',')[4])
    model['KGE_volume_C32'] = np.float(lines[10].strip().split(',')[4])
    model['KGE_W_stem_tree_C31'] = np.float(lines[11].strip().split(',')[4])
    model['KGE_W_stem_tree_C32'] = np.float(lines[12].strip().split(',')[4])
    model['KGE_W_Leaf_C31'] = np.float(lines[13].strip().split(',')[4])
    model['KGE_W_Leaf_C32'] = np.float(lines[14].strip().split(',')[4])

    print("KGE_LAI_C31:            ",model['KGE_LAI_C31'])
    print("")
    print("KGE_LAI_C32:            ",model['KGE_LAI_C32'])
    print("")
    print("KGE_stocking_C31:            ",model['KGE_stocking_C31'])
    print("")
    print("KGE_stocking_C32:            ",model['KGE_stocking_C32'])
    print("")
    print("KGE_W_stem_C31:            ",model['KGE_W_stem_C31'])
    print("")
    print("KGE_W_stem_C32:            ",model['KGE_W_stem_C32'])
    print("")
    print("KGE_basalarea_C31:           ",model['KGE_basalarea_C31'])
    print("")
    print("KGE_basalarea_C32:           ",model['KGE_basalarea_C32'])
    print("")
    print("KGE_volume_C31:           ",model['KGE_volume_C31'])
    print("")
    print("KGE_volume_C32:           ",model['KGE_volume_C32'])
    print("")
    print("KGE_W_stem_tree_C31:           ",model['KGE_W_stem_tree_C31'])
    print("")
    print("KGE_W_stem_tree_C32:           ",model['KGE_W_stem_tree_C32'])
    print("")
    print("KGE_W_Leaf_C31:           ",model['KGE_W_Leaf_C31'])
    print("")
    print("KGE_W_Leaf_C31:           ",model['KGE_W_Leaf_C32'])

    ## ---------------
    ## extract model output: Hydrographs: simulated Q
    ## ---------------
    #model['Q']  = 0.0
    #warmup = 2*365  # 1 # model timestep 1 day and want to skip 2 years  # first day 1991-01-01 00:00:00.00 (checked)
    #model['Q']  = np.transpose(fread(str(Path(tmp_folder,"output","Hydrographs.csv")),skip=warmup+1,cskip=4,nc=1))[0]

    #print("Q:              ",model['Q'][0:4],"...",model['Q'][-4:])
    #print("Q_range:         [",np.min(model['Q']),",",np.max(model['Q']),"]")
    #print("shape Q:        ",np.shape(model['Q']))
    #print("")

    ## ---------------
    ## extract model output: BETWEEN_PONDED_WATER_AND_SOIL[0]_Daily_Average_BySubbasin.csv: accumulated infiltration volume
    ## ---------------
    #model['infiltration']  = 0.0
    #warmup = 2*365  # 1 # model timestep 1 day and want to skip 2 years  # first day 1990-12-31 00:00:00.00 (checked) But all timesteps are shifted by 1 day...
    ##
    ## de-accumulated infiltration volume
    #model['infiltration'] = np.transpose(fread(str(Path(tmp_folder,"output","BETWEEN_PONDED_WATER_AND_SOIL[0]_Daily_Average_BySubbasin.csv")),skip=warmup,cskip=2,nc=1))[0]
    #model['infiltration'] = np.diff(model['infiltration'])

    #print("Infiltration I: ",model['infiltration'][0:4],"...",model['infiltration'][-4:])
    #print("I_range:         [",np.min(model['infiltration']),",",np.max(model['infiltration']),"]")
    #print("shape I:        ",np.shape(model['infiltration']))
    #print("")

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
