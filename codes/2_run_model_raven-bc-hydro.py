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
# python 2_run_model_raven-hmets.py \
#                       -i parameter_sets_1_scaled_para21_M.dat \
#                       -s XXX
#                       -o model_output.pkl
#                       -p Alouette

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
sys.path.append(os.path.abspath(dir_path+'/../examples/raven-bc-hydro/model'))

import argparse
import numpy as np
import scipy.stats as stats
import copy
import pickle
from   pathlib2        import Path
import subprocess
import shutil
import glob as glob

#from   raven_templates import RVI, RVT, RVP, RVH, RVC          # in examples/raven-bc-hydro/model/
#from   raven_common    import writeString, makeDirectories     # in examples/raven-bc-hydro/model/
from   fread           import fread                            # in lib/

infile      = 'example_raven-bc-hydro/parameter_sets_1_scaled_para15_M.dat'     # name of file containing sampled parameter sets to run the model
outfile     = 'example_raven-bc-hydro/model_output.pkl'                         # name of file used to save (scalar) model outputs
skip        = None                                                              # number of lines to skip in input file
subproject  = '20250404-v5-Alouette'

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
parser.add_argument('-p', '--subproject', action='store',
                    default=subproject, dest='subproject', metavar='subproject',
                    help="Name of subproject (e.g. subbasin) (default: 'Alouette').")

args       = parser.parse_args()
infile     = args.infile
outfile    = args.outfile
skip       = args.skip
subproject = args.subproject


print(subproject.split('-'))
setupdate = subproject.split('-')[0]  # '20250404'
version   = subproject.split('-')[1]  #'v5'
basin     = subproject.split('-')[2]  # 'Alouette'


del parser, args

def writeString(fname, string):
    """
    Arguments
    ---------
    fname (Union[Text, pathlib2.Path]) : file name
    string (Text)                      : file's content to write

    Return
    ------
    None

    Purpose
    -------
    Write given string to file. All necessary directories
    will be created.
    """

    makeDirectories([Path(fname)])
    with open(str(fname), "w") as f:
        f.write(string)
    f.close()

def makeDirectories(fnames):
    """
    Arguments
    ---------
    fnames (List[pathlib2.Path]): file or directory names.

    Return
    ------
    None


    Purpose
    -------
    Create all directories necessary to be able to access
    the given filenames or directory
    
    Note
    ----
    Anything in 'fnames' with a file extension in the form of '.*' is
    considered to be file, anything else to represent a directory.
    """

    dirs = set([f.parent if f.suffix else f for f in fnames])
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)

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

    dict_paras = {}

    if version in ['v5','v6']:
        dict_paras['dar_x02'] = paras[0] + paras[1]
        dict_paras['dar_x05'] = 10.0**paras[4]
        dict_paras['dar_x06'] = 10.0**paras[5]
        dict_paras['dar_x12'] = paras[11] + paras[14]
    elif version in ['v14','v15']:
        dict_paras['dar_x02'] = paras[0] + paras[1]
        dict_paras['dar_x05'] = 10.0**paras[3]
        dict_paras['dar_x06'] = 10.0**paras[4]
        dict_paras['dar_x12'] = paras[9] + paras[12]
    else:
        raise ValueError('Version not implemented')

    #print(paras)

    # ---------------
    # paste all paras
    # ---------------
    #   ex.:    string = "parameter x01 = {par[x01]} and another parameter x02 = {par[x02]}"
    #           keys   = ['x01','x02']
    #           vals   = [1.0,3.0]
    #           string.format(par=dict(zip(keys,vals)))
    #
    #           --> 'parameter x01 = 1.0 and another parameter x02 = 3.0'
    #
    # to replace patterns: par_x01 by parameter value paras[0]
    #                      par_x02 by parameter value paras[1]
    #                      ...
    if version in ['v5','v6']:
        # parameters are named in order and without gaps
        if len(paras) > 1 and len(paras) < 100:    # i want also <9 paras to be with leading zero here (otherwise replace ">1" with ">9")
            for ii in range(1,len(paras)+1):
                dict_paras['par_x{:02d}'.format(ii)] = paras[ii-1]
        elif len(paras) > 99 and len(paras) < 1000:
            for ii in range(1,len(paras)+1):
                dict_paras['par_x{:03d}'.format(ii)] = paras[ii-1]
        elif len(paras) <= 9:
            for ii in range(1,len(paras)+1):
                dict_paras['par_x{:01d}'.format(ii)] = paras[ii-1]
        else:
            raise ValueError("More than 999 parameters are not implemented yet!")
    elif version in ['v14','v15']:
        # parameters are named in order but with gaps (par_x03 and par_x10 missing)
        dict_paras['par_x01'] = paras[0]
        dict_paras['par_x02'] = paras[1]
        dict_paras['par_x04'] = paras[2]
        dict_paras['par_x05'] = paras[3]
        dict_paras['par_x06'] = paras[4]
        dict_paras['par_x07'] = paras[5]
        dict_paras['par_x08'] = paras[6]
        dict_paras['par_x09'] = paras[7]
        dict_paras['par_x11'] = paras[8]
        dict_paras['par_x12'] = paras[9]
        dict_paras['par_x13'] = paras[10]
        dict_paras['par_x14'] = paras[11]
        dict_paras['par_x15'] = paras[12]
        dict_paras['par_x16'] = paras[13]
        dict_paras['par_x17'] = paras[14]
        dict_paras['par_x18'] = paras[15]
        dict_paras['par_x19'] = paras[16]
        dict_paras['par_x20'] = paras[17]
        dict_paras['par_x21'] = paras[18]
        
    else:
        raise ValueError('Version not implemented')

    #print('dict_paras = ',dict_paras)

    # ---------------
    # create a run folder
    # ---------------
    tmp_folder = "/tmp/eee-analysis/"+str(run_id) # "/tmp/juletest" #  TODO a generic folder name in /tmp
    raven_exe_name   = os.path.abspath(dir_path+"/../../bc-hydro-raven/data_in/raven/Raven_v3.8.exe") #Raven.exe")
    raven_obs_folder = os.path.abspath(dir_path+"/../../bc-hydro-raven/data_in/ostrich_template_"+setupdate+"_"+version+"/"+basin+"/model")

    if os.path.exists(tmp_folder):
        shutil.rmtree(tmp_folder)

    # all Raven setup files
    tpl_files = glob.glob(str(Path(raven_obs_folder).parent)+'/*.rv*.tpl')

    RV_filebase = str( Path(Path(tpl_files[0]).stem).stem ).replace('_solution','')
    print('RV_filebase = ',RV_filebase)
    for tpl_file in tpl_files:
        
        with open(tpl_file, 'r') as ff: 
  
            # Reading the content of the file 
            # using the read() function and storing 
            # them in a new variable 
            data = ff.read() 
    
        # Searching and replacing the text 
        # using the replace() function
        for search_text in dict_paras.keys():
            replace_text = str(dict_paras[search_text])
            data = data.replace(search_text, replace_text)
            
            #print('replaced "{}" with "{}"'.format(search_text,replace_text))
  
        # Opening our text file in write only 
        # mode to write the replaced content
        out_file = Path(tmp_folder,Path(tpl_file).stem)
        makeDirectories([out_file])
        with open(out_file, 'w') as ff: 
  
            # Writing the replaced data in our 
            # text file 
            ff.write(data)

        print('Wrote: {}'.format(out_file))
            

    # link executable
    if not(os.path.exists(str(Path(tmp_folder,os.path.basename(raven_exe_name))))):
        print("from: ",os.path.realpath(raven_exe_name))
        print("to:   ",str(Path(tmp_folder,os.path.basename(raven_exe_name))))
        os.symlink(os.path.realpath(raven_exe_name), str(Path(tmp_folder,os.path.basename(raven_exe_name))))

    # link observations folder
    if not(os.path.exists(str(Path(tmp_folder,os.path.basename(raven_obs_folder))))):
        os.symlink(os.path.realpath(raven_obs_folder), str(Path(tmp_folder,os.path.basename(raven_obs_folder))))

    # create ouput folder
    out_folder = str(Path(tmp_folder,"output"))
    os.makedirs(out_folder)

    # ---------------
    # run the model with these input rv* files
    # ---------------
    cmd = [str(Path(tmp_folder,os.path.basename(raven_exe_name))),str(Path(tmp_folder,RV_filebase)),"-o",str(Path(tmp_folder,"output"))+'/']
    print("run cmd: ",' '.join(cmd))

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    print("")
    print("Raven standard output:")
    for line in process.stdout:
        print(">>> ",line.rstrip()) # rstrip removes trailing \n

    if not(os.path.exists(str(Path(tmp_folder,"output",RV_filebase+"_Diagnostics.csv")))):
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
    # extract model output: Diagnostics: NSE   (HYDROGRAPH_POST2K)
    # ---------------
    model['nse'] = 0.0
    ff = open(str(Path(tmp_folder,"output",RV_filebase+"_Diagnostics.csv")), "r")
    lines = ff.readlines()
    ff.close()

    nse = float(lines[-2].strip().split(',')[3])
    print("NSE:            ",nse)
    model['nse'] = nse
    print("")

    # ---------------
    # extract model output: Diagnostics: KGE   (HYDROGRAPH_POST2K)
    # ---------------
    model['kge'] = 0.0
    ff = open(str(Path(tmp_folder,"output",RV_filebase+"_Diagnostics.csv")), "r")
    lines = ff.readlines()
    ff.close()

    kge = float(lines[-2].strip().split(',')[7])
    print("KGE:            ",kge)
    model['kge'] = kge
    print("")

    # ---------------
    # extract model output: Hydrographs: simulated Q
    # ---------------
    model['Q']  = 0.0
    warmup = 2*365  # 1 # model timestep 1 day and want to skip 2 years  # first day 1991-01-01 00:00:00.00 (checked)
    model['Q']  = np.transpose(fread(str(Path(tmp_folder,"output",RV_filebase+"_Hydrographs.csv")),skip=warmup+1,cskip=4,nc=1))[0]

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
    skip = int(parasets[0].strip().split(':')[1])
else:
    skip = int(skip)
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
