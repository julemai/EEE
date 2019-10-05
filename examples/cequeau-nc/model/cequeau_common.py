#!/usr/bin/env python

# Copyright 2019 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of Juliane Mai's personal code library.
#
# Juliane Mai's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Juliane Mai's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public Licensefstop
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#

from __future__ import print_function

"""
Helper functions for Sobol' sensitivity analysis of CEQUEAU

History
-------
Written,  JM, Jun 2019
"""

from pathlib2   import Path
import netcdf4  as nc4     # in lib/
import numpy as np

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

def get_discharge(start_day, ncfile, duration=1, group=None, var=None, ibasin=None, ilag=None):

    # ibasin ... index of basin (starts with 0)

    ncdata = nc4.NcDataset(ncfile, "r")
    ncdata_time = ncdata.variables["t"][:]
    ncdata_time_unit = ncdata.variables['t'].units
    ncdata_time_cal  = ncdata.variables['t'].calendar

    # -------------------
    # below takes 0.1925 sec  (total runtime = 0.2065 sec)
    # -------------------
    # ncdata_time = nc4.netcdf4.num2date(ncdata_time,units=ncdata_time_unit,calendar=ncdata_time_cal)          # <<<< this is awefully slow
    
    # -------------------
    # below takes 0.000126 sec (total runtime = 0.00997 sec  --> reduced by 95.2%)
    # -------------------
    start_day = nc4.netcdf4.date2num(start_day,units=ncdata_time_unit,calendar=ncdata_time_cal)

    if ( len(np.where( ncdata_time == start_day )[0]) == 1):

        idx = np.where( ncdata_time == start_day )[0][0]

        if ( not(group is None) and not(var is None) and not(ibasin is None) and not(ilag is None) ):
            ncdata_var = ncdata.groups[group].variables[var][idx:idx+duration,ibasin,ilag]
        elif ( not(var is None) and not(ibasin is None) ):
            ncdata_var = ncdata.variables[var][ibasin,idx:idx+duration]
        else:
            raise ValueError('common: get_discharge: this usage of get_discharge() is not implemented yet')
            
    else:

        print('ncdata_time = ',ncdata_time)
        print('start_day=',start_day,' found ',len(np.where( ncdata_time == start_day )[0]),' times')
        raise ValueError('common: get_discharge: start_day either found multiple times or not at all')

    if len(ncdata_var) == 1:
        return ncdata_var[0]
    else:
        return ncdata_var
