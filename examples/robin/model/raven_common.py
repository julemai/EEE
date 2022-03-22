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
Helper functions for Sobol' sensitivity analysis of RAVEN

History
-------
Written,  JM, Jun 2019
"""

from pathlib2 import Path

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
