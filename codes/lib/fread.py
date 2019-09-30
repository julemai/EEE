#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

__all__ = ['fread']

def fread(infile, nc=0, cname=None, skip=0, cskip=0, hskip=0, hstrip=True, separator=None,
          squeeze=False, reform=False, skip_blank=False, comment=None,
          fill=False, fill_value=0, strip=None, encoding='ascii', errors='ignore',
          header=False, full_header=False,
          transpose=False, strarr=False):
    """
        Read numbers into 2D-array from a file.

        This routines is exactly the same as sread but transforms
        everything to floats, handling NaN and Inf.


        Definition
        ----------
        def fread(infile, nc=0, cname=None, skip=0, cskip=0, hskip=0, hstrip=True, separator=None,
                  squeeze=False, reform=False, skip_blank=False, comment=None,
                  fill=False, fill_value=0, strip=None, encoding='ascii', errors='ignore',
                  header=False, full_header=False,
                  transpose=False, strarr=False):


        Input
        -----
        infile         source file name


        Optional Input Parameters
        -------------------------
        nc           number of columns to be read (default: all (nc<=0))
                     nc can be a vector of column indexes,
                     starting with 0; cskip will be ignored then.
        cname        columns can alternatively be chosen by the values in the first header line;
                     must be iterable with strings.
        skip         number of lines to skip at the beginning of file (default: 0)
        cskip        number of columns to skip at the beginning of each line (default: 0)
        hskip        number of lines in skip that do not belong to header (default: 0)
        hstrip       If true strip header cells to match with cname (default: True)
        separator    column separator
                     If not given, columns separator are (in order):
                     comma, semi-colon, whitespace
        comment      line gets excluded if first character of line is in comment sequence
                     sequence can be e.g. string, list or tuple
        fill_value   value to fill in if not enough columns in line
                     and fill=True (default: 0, and '' for header)
        strip        Strip strings with str.strip(strip).
                     If None then strip quotes " and ' (default).
                     If False then no strip (30% faster).
                     Otherwise strip character given by strip.
        encoding     Specifies the encoding which is to be used for the file.
                     Any encoding that encodes to and decodes from bytes is allowed.
                     (default: ascii)
        errors       Errors may be given to define the error handling during encoding of the file.
                     Possible values: strict, replace, ignore (default: ignore).


        Options
        -------
        squeeze      True:  2-dim array will be cleaned of degenerated
                            dimension, i.e. results in vector
                     False: array will be two-dimensional as read (default)
        reform       Same as squeeze.
        skip_blank   True:  continues reading after blank line
                     False: stops reading at first blank line (default)
        fill         True:  fills in fill_value if not enough columns in line
                     False: stops execution and returns None if not enough
                            columns in line (default)
        header       True:  header strings will be returned
                     False  numbers in file will be returned (default)
        full_header  True:  header is a string vector of the skipped rows
                     False: header will be split in columns,
                            exactly as the data, and will hold only the
                            selected columns (default)
        transpose    True:  column-major format output(0:ncolumns,0:nlines)
                     False: row-major format output(0:nlines,0:ncolumns) (default)
        strarr       True:  return header as numpy array of strings
                     False: return header as list


        Output
        ------
        Depending on options:
            2D-array of floats if header=False
            String array of file header if header=True
            String vector of file header if header=True and full_header=True


        Restrictions
        ------------
        If header=True then skip is counterintuitive because it is
          actally the number of header rows to be read. This is to
          be able to have the exact same call of the function, once
          with header=False and once with header=True.
        If fill=True, blank lines are not filled but are expected end of file.
        transpose=True has no effect on 1D output such as 1 header line


        Examples
        --------
        >>> # Create some data
        >>> filename = 'test.dat'
        >>> ff = open(filename,'w')
        >>> ff.writelines('head1 head2 head3 head4\\n')
        >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
        >>> ff.writelines('2.1 2.2 2.3 2.4\\n')
        >>> ff.close()

        >>> # Read sample file in different ways
        >>> # header
        >>> print(fread(filename, nc=2, skip=1, header=True))
        ['head1', 'head2']
        >>> print(fread(filename, nc=2, skip=1, header=True, full_header=True))
        ['head1 head2 head3 head4']
        >>> print(fread(filename, nc=1, skip=2, header=True))
        [['head1'], ['1.1']]
        >>> print(fread(filename, nc=1, skip=2, header=True, squeeze=True))
        ['head1', '1.1']
        >>> print(fread(filename, nc=1, skip=2, header=True, strarr=True))
        [['head1']
         ['1.1']]

        >>> # data
        >>> from autostring import astr
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=2), 1, pp=True))
        [['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, cskip=1), 1, pp=True))
        [['1.2' '1.3' '1.4']
         ['2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, nc=2, skip=1, cskip=1), 1, pp=True))
        [['1.2' '1.3']
         ['2.2' '2.3']]
        >>> print(astr(fread(filename, nc=[1,3], skip=1), 1, pp=True))
        [['1.2' '1.4']
         ['2.2' '2.4']]
        >>> print(astr(fread(filename, nc=1, skip=1), 1, pp=True))
        [['1.1']
         ['2.1']]
        >>> print(astr(fread(filename, nc=1, skip=1, reform=True), 1, pp=True))
        ['1.1' '2.1']

        >>> # skip blank lines
        >>> ff = open(filename, 'a')
        >>> ff.writelines('\\n')
        >>> ff.writelines('3.1 3.2 3.3 3.4\\n')
        >>> ff.close()
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True, comment='#!'), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']
         ['3.1' '3.2' '3.3' '3.4']]

        >>> # skip comment lines
        >>> ff = open(filename, 'a')
        >>> ff.writelines('# First comment\\n')
        >>> ff.writelines('! Second 2 comment\\n')
        >>> ff.writelines('4.1 4.2 4.3 4.4\\n')
        >>> ff.close()
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, nc=[2], skip_blank=True, comment='#'), 1, pp=True))
        [['1.3']
         ['2.3']
         ['3.3']
         ['2.0']
         ['4.3']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True, comment='#!'), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']
         ['3.1' '3.2' '3.3' '3.4']
         ['4.1' '4.2' '4.3' '4.4']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True, comment=('#','!')), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']
         ['3.1' '3.2' '3.3' '3.4']
         ['4.1' '4.2' '4.3' '4.4']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True, comment=['#','!']), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']
         ['3.1' '3.2' '3.3' '3.4']
         ['4.1' '4.2' '4.3' '4.4']]

        >>> # fill missing columns
        >>> ff = open(filename, 'a')
        >>> ff.writelines('5.1 5.2\\n')
        >>> ff.close()
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, skip_blank=True, comment='#!', fill=True, fill_value=-1), 1, pp=True))
        [[' 1.1' ' 1.2' ' 1.3' ' 1.4']
         [' 2.1' ' 2.2' ' 2.3' ' 2.4']
         [' 3.1' ' 3.2' ' 3.3' ' 3.4']
         [' 4.1' ' 4.2' ' 4.3' ' 4.4']
         [' 5.1' ' 5.2' '-1.0' '-1.0']]

        >>> # transpose
        >>> print(astr(fread(filename, skip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(astr(fread(filename, skip=1, transpose=True), 1, pp=True))
        [['1.1' '2.1']
         ['1.2' '2.2']
         ['1.3' '2.3']
         ['1.4' '2.4']]

        >>> # Create some more data with Nan and Inf
        >>> filename1 = 'test1.dat'
        >>> ff = open(filename1, 'w')
        >>> ff.writelines('head1 head2 head3 head4\\n')
        >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
        >>> ff.writelines('2.1 nan Inf "NaN"\\n')
        >>> ff.close()

        >>> # Treat Nan and Inf with automatic strip of " and '
        >>> print(astr(fread(filename1, skip=1, transpose=True), 1, pp=True))
        [['1.1' '2.1']
         ['1.2' 'nan']
         ['1.3' 'inf']
         ['1.4' 'nan']]

        >>> # Create some more data with escaped numbers
        >>> filename2 = 'test2.dat'
        >>> ff = open(filename2, 'w')
        >>> ff.writelines('head1 head2 head3 head4\\n')
        >>> ff.writelines('"1.1" "1.2" "1.3" "1.4"\\n')
        >>> ff.writelines('2.1 nan Inf "NaN"\\n')
        >>> ff.close()

        >>> # Strip
        >>> print(astr(fread(filename2,  skip=1,  transpose=True,  strip='"'), 1, pp=True))
        [['1.1' '2.1']
         ['1.2' 'nan']
         ['1.3' 'inf']
         ['1.4' 'nan']]

        >>> # Create some more data with an extra (shorter) header line
        >>> filename3 = 'test3.dat'
        >>> ff = open(filename3, 'w')
        >>> ff.writelines('Extra header\\n')
        >>> ff.writelines('head1 head2 head3 head4\\n')
        >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
        >>> ff.writelines('2.1 2.2 2.3 2.4\\n')
        >>> ff.close()

        >>> print(astr(fread(filename3, skip=2, hskip=1), 1, pp=True))
        [['1.1' '1.2' '1.3' '1.4']
         ['2.1' '2.2' '2.3' '2.4']]
        >>> print(fread(filename3, nc=2, skip=2, hskip=1, header=True))
        ['head1', 'head2']

        >>> # cname
        >>> print(astr(fread(filename, cname='head2', skip=1, skip_blank=True, comment='#!', squeeze=True), 1, pp=True))
        ['1.2' '2.2' '3.2' '4.2' '5.2']
        >>> print(astr(fread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!'), 1, pp=True))
        [['1.1' '1.2']
         ['2.1' '2.2']
         ['3.1' '3.2']
         ['4.1' '4.2']
         ['5.1' '5.2']]
        >>> print(fread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True))
        ['head1', 'head2']
        >>> print(fread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True, full_header=True))
        ['head1 head2 head3 head4']
        >>> print(astr(fread(filename, cname=['  head1','head2'], skip=1, skip_blank=True, comment='#!', hstrip=False), 1, pp=True))
        [['1.2']
         ['2.2']
         ['3.2']
         ['4.2']
         ['5.2']]

        >>> # Clean up doctest
        >>> import os
        >>> os.remove(filename)
        >>> os.remove(filename1)
        >>> os.remove(filename2)
        >>> os.remove(filename3)


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2009-2017 Matthias Cuntz - mc (at) macu (dot) de

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.


        History
        -------
        Written,  MC, Jul 2009
        Modified, MC, Feb 2012 - transpose
                  MC, Feb 2013 - ported to Python 3
                  MC, Nov 2014 - bug when nc is list and contains 0
                  MC, Nov 2014 - hskip
                  MC, Feb 2015 - speed: everything list until very end
                  MC, Nov 2017 - use range instead of np.arange for producing indexes
                  MC, Nov 2017 - cname, sname, file->infile, hstrip
                  MC, Jun 2019 - open(errors='ignore') to ignore unicode characters, for example, on read
                  MC, Jul 2019 - errors='ignore' compatible with Python2 and Python3
                                 -> returns header in unicode in Python2
                  MC, Aug 2019 - use codecs module and allow user encoding and error handling
    """
    #
    # Open file
    import codecs
    f = codecs.open(infile, 'r', encoding=encoding, errors=errors)
    #
    # Read header and Skip lines
    if hskip > 0:
        ihskip = 0
        while ihskip < hskip:
            tmp = f.readline().rstrip()
            ihskip += 1
    if skip > 0:
        head = ['']*(skip-hskip)
        iskip = 0
        while iskip < (skip-hskip):
            head[iskip] = str(f.readline().rstrip())
            iskip += 1
    #
    # read first line to determine nc and separator (if not set)
    split = -1
    while True:
        s = f.readline().rstrip()
        if len(s) == 0:
            if skip_blank:
                continue
            else:
                break
        if comment is not None:
            if (s[0] in comment): continue
        break
    if separator is None:
        sep = ','
        res = s.split(sep)
        nres = len(res)
        if nres == 1:
            sep = ';'
            res = s.split(sep)
            nres = len(res)
            if nres == 1:
                sep = None
                res = s.split(sep)
                nres = len(res)
    else:
        sep = separator
        res = s.split(sep)
        nres = len(res)
    #
    # Determine indices
    if nc != 0 and cname is not None:
        f.close()
        raise ValueError('nc and cname are mutually exclusive.')
    if cname is not None:
        # from first header line
        if (skip-hskip) <= 0:
            f.close()
            raise IOError('No header line left for choosing columns by name.')
        if not isinstance(cname, (list, tuple, np.ndarray)): cname = [cname]
        if hstrip: cname = [ h.strip() for h in cname ]
        hres = head[0].split(sep)
        if hstrip: hres = [ h.strip() for h in hres ]
        iinc = []
        for k in range(len(hres)):
            if hres[k] in cname: iinc.append(k)
        nnc = len(iinc)
    else:
        # from nc keyword
        if isinstance(nc, (list, tuple, np.ndarray)):
            nnc  = len(nc)
            iinc = tuple(nc)
        else:
            if nc <= 0:
                iinc = range(cskip,nres)
                nnc  = nres-cskip
            else:
                iinc = range(cskip,cskip+nc)
                nnc = nc
    miinc = max(iinc)
    #
    # Header
    if header:
        # Split header
        var = None
        if (skip-hskip) > 0:
            if full_header:
                var = head
            else:
                var = list()
                k = 0
                while k < (skip-hskip):
                    hres = head[k].split(sep)
                    nhres = len(hres)
                    if (miinc >= nhres) and (not fill):
                        f.close()
                        raise ValueError('Line has not enough columns to index - 01: '+head[k])
                    null = line2var(hres, var, iinc, strip)
                    k += 1
                if (skip-hskip) == 1: var = var[0]
            f.close()
        else:
            var = None
            f.close()
            return var
        if strarr:
            var = np.array(var, dtype=np.str)
            if transpose: var = var.T
            if squeeze or reform: var = var.squeeze()
            if fill: var = np.where(var=='', fill_value, var)
        else:
            if fill:
                var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
            if squeeze or reform:
                maxi = max([ len(i) for i in var])
                if maxi==1: var = [ i[0] for i in var ]
            if transpose and isinstance(var[0], list):
                var = [ list(i) for i in zip(*var) ] # transpose
        return var
    #
    # Values - first line
    if (miinc >= nres) and (not fill):
        f.close()
        raise ValueError('Line has not enough columns to index - 02: '+s)
    var = list()
    null = line2var(res, var, iinc, strip)
    #
    # Values - rest of file
    for line in f:
        s = str(line.rstrip())
        if len(s) == 0:
            if skip_blank:
                continue
            else:
                break
        if comment is not None:
            if (s[0] in comment): continue
        res = s.split(sep)
        nres = len(res)
        if (miinc >= nres) and (not fill):
            f.close()
            raise ValueError('Line has not enough columns to index - 03: '+s)
        null = line2var(res, var, iinc, strip)

    f.close()
    # list -> array
    if fill:
        var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
    var = np.array(var, dtype=np.float)
    if squeeze or reform: var = var.squeeze()
    if transpose: var = var.T

    return var


# Helper for append var with current line already splitted into list
def line2var(res, var, iinc, strip):
    nres = len(res)
    if strip is None:
        tmp = [res[i].strip('"').strip("'") for i in iinc if i < nres]
    elif not strip:
        tmp = [res[i] for i in iinc if i < nres]
    else:
        tmp = [res[i].strip(strip) for i in iinc if i < nres]
    rest = len([ i for i in iinc if i >= nres ])
    if rest > 0:
        tmp.extend(['']*rest)
    var.append(tmp)
    return


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from autostring import astr
    # filename = 'test.dat'
    # ff = open(filename,'w')
    # ff.writelines('head1 head2 head3 head4\n')
    # ff.writelines('1.1 1.2 1.3 1.4\n')
    # ff.writelines('2.1 2.2 2.3 2.4\n')
    # ff.writelines('\n')
    # ff.writelines('3.1 3.2 3.3 3.4\n')
    # ff.writelines('# First comment\n')
    # ff.writelines('! Second 2 comment\n')
    # ff.writelines('4.1 4.2 4.3 4.4\n')
    # ff.close()
    # print(astr(fread(filename, skip=1), 1, pp=True))
    # print(astr(fread(filename, skip=1, nc=[2], skip_blank=True, comment='#'), 1, pp=True))
    # print(astr(fread(filename, skip=1, skip_blank=True, comment='#!'), 1, pp=True))
    # print(astr(fread(filename, cname='head2', skip=1, skip_blank=True, comment='#!'), 1, pp=True))
    # print(astr(fread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!'), 1, pp=True))
    # print(fread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True))
    # print(fread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True, full_header=True))
