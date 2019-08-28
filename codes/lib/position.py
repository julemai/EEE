#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

__all__ = ['position']

def position(row=1, col=1, num=1,
             left=0.125, right=0.9, bottom=0.1, top=0.9,
             hspace=0.1, vspace=None, wspace=None,
             width=None, height=None,
             sortcol=False, golden=False, inversegolden=False,
             figsize=(1.,1.)):
    """
        Gives positions of subplots.
        To be used with add_axes instead of subplot.

        All dimensions are fractions of the figure width or height.
        Figure and subplot spaces are the same as for figure.subplotparams
        except for hspace and vspace, which are halved.

        If the figsize keyword is given, a rectangular section of the figure
        will be used.


        Definition
        ----------
        def position(row=1, col=1, num=1,
                     left=0.125, right=0.9, bottom=0.1, top=0.9,
                     hspace=0.1, vspace=None, wspace=None,
                     width=None, height=None,
                     sortcol=False, golden=False, inversegolden=False,
                     figsize=(1.,1.)):


        Optional Input
        --------------
        row            number of subplot rows (default 1)
        col            number of subplot columns (default 1)
        num            subplot number (default 1)
        left           left border of plot (default 0.125)
        right          right border of plot (default 0.9)
        bottom         bottom border of plot (default 0.1)
        top            top border of plot (default 0.9)
        hspace         space between columns (default 0.1)
        vspace         space between rows (default 0.1)
        wspace         historical, same as vspace; will be overwritten by vspace
        width          prescribe width of plots (default None)
        height         prescribe height of plots (default None)
        sortcol        fill columns then rows (default False)
        golden         golden ratio of width/height = (1+sqrt(5))/2
                       (default False)
        inversegolden  golden ratio of height/width
                       (overwritten by golden) (default False)
        figsize        (width, height) of figure as given by e.g.
                       matplotlib.rcParams['figure.figsize'].
                       Scales everything to rectangular section
                       (default (1,1))


        Output
        ------
        position array with [left, bottom, width, height)
        to be used with fig.add_axes.


        Examples
        --------
        # Use, for example, as follows
        # fig1 = figure(1)
        # sub1 = fig1.add_axes(position(2,2,1))
        # sub2 = fig1.add_axes(position(2,2,2))

        # if you want to have a true rectangle
        # figsize = matplotlib.rcParams['figure.figsize']
        # sub = fig1.add_axes(position(1,1,1,figsize=figsize,left=0.1))

        # if you want to have a true golden ratio
        # sub = fig1.add_axes(position(1,1,1,figsize=figsize,golden=True))

        # Doctest examples
        >>> from autostring import astr
        >>> print(astr(position(2,2,1),3,pp=True))
        ['0.125' '0.550' '0.338' '0.350']
        >>> print(astr(position(2,2,1,sortcol=True),3,pp=True))
        ['0.125' '0.550' '0.338' '0.350']
        >>> print(astr(position(2,2,1,golden=True),3,pp=True))
        ['0.125' '0.409' '0.338' '0.209']
        >>> print(astr(position(2,2,1,inversegolden=True),3,pp=True))
        ['0.125' '0.550' '0.216' '0.350']
        >>> print(astr(position(2,2,1,golden=True,sortcol=True),3,pp=True))
        ['0.125' '0.409' '0.338' '0.209']
        >>> print(astr(position(2,2,1,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,vspace=0.),3,pp=True))
        ['0.000' '0.500' '0.500' '0.500']
        >>> print(astr(position(2,2,2,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,vspace=0.),3,pp=True))
        ['0.500' '0.500' '0.500' '0.500']
        >>> print(astr(position(2,2,3,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,vspace=0.),3,pp=True))
        ['0.000' '0.000' '0.500' '0.500']
        >>> print(astr(position(2,2,4,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,vspace=0.),3,pp=True))
        ['0.500' '0.000' '0.500' '0.500']
        >>> print(astr(position(2,2,1,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,vspace=0.,golden=True),3,pp=True))
        ['0.000' '0.309' '0.500' '0.309']
        >>> print(astr(position(2,2,2,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,vspace=0.,golden=True),3,pp=True))
        ['0.500' '0.309' '0.500' '0.309']
        >>> print(astr(position(2,2,3,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,vspace=0.,golden=True),3,pp=True))
        ['0.000' '0.000' '0.500' '0.309']
        >>> print(astr(position(2,2,4,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,vspace=0.,golden=True),3,pp=True))
        ['0.500' '0.000' '0.500' '0.309']
        >>> figsize=[8,11]
        >>> print(astr(position(2,2,1,golden=True,sortcol=True,figsize=figsize),3,pp=True))
        ['0.125' '0.324' '0.338' '0.152']
        >>> print(astr(position(2,2,1,figsize=figsize,left=0.1),3,pp=True))
        ['0.100' '0.427' '0.350' '0.255']
        >>> print(astr(position(2,2,1,figsize=figsize,left=0.1,golden=True),3,pp=True))
        ['0.100' '0.330' '0.350' '0.157']


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2009-2016 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Aug 2009
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Jul 2013 - vspace, wspace obsolete
                  MC, Apr 2014 - assert
                  ST, Feb 2016 - added height and width
    """
    #
    # Check
    nplots = row*col
    assert num <= nplots, 'num > number of plots: '+str(num)+' > '+str(nplots)
    assert right-left > 0., 'right > left: '+str(right)+' > '+str(left)
    assert top-bottom > 0., 'top < bottom: '+str(top)+' < '+str(bottom)
    if vspace != None:
        ivspace = vspace
    elif wspace != None:
        ivspace = wspace
    else:
        ivspace = 0.1
    #
    # Scaling to figsize
    scalex = figsize[1]/float(max(figsize))
    scaley = figsize[0]/float(max(figsize))
    #
    # width, height
    if width is None:
        dx = (right-left-(col-1)*hspace)/col
    else:
        dx = width
    if height is None:
        dy = (top-bottom-(row-1)*ivspace)/row
    else:
        dy = height
    #
    # golden ratio
    ratio = (1.+np.sqrt(5.))/2.
    if golden:
        width  = dx
        height = dx / ratio
        checkheight = (top-bottom-row*height) - (row-1)*ivspace
        if checkheight < 0.:
            height = dy
            width  = dy * ratio
            checkwidth = (right-left-col*width) - (col-1)*hspace
            if checkwidth < 0.:
                raise ValueError('golden ratio does not work. Have to recode.')
    else:
        if inversegolden:
            height = dy
            width  = dy / ratio
            checkwidth = (right-left-col*width) - (col-1)*hspace
            if checkwidth < 0.:
                width  = dx
                height = dx * ratio
                checkheight = (top-bottom-row*height) - (row-1)*ivspace
                if checkheight < 0.:
                    raise ValueError('inverse golden ratio does not work. Have to recode.')
        else:
            width  = dx
            height = dy
    #
    # order row/colmn, column/row
    if sortcol:
        irow = (num-1) % row
        icol = (num-1) // row
    else:
        irow = (num-1) // col
        icol = (num-1) % col
    #
    # position
    pos    = np.empty(4)
    pos[0] = left   + icol*(width+hspace)           *scalex
    pos[1] = bottom + (row-1-irow)*(height+ivspace) *scaley
    pos[2] = width  *scalex
    pos[3] = height *scaley
    #
    return pos

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

