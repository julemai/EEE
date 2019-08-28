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


# An example calling sequence to derive the cutoff threshold is given below. The Elementary Effects need to
# be determined before and stored in a file (option -e). The model parameters need to be specified (option -m).
# Results can be plotted as a PDF (option -p) or suppressed (option -n). LaTeX fonts is used in plot if
# option '-t' is set. In case the cutoff is already known, it can be specified (option -c). In that case only
# the parameters that are below the threshold are determined, plots are created, and a new parameter info file
# is written (last column will change constantly).
#
# python 4_derive_threshold.py \
#                -e example_ishigami-homma/eee_results.dat
#                -m example_ishigami-homma/parameters.dat
#                -p example_ishigami-homma/eee_results.pdf
#                -c -1
#                -t
#

"""
This script is to derive the cutoff threshold is given below. The Elementary Effects need to 
be determined before and stored in a file (option -e). The model parameters need to be specified (option -m). 
Results can be plotted as a PDF (option -p) or suppressed (option -n). LaTeX fonts is used in plot if 
option '-t' is set. In case the cutoff is already known, it can be specified (option -c). In that case only 
the parameters that are below the threshold are determined, plots are created, and a new parameter info file 
is written (last column will change constantly).

In case multiple model outputs are considered, a multi-objective approach is applied.

History
-------
Written,  JM, Mar 2019
"""
        
# -------------------------------------------------------------------------
# Command line arguments
# -------------------------------------------------------------------------

cutoff             = '-1'
eefile             = 'eee_results.dat'
maskfile           = 'parameters.dat'
pdffile            = 'eee_results.pdf'
usetex             = False
noplot             = False
multi_obj_approach = 'triangle'    # rectangle :: one of the cutoff values must be exceeded to be tagged as informative
                                   # triangle  :: combination of all cutoffs lies below the hyperplane through the single cuttoffs

import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="Plotting sorted Elementary Effects and mark cut-off for discarding parameters for next iteration of MUCM.")

# run plot_ee.py 
parser.add_option('-c', '--cutoff', action='store', dest='cutoff', type='string',
                  default=cutoff, metavar='Cutoff value',
                  help='Cut-off value, i.e. all parameters with Elementary Effect above this threshold will be discarded for next iteration. Multiple cutoffs are separeted by colons (default: cutoff=-1).')
parser.add_option('-e', '--eefile', action='store', dest='eefile', type='string',
                  default=eefile, metavar='File with elementary effects',
                  help='File where Elementary Effects are stored. (default: eefile=ee_results.dat).')
parser.add_option('-m', '--maskfile', action='store', dest='maskfile', type='string',
                  default=maskfile, metavar='File with 0/1 for parameters',
                  help='Name of file where masked parameters are specified (default: maskfile=mask_para.dat).')
parser.add_option('-p', '--pdffile', action='store', dest='pdffile', type='string',
                  default=pdffile, metavar='PDF output file',
                  help='Name of pdf output file (default: open X-window).')
parser.add_option('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                  help="Use LaTeX to render text in pdf.")
parser.add_option('-n', '--noplot', action='store_true', default=noplot, dest="noplot",
                  help="No plot will be produced (for running on GridEngine).")
(opts, args) = parser.parse_args()

cutoff    = opts.cutoff
eefile    = opts.eefile    # file containing Elementary Effects
maskfile  = opts.maskfile  # mask_para.dat
pdffile   = opts.pdffile   # pdf file for plots
usetex    = opts.usetex    # if tex modus should be used
noplot    = opts.noplot    # if plot will be produced

del parser, opts, args

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/lib')

import numpy as np
import scipy.optimize as opt

from general_functions import curvature, logistic_offset_p, dlogistic, d2logistic  # in lib/
from fit_functions     import cost_square                   # in lib/
from fsread            import fsread                        # in lib/
from autostring        import astr                          # in lib/
from position          import position                      # in lib/
from str2tex           import str2tex                       # in lib/

# Func to minimize curvature
def mcurvature(*args, **kwargs):
    return -curvature(*args, **kwargs)

# -------------------------------------------------------------------------
# Read parameter info file
# -------------------------
# parameter info file has following header:
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



# -------------------------------------------------------------------------
# Read Elementary Effects file
# -------------------------------------------------------------------------
ee        = fsread(eefile,cskip=2,comment='#')
nobj      = np.int(np.shape(ee)[1]/2)
ee        = ee[:,0:nobj]
ee_masked = ee[mask_para]
obj_names = []
f = open(eefile, 'r')
for iobj in range(nobj):
    line = f.readline()
    obj_names.append(line.split(':')[1].strip())
f.close()
print('Number of objectives :: ',nobj)
print('Names of objectives  :: ',obj_names)
print('')

# -------------------------------------------------------------------------
# Customize plots
# -------------------------------------------------------------------------
if (not(noplot)):
    # Plot - paper_plots, but also all if not otherwise defined
    nrow        = 1           # # of rows per figure
    ncol        = 1           # # of columns per figure
    hspace      = 0.12        # x-space between plots
    vspace      = 0.12        # y-space between plots
    textsize    = 12          # Standard text size
    dt          = 4           # # of hours between tick marks on plots
    dxabc       = 0.85        # % shift from left y-axis of a,b,c,... labels
    dyabc       = 0.85        # % shift from lower x-axis of a,b,c,... labels
    dyabcdown   = 0.05        # y-shift if abc in lower right corner

    lwidth      = 1.5         # linewidth
    elwidth     = 1.0         # errorbar line width
    alwidth     = 1.0         # axis line width
    msize       = 5.0         # marker size
    mwidth      = 1.5         # marker edge width
    # color: 'b'|'g'|'r'|'c'|'m'|'y'|'k'|'w'
    #        'blue'|'green'|'red'|'cyan'|'magenta'|'yellow'|'black'|'white'
    #        hex string '#eeefff' | RGB tuple (1,0.5,1) | html names 'burlywod', 'chartreuse', ...
    #        grayscale intensity, e.g. '0.7', 'k'='0.0'
    mcol1       = 'black'       # primary marker colour
    mcol2       = 'gray'       # color of second markers
    mcol3       = '0.0'       # color of third markers
    lcol1       = 'black'      # primary line colour
    lcol2       = 'gray'
    lcol3       = '0.0'       # color of third lines

    llxbbox     = 0.5          # y-anchor legend bounding box
    llybbox     = 0.95        # y-anchor legend bounding box
    llrspace    = 0.          # spacing between rows in legend
    llcspace    = 1.0         # spacing between columns in legend
    llhtextpad  = 0.4         # the pad between the legend handle and text
    llhlength   = 1.5         # the length of the legend handles
    frameon     = False       # if True, draw a frame around the legend. If None, use rc
    llxbbox2    = 0.60       # Tight bounding of symbol and text (w/o lines)
    llhtextpad2 = 0.         #                   "
    llhlength2  = 1.0        #                   "

    if (pdffile == ''):
        outtype = 'x'
    else:
        outtype = 'pdf'

    import matplotlib as mpl
    if (outtype == 'pdf'):
        mpl.use('PDF') # set directly after import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        # Customize: http://matplotlib.sourceforge.net/users/customizing.html
        mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
        mpl.rc('figure', figsize=(8.27,11.69/5)) # a fifth of a4 portrait
        if usetex:
            mpl.rc('text', usetex=True)
            # mpl.rc('text.latex', unicode=True)
            mpl.rcParams['text.latex.preamble']=r'\usepackage{wasysym}'
        else:
            #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
            #mpl.rc('font',**{'family':'serif','serif':['times']})
            mpl.rcParams['font.family'] = 'serif'
            mpl.rcParams['font.serif']  = 'Times'
        mpl.rc('font', size=textsize)
    else:
        import matplotlib.pyplot as plt
        #mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
        mpl.rc('figure', figsize=(4./5.*8.27,4./5.*4.69/5)) # a fifth of a4 portrait
        mpl.rc('font', size=textsize)
    mpl.rc('lines', linewidth=lwidth, color='black')
    mpl.rc('axes', linewidth=alwidth, labelcolor='black')
    mpl.rc('path', simplify=False) # do not remove

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------
if (not(noplot)):
    if (outtype == 'pdf'):
        print('Plot PDF ', pdffile)
        pdf_pages = PdfPages(pdffile)
    else:
        print('Plot X')
    figsize = mpl.rcParams['figure.figsize']

    ifig = 0

# -------------------------------------------------------------------------
# Fig 1 - sorted Elementary Effects
#
if (multi_obj_approach == 'rectangle'):
    keepit = []
cutoff_obj = np.ones(nobj) * -9999.
for iobj in range(nobj):

    sort_idx  = np.argsort(ee_masked[:,iobj])
    print('')
    print('OBJECTIVE #',iobj+1)
    print('sorted ee:              ',astr(ee_masked[sort_idx,iobj],prec=4))
    print('correspond to para:     ',idx_para[sort_idx]+1)
    
    if (not(noplot)):
        ifig += 1
        iplot = 0
        # print('      Plot - Fig ', ifig)
        fig = plt.figure(ifig)

    if (not(noplot)):
        iplot += 1
        ylab  = r'$\mathrm{EE}$'
        sub   = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

    if cutoff == '-1':
        # normalise so that we can search for y' = 1
        xx = np.arange(dims)/np.float(dims-1)
        yy = ee_masked[sort_idx,iobj]/np.amax(ee_masked[:,iobj])   # yy = np.sort(mustar) / mumax

        # --------------------------
        # calculate cutoff
        # --------------------------
        pini=np.array([1.0, 10.0, 0.5, 0.0]) # [y-max, steepness, inflection point, offset]
        popti, f, d = opt.fmin_l_bfgs_b(cost_square, pini,
                                        args=(logistic_offset_p,xx,yy),
                                        approx_grad=1,
                                        bounds=[(None,None),(None,None),(None,None),(None,None)],#,(0.,0.2)],
                                        iprint=0, disp=0)
        print('Params of logistic function: ', astr(popti,prec=4))

        # fitted line
        yy2 = logistic_offset_p(xx,popti)
        if (not(noplot)):
            nxx = 100
            xx2 = np.arange(nxx) / np.float(nxx-1)
            yy2 = logistic_offset_p(xx2, popti)
            line2 = plt.plot(xx2, yy2)
            plt.setp(line2, linestyle='-', linewidth=lwidth, color=lcol1, marker='None', label=str2tex('$L$',usetex=usetex))
        
        # steepest curvature
        x_scaled = opt.brent(mcurvature, # minimizes
                             args=(dlogistic, d2logistic, popti[0], popti[1], popti[2]),
                             brack=(xx[0],xx[-1]))
        curvatures = logistic_offset_p(x_scaled, popti)
        gx_scaled = curvatures         # g(n_thresh)
        mumax = np.amax(ee_masked[:,iobj])
        if (logistic_offset_p(x_scaled, popti) > 0.2) or (x_scaled < xx[0]):
            x_scaled  = xx[0]
            gx_scaled = np.min(ee_masked[:,iobj]) / mumax
        cutoff1 = gx_scaled                                                # in 0-1 range # g(n_thresh)
        cutoff_obj[iobj]  = cutoff1*mumax                                  # in EE range  # g(n_thresh)*mu_thresh
        x33 = np.arange(0.0125,1.03,0.0345)
        y33 = np.array([curvatures for ii in range(np.shape(x33)[0])])
        line6 = plt.plot(x33, y33)
        plt.setp(line6, linestyle='None', linewidth=lwidth, color=lcol1, marker='x', markeredgecolor=mcol1, markerfacecolor='None',
                 markersize=msize/2, markeredgewidth=mwidth/2, label=str2tex('$L(x_k)$',usetex=usetex))

        print('Cutoff(s): ', astr(cutoff_obj[iobj],prec=4))

        # This selects parameters where ORIGINAL ee is above threshold
        if (multi_obj_approach == 'rectangle'):
            new = list( idx_para[sort_idx[ee_masked[sort_idx,iobj] < cutoff_obj[iobj]]] )
            if (iobj == 0):
                keepit = new
            else:
                # intersection of parameters
                keepit = list(set(keepit) & set(new))
            print('')
            print('Parameters for next iteration: ', astr(np.sort(keepit)+1),'   --> ',np.size(keepit),' parameters')

        if (not(noplot)):
            # threshold
            xmin, xmax = sub.get_xlim()
            line3 = plt.plot([xmin,xmax], [cutoff1,cutoff1])
            plt.setp(line3, linestyle='--', linewidth=lwidth, color=lcol2, marker='None', label=str2tex('$\eta^*_{thres}$',usetex=usetex))

        if (not(noplot)):
            sub.text(1.02, 0.5, str2tex(obj_names[iobj],usetex=usetex),
                         rotation=90, fontsize='large',
                         horizontalalignment='left', verticalalignment='center',
                         transform=sub.transAxes)

        if (not(noplot)):
            xnames     = str2tex([ '$'+iparname+'$' for iparname in para_name_mask[sort_idx] ],usetex=usetex)
            xlabel     = str2tex('Parameter Name',usetex=usetex)
            ylabel     = str2tex('$\eta^*$',usetex=usetex)
            plt.setp(sub, xticks=xx, xticklabels=xnames, xlabel=xlabel, ylabel=ylabel)

            # Elementary Effects (informative)
            iidx = np.where(yy >= cutoff1)
            mark1 = plt.plot(xx[iidx], yy[iidx])
            plt.setp(mark1, linestyle='None', marker='o', markeredgecolor=mcol1, markerfacecolor='None',
                 markersize=msize, markeredgewidth=mwidth, label=str2tex('retained',usetex=usetex))
            #linestyle='-', linewidth=lwidth, color=lcol1, marker='None')

            # Elementary Effects (non-informative)
            iidx = np.where(yy < cutoff1)
            mark2 = plt.plot(xx[iidx], yy[iidx])
            plt.setp(mark2, linestyle='None', marker='o', markeredgecolor=mcol2, markerfacecolor='None',
                 markersize=msize, markeredgewidth=mwidth, label=str2tex('non-influential',usetex=usetex))
            #linestyle='-', linewidth=lwidth, color=lcol1, marker='None')

            xlim = [0-1./(dims-1),1+1./(dims-1)]
            plt.setp(sub,  xlim=xlim)
            
    else:
        # Split the given string
        cutoff_obj[iobj] = np.float(cutoff.split(':')[iobj])    
        
        xx = np.arange(dims)
        yy = ee_masked[sort_idx,iobj]

        if (not(noplot)):
            xnames     = str2tex([ '$'+iparname+'$' for iparname in para_name_mask[sort_idx] ],usetex=usetex)
            xlabel     = str2tex('Parameter Name',usetex=usetex)
            ylabel     = str2tex('$\mu^*$',usetex=usetex)
            plt.setp(sub, xticks=xx, xticklabels=xnames, xlabel=xlabel, ylabel=ylabel)

            # Elementary Effects (informative)
            iidx = np.where(yy >= cutoff_obj[iobj])
            mark1 = plt.plot(xx[iidx], yy[iidx])
            plt.setp(mark1, linestyle='None', marker='o', markeredgecolor=mcol1, markerfacecolor='None',
                 markersize=msize, markeredgewidth=mwidth, label=str2tex('retained',usetex=usetex))
            #linestyle='-', linewidth=lwidth, color=lcol1, marker='None')

            # Elementary Effects (non-informative)
            iidx = np.where(yy < cutoff_obj[iobj])
            mark2 = plt.plot(xx[iidx], yy[iidx])
            plt.setp(mark2, linestyle='None', marker='o', markeredgecolor=mcol2, markerfacecolor='None',
                 markersize=msize, markeredgewidth=mwidth, label=str2tex('non-influential',usetex=usetex))
            #linestyle='-', linewidth=lwidth, color=lcol1, marker='None')

        if (multi_obj_approach == 'rectangle'):
            new = list( idx_para[sort_idx[ee_masked[sort_idx,iobj] < cutoff_obj[iobj]]] )
            if (iobj == 0):
                keepit = new
            else:
                # intersection of parameters
                keepit = list(set(keepit) & set(new))
            print('Kept parameters (start with 0): ', astr(np.sort(keepit)+1),'   --> ',np.size(keepit),' parameters')

        if (not(noplot)):
            # threshold
            xmin, xmax = sub.get_xlim()
            line3 = plt.plot([xmin,xmax], [cutoff_obj[iobj],cutoff_obj[iobj]])
            plt.setp(line3, linestyle='--', linewidth=lwidth, color=lcol2, marker='None', label=str2tex('$\mu^*_{thres}$',usetex=usetex))

        if (not(noplot)):
            sub.text(1.02, 0.5, str2tex(obj_names[iobj],usetex=usetex),
                         rotation=90, fontsize='large',
                         horizontalalignment='left', verticalalignment='center',
                         transform=sub.transAxes)

        if (not(noplot)):
            xlim = [-1,dims]
            plt.setp(sub,  xlim=xlim)

    if (not(noplot)):
        # legend
        ll = plt.legend(frameon=frameon, ncol=6, bbox_to_anchor=(llxbbox,llybbox), loc='lower center',
                    scatterpoints=1, numpoints=1,
                   labelspacing=llrspace, columnspacing=llcspace, handletextpad=llhtextpad, handlelength=llhlength)
        plt.setp(ll.get_texts(), fontsize='small')

    if (not(noplot)):
        if (outtype == 'pdf'):
            pdf_pages.savefig(fig)
            plt.close(fig)
        elif (outtype == 'png'):
            pngfile = pngbase+"{0:04d}".format(ifig)+".png"
            fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
            plt.close(fig)

    # Write cutoff file
    splits=maskfile.split('/')
    ofile=''
    for i in range(0,len(splits)-1):
        ofile=ofile+splits[i]+'/'
    ofile=ofile+'cutoff_'+astr(iobj+1)+'.dat'
    f = open(ofile,'w')
    print('cutoff', file=f)
    if (type(cutoff_obj[iobj])==np.float) | (type(cutoff_obj[iobj])==np.float64):
        print(cutoff_obj[iobj], file=f)
    else:
        print(cutoff_obj[iobj][0], file=f)
    f.close()
    print("wrote:   '"+ofile+"'")

if (multi_obj_approach == 'triangle'):
    keepit = list(idx_para[np.where(np.sum(ee_masked[sort_idx]/cutoff_obj,axis=1)<1)[0]])
    keepit = list(idx_para[sort_idx[np.sum(ee_masked[sort_idx]/cutoff_obj,axis=1)<1]])
    print('')
    if keepit != []:
        print('Parameters for next iteration: ', astr(np.sort(keepit)+1),'   --> ',np.size(keepit),' parameters')
    else:
        print('Parameters for next iteration: ', np.sort(keepit)+1,'   --> ',np.size(keepit),' parameters')

# Write masked parameter file
ofile=maskfile+'.new'
#print('')
#print('Write ascii data ', ofile)
f = open(ofile,'w')
print('# para   dist       lower     upper     default   informative(0)_or_noninformative(1)', file=f)
print('#                   mean      stddev                                                 ', file=f)
for ii in range(mask_para.shape[0]):
    kk = '0'
    if ii in keepit: kk = '1'
    print(para_name[ii], para_dist[ii], lower_bound[ii], upper_bound[ii], initial[ii], kk, file=f)
f.close()
print("wrote:   '"+ofile+"'")

if (not(noplot)):
    if (outtype == 'pdf'):
        pdf_pages.close()
    elif (outtype == 'png'):
        pass
    else:
        plt.show()

