#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Defines common functions that are used in curve_fit or fmin parameter estimations.

    Defines the functions in two forms (ex. of 3 params):
        1. func(x, p1, p2, p3)
        2. func_p(x, p)  with p[0:3]
    These cost functions can be used for example with curve_fit
        p, cov  = opt.curve_fit(jams.functions.f1x, x, y, p0=[p0,p1])

    Defines also two cost functions, one with absolute sum, one with squared sum of deviations.
        3. cost_func    sum(abs(obs-func))
        4. cost2_func   sum((obs-func)**2)
    These cost functions can be used for example with fmin
        p = opt.fmin(jams.functions.cost_f1x, np.array([p1,p2]), args=(x,y), disp=False)
    or
        p, nfeval, rc = opt.fmin_tnc(jams.functions.cost_f1x, [p1,p2], bounds=[[None,None],[None,None]],
                                     args=(x,y), approx_grad=True, disp=False)

    Note the different argument orders:
    curvefit wants f(x,*args) with the independent variable as the first argument
             and the parameters to fit as separate remaining arguments.
    fmin is a general minimiser with respect to the first argument, i.e. func(p,*args).

    There are also two common cost functions (absolute and squared deviations) where any function
    in the form func(x, p) can be used as second argument:
        5. cost_abs(p, func, x, y)
        6. cost_square(p, func, x, y)
    Used for example as
        p = opt.fmin(jams.functions.cost_abs, np.array([p1,p2]), args=(jams.functions.f1x_p,x,y), disp=False)
    or
        p, nfeval, rc = opt.fmin_tnc(jams.functions.cost_square, [p1,p2], bounds=[[None,None],[None,None]],
                                     args=(jams.functions.f1x_p,x,y), approx_grad=True, disp=False)


    Definition
    ----------
    Current functions are (there is always the second form with the name appended by _p;
    these are used in the cost functions).

        arrhenius         1 param:  Arrhenius temperature dependence of biochemical rates:
                                    exp((T-TC25)*E/(T25*R*(T+T0))), parameter: E
        f1x               2 params: General 1/x function: a + b/x
        fexp              3 params: General exponential function: a + b * exp(c*x)
        gauss             2 params: Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2)), parameter: mu, sig
        lasslop           6 params: Lasslop et al. (2010) a rectangular, hyperbolic light-response GPP
                                    with Lloyd & Taylor (1994) respiration and the maximum canopy uptake
                                    rate at light saturation decreases exponentially with VPD as in Koerner (1995)
        line0             1 params: Straight line: a*x
        line              2 params: Straight line: a + b*x
        lloyd_fix         2 params: Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC
        lloyd_only_rref   1 param:  Lloyd & Taylor (1994) Arrhenius type with fixed exponential term
        logistic          3 params: Logistic function: a/(1+exp(-b(x-c)))
        logistic_offset   4 params: Logistic function with offset: a/(1+exp(-b(x-c))) + d
        logistic2_offset  7 params: Double logistic function with offset L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a2

        poly              n params: General polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        sabx              2 params: sqrt(f1x), i.e. general sqrt(1/x) function: sqrt(a + b/x)
        see               3 params: Sequential Elementary Effects fitting function: a*(x-b)**c


    Input / Output
    --------------
    See the help of the individual functions for explanations of in/out, etc.


    Examples
    --------
    >>> Rref = 1.0
    >>> E0   = 126.
    >>> T    = 293.15
    >>> resp = 2.0
    >>> from autostring import astr
    >>> print(astr(lloyd_fix(T, Rref, E0),3,pp=True))
    1.406
    >>> print(astr(lloyd_fix_p(T, [Rref, E0]),3,pp=True))
    1.406
    >>> print(astr(cost_lloyd_fix([Rref, E0], T, resp),3,pp=True))
    0.594
    >>> print(astr(cost2_lloyd_fix([Rref, E0], T, resp),3,pp=True))
    0.353

    >>> print(astr(poly(T,2,1),3,pp=True))
    295.150
    >>> print(astr(poly_p(T,[2,1]),3,pp=True))
    295.150


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT License.

    Copyright (c) 2012-2015 Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  MC, Dec 2012
    Modified, MC, Feb 2013 - ported to Python 3
              MC, May 2013 - general cost function cost_abs, cost_square
              MC, Oct 2013 - test functions such as Rosenbrock, Griewank, etc.
              MC, Feb 2014 - line0
              MC, Mar 2015 - separate file
"""
import numpy as np
import scipy.special as sp
from general_functions import logistic_p, logistic_offset_p, logistic2_offset_p

__all__ = ['cost_abs', 'cost_square',
           'f1x', 'f1x_p', 'cost_f1x', 'cost2_f1x',
           'fexp', 'fexp_p', 'cost_fexp', 'cost2_fexp',
           'gauss', 'gauss_p', 'cost_gauss', 'cost2_gauss',
           'lasslop', 'lasslop_p', 'cost_lasslop', 'cost2_lasslop',
           'line', 'line_p', 'cost_line', 'cost2_line',
           'line0', 'line0_p', 'cost_line0', 'cost2_line0',
           'lloyd_fix', 'lloyd_fix_p', 'cost_lloyd_fix', 'cost2_lloyd_fix',
           'lloyd_only_rref', 'lloyd_only_rref_p', 'cost_lloyd_only_rref', 'cost2_lloyd_only_rref',
           'sabx', 'sabx_p', 'cost_sabx', 'cost2_sabx',
           'poly', 'poly_p', 'cost_poly', 'cost2_poly',
           'cost_logistic', 'cost2_logistic',
           'cost_logistic_offset', 'cost2_logistic_offset',
           'cost_logistic2_offset', 'cost2_logistic2_offset',
           'see', 'see_p', 'cost_see', 'cost2_see']

# -----------------------------------------------------------
# general cost functions
def cost_abs(p, func, x, y):
    """ General cost function for robust optimising func(p, x) vs y with sum of absolute deviations"""
    return np.sum(np.abs(y-func(x,p)))

def cost_square(p, func, x, y):
    """ General cost function for least square optimising func(p, x) vs y"""
    return np.sum((y-func(x,p))**2)


# -----------------------------------------------------------
# a+b/x
def f1x(x,a,b):
  ''' General 1/x function: a + b/x
        x    independent variable
        a    1. parameter
        b    2. parameter
  '''
  return a+b/x

def f1x_p(x,p):
  ''' General 1/x function: a + b/x
        x    independent variable
        p    array of size 2, parameters
  '''
  return p[0]+p[1]/x

def cost_f1x(p,x,y):
  ''' Sum of absolut errors between obs and general 1/x function: a + b/x
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-f1x_p(x,p)))

def cost2_f1x(p,x,y):
  ''' Sum of squared errors between obs and general 1/x function: a + b/x
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-f1x_p(x,p))**2)


# -----------------------------------------------------------
# a+b*exp(c*x)
def fexp(x,a,b,c):
  ''' General exponential function: a + b * exp(c*x)
        x    independent variable in exp
        a    1. parameter
        b    2. parameter
        c    3. parameter
  '''
  return a+b*np.exp(c*x)

def fexp_p(x,p):
  ''' General exponential function: a + b * exp(c*x)
        x    independent variable in exp
        p    array of size 3, parameters
  '''
  return p[0]+p[1]*np.exp(p[2]*x)

def cost_fexp(p,x,y):
  ''' Sum of absolut errors between obs and general exponential function: a + b * exp(c*x)
        p    array of size 3, parameters
        x    independent variable in exp
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-fexp_p(x,p)))

def cost2_fexp(p,x,y):
  ''' Sum of squared errors between obs and general exponential function: a + b * exp(c*x)
        p    array of size 3, parameters
        x    independent variable in exp
        y    dependent variable to optimise
  '''
  return np.sum((y-fexp_p(x,p))**2)


# -----------------------------------------------------------
# Gauss: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
def gauss(x,mu,sig):
  ''' Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
        x    independent variable
        mu   mean
        sig  width
  '''
  return np.exp(-(x-mu)**2/(2.*sig**2))/(sig*np.sqrt(2.*np.pi))

def gauss_p(x,p):
  ''' Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
        x    independent variable in exp
        p    array of size 2, parameters mu and sig
  '''
  return np.exp(-(x-p[0])**2/(2.*p[1]**2))/(p[1]*np.sqrt(2.*np.pi))

def cost_gauss(p,x,y):
  ''' Sum of absolut errors between obs and Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
        p    array of size 2, parameters mu and sig
        x    independent variable in exp
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-gauss_p(x,p)))

def cost2_gauss(p,x,y):
  ''' Sum of squared errors between obs and Gauss function: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
        p    array of size 2, parameters mu and sig
        x    independent variable in exp
        y    dependent variable to optimise
  '''
  return np.sum((y-gauss_p(x,p))**2)

# -----------------------------------------------------------
# lasslop
def lasslop(Rg, et, VPD, alpha, beta0, k, Rref):
    """ Lasslop et al. (2010) is basically the rectangular, hyperbolic
        light-response of NEE as by Falge et al. (2001), where the
        respiration is calculated with Lloyd & Taylor (1994), and the
        maximum canopy uptake rate at light saturation decreases
        exponentially with VPD as in Koerner (1995).
        Rg      Global radiation [W m-2]
        et      Exponential in Lloyd & Taylor: np.exp(E0*(1./(Tref-T0)-1./(T-T0))) []
        VPD     Vapour Pressure Deficit [Pa]
        alpha   Light use efficiency, i.e. initial slope of light response curve [umol(C) J-1]
        beta0   Maximum CO2 uptake rate at VPD0=10 hPa [umol(C) m-2 s-1]
        k       e-folding of exponential decrease of maximum CO2 uptake with VPD increase [Pa-1]
        Rref    Respiration at Tref (10 degC) [umol(C) m-2 s-1]
    """
    # Lloyd & Taylor (1994)
    gamma = Rref*et
    # Koerner (1995)
    VPD0  = 1000. # 10 hPa
    kk    = np.maximum(np.minimum(-k*(VPD-VPD0), 600.), -600.)
    beta  = np.where(VPD > VPD0, beta0*np.exp(kk), beta0)
    return -alpha*beta*Rg/(alpha*Rg+beta) + gamma

def lasslop_p(Rg, p):
    """ Lasslop et al. (2010) is basically the rectangular, hyperbolic
        light-response of NEE as by Falge et al. (2001), where the
        respiration is calculated with Lloyd & Taylor (1994), and the
        maximum canopy uptake rate at light saturation decreases
        exponentially with VPD as in Koerner (1995).
        Rg     Global radiation [W m-2]
        p[0]   Exponential in Lloyd & Taylor: np.exp(E0*(1./(Tref-T0)-1./(T-T0))) []
        p[1]   Vapour Pressure Deficit [Pa]
        p[2]   Light use efficiency, i.e. initial slope of light response curve [umol(C) J-1]
        p[3]   Maximum CO2 uptake rate at VPD0=10 hPa [umol(C) m-2 s-1]
        p[4]   e-folding of exponential decrease of maximum CO2 uptake with VPD increase [Pa-1]
        p[5]   Respiration at Tref (10 degC) [umol(C) m-2 s-1]
    """
    # Lloyd & Taylor (1994)
    gamma = p[5]*p[0]
    # Koerner (1995)
    VPD0  = 1000. # 10 hPa
    kk    = np.maximum(np.minimum(-p[4]*(p[1]-VPD0), 600.), -600.)
    beta  = np.where(p[1] > VPD0, p[3]*np.exp(kk), p[3])
    return -p[2]*beta*Rg/(p[2]*Rg+beta) + gamma

def cost_lasslop(p, Rg, et, VPD, NEE):
    """ Cost function for Lasslop with sum of absolute deviations """
    return np.sum(np.abs(NEE-lasslop(Rg, et, VPD, p[0], p[1], p[2], p[3])))

def cost2_lasslop(p, Rg, et, VPD, NEE):
    """ Cost function for Lasslop with sum of squared deviations """
    return np.sum((NEE-lasslop(Rg, et, VPD, p[0], p[1], p[2], p[3]))**2)


# -----------------------------------------------------------
# a+b*x
def line(x,a,b):
  ''' Straight line: a + b*x
        x    independent variable
        a    1. parameter
        b    2. parameter
  '''
  return a+b*x

def line_p(x,p):
  ''' Straight line: a + b*x
        x    independent variable
        p    array of size 2, parameters
  '''
  return p[0]+p[1]*x

def cost_line(p,x,y):
  ''' Sum of absolut errors between obs and straight line: a + b*x
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-line_p(x,p)))

def cost2_line(p,x,y):
  ''' Sum of squared errors between obs and straight line: a + b*x
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-line_p(x,p))**2)


# -----------------------------------------------------------
# b*x
def line0(x,a):
  ''' Straight line through origin: a*x
        x    independent variable
        a    parameter
  '''
  return a*x

def line0_p(x,p):
  ''' Straight line through origin: p[0]*x
        x    independent variable
        p    parameter
  '''
  return p*x

def cost_line0(p,x,y):
  ''' Sum of absolut errors between obs and straight line though origin: a*x
        p    parameter
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-line0_p(x,p)))

def cost2_line0(p,x,y):
  ''' Sum of squared errors between obs and straight line through origin: a*x
        p    parameter
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-line0_p(x,p))**2)


# -----------------------------------------------------------
# lloyd_fix
def lloyd_fix(T, Rref, E0):
    """ Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC
        T       Temperature [k]
        Rref    Respiration at Tref=10 degC [umol(C) m-2 s-1]
        E0      Activation energy [K]
    """
    Tref = 283.15 #  10    [degC]
    T0   = 227.13 # -46.02 [degC]
    return Rref*np.exp(E0*(1./(Tref-T0)-1./(T-T0)))

def lloyd_fix_p(T, p):
    """ Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC
        T       Temperature [k]
        p[0]    Respiration at Tref=10 degC [umol(C) m-2 s-1]
        p[1]    Activation energy [K]
    """
    Tref = 283.15 #  10    [degC]
    T0   = 227.13 # -46.02 [degC]
    return p[0]*np.exp(p[1]*(1./(Tref-T0)-1./(T-T0)))

def cost_lloyd_fix(p, T, resp):
    """ Cost function for Lloyd with sum of absolute deviations """
    return np.sum(np.abs(resp-lloyd_fix_p(T,p)))

def cost2_lloyd_fix(p, T, resp):
    """ Cost function for Lloyd with sum of squared deviations """
    return np.sum((resp-lloyd_fix_p(T,p))**2)


# -----------------------------------------------------------
# lloyd_only_rref
def lloyd_only_rref(et, Rref):
    """ If E0 is know in Lloyd & Taylor (1994) then one can calc
        the exponential term outside the routine and the fitting
        becomes linear.
        et      exp-term in Lloyd & Taylor
        Rref    Respiration at Tref=10 degC [umol(C) m-2 s-1]
    """
    return Rref*et

def lloyd_only_rref_p(et, p):
    """ If E0 is know in Lloyd & Taylor (1994) then one can calc
        the exponential term outside the routine and the fitting
        becomes linear.
        et      exp-term in Lloyd & Taylor
        p[0]    Respiration at Tref=10 degC [umol(C) m-2 s-1]
    """
    return p[0]*et

def cost_lloyd_only_rref(p, et, resp):
    """ Cost function for rref with sum of absolute deviations """
    return np.sum(np.abs(resp-lloyd_only_rref_p(et,p)))

def cost2_lloyd_only_rref(p, et, resp):
    """ Cost function for rref  with sum of squared deviations """
    return np.sum((resp-lloyd_only_rref_p(et,p))**2)


# -----------------------------------------------------------
# sqrt(a + b/x) - theoretical form of Jackknife-after-bootstrap

def sabx(x, a, b):
  """ sqrt(a + b/x)
        a, b  parameters
        x     independent variable
  """
  return np.sqrt(a+b/x)

def sabx_p(x, p):
  """ sqrt(p[0] + p[1]/x)
        p    array of size 2, parameters
        x    independent variable
  """
  return np.sqrt(p[0]+p[1]/x)

def cost_sabx(p,x,y):
  ''' Cost function for sqrt of general 1/x function with sum of absolute deviations
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-sabx_p(x,p)))

def cost2_sabx(p,x,y):
  ''' Cost function for sqrt of general 1/x function with sum of squared deviations
        p    array of size 2, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-sabx_p(x,p))**2)


# -----------------------------------------------------------
# c0 + c1*x + c2*x**2 + ... + cn*x**n
def poly(x,*args):
  ''' General polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        x     independent variable
        c0    1. parameter
        c1    2. parameter
        ...
  '''
  return np.polynomial.polynomial.polyval(x, list(args))

def poly_p(x,p):
  ''' General polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        x    independent variable
        p    array of size n, parameters
  '''
  return np.polynomial.polynomial.polyval(x, p)

def cost_poly(p,x,y):
  ''' Sum of absolut errors between obs and general polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        p    array of size n, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-poly_p(x,p)))

def cost2_poly(p,x,y):
  ''' Sum of squared errors between obs and general polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n
        p    array of size n, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-poly_p(x,p))**2)


# -----------------------------------------------------------
# a/(1+exp(-b(x-c))) - logistic function
def cost_logistic(p, x, y):
  ''' Cost function for logistic function fitting function with sum of absolute deviations
        p    array of size 3, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-logistic_p(x,p)))

def cost2_logistic(p,x,y):
  ''' Cost function for logistic function fitting function with sum of squared deviations
        p    array of size 3, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-logistic_p(x,p))**2)


# -----------------------------------------------------------
# a/(1+exp(-b(x-c))) + d - logistic function with offset
def cost_logistic_offset(p, x, y):
  ''' Cost function for logistic function with offset fitting function with sum of absolute deviations
        p    4D-array of parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-logistic_offset_p(x,p)))

def cost2_logistic_offset(p,x,y):
  ''' Cost function for logistic function with offset fitting function with sum of squared deviations
        p    4D-array of parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-logistic_offset_p(x,p))**2)


# -----------------------------------------------------------
# L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a2 - double logistic function with offset
def cost_logistic2_offset(p, x, y):
  ''' Cost function for double logistic function with offset fitting function with sum of absolute deviations
        p    4D-array of parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-logistic2_offset_p(x,p)))

def cost2_logistic2_offset(p,x,y):
  ''' Cost function for double logistic function with offset fitting function with sum of squared deviations
        p    4D-array of parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-logistic2_offset_p(x,p))**2)


# -----------------------------------------------------------
# a*(x-b)**c - Sequential Elementary Effects fitting function
def see(x, a, b, c):
  """ a*(x-b)**c
        a, b, c  parameters
        x        independent variable
  """
  return np.where((x-b)<0., 0., a*(x-b)**c)

def see_p(x, p):
  """ p[0]*(x-p[1])**p[2]
        a, b, c  parameters
        x        independent variable
  """
  return np.where((x-p[1]) < 0., 0., p[0] * (x-p[1])**p[2])

def cost_see(p, x, y):
  ''' Cost function for Sequential Elementary Effects fitting function with sum of absolute deviations
        p    array of size 3, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum(np.abs(y-see_p(x,p)))

def cost2_see(p,x,y):
  ''' Cost function for Sequential Elementary Effects fitting function with sum of squared deviations
        p    array of size 3, parameters
        x    independent variable
        y    dependent variable to optimise
  '''
  return np.sum((y-see_p(x,p))**2)


# -----------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # Rref = 1.0
    # E0   = 126.
    # T    = 293.15
    # resp = 2.0
    # print(lloyd_fix(T, Rref, E0))
    # #1.40590910521
    # print(lloyd_fix_p(T, [Rref, E0]))
    # #1.40590910521
    # print(cost_lloyd_fix([Rref, E0], T, resp))
    # #0.59409089479
    # print(cost2_lloyd_fix([Rref, E0], T, resp))
    # #0.352943991272

    # print(poly(T,2,1))
    # #295.15
    # print(poly_p(T,[2,1]))
    # #295.15
