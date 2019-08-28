#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Defines general functions, derivatives, etc.


    Definition
    ----------
    Current functions are:

    curvature             Curvature of function f: f''/(1+f'^2)^3/2
    logistic              logistic function L/(1+exp(-k(x-x0)))
    logistic_p
    dlogistic             First derivative of logistic function
    d2logistic            Second derivative of logistic function
    logistic_offset       logistic function with offset L/(1+exp(-k(x-x0))) + a
    logistic_offset_p
    dlogistic_offset      First derivative of logistic function with offset
    d2logistic_offset     Second derivative of logistic function with offset
    logistic2_offset      Double logistic function with offset L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a2
    logistic2_offset_p
    dlogistic2_offset     First derivative of double logistic function with offset
    d2logistic2_offset    Second derivative of double logistic function with offset


    Input / Output
    --------------
    See the help of the individual functions for explanations of in/out, etc.


    Examples
    --------
    ToDo.


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT License.

    Copyright (c) 2015-2017 Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  MC, Mar 2015
    Modified, MC, Dec 2017 - logistic_p, logistic_offset_p
"""
import numpy as np
import scipy.special as sp

__all__ = ['curvature',
           'logistic', 'dlogistic', 'd2logistic', 'logistic_p',
           'logistic_offset', 'dlogistic_offset', 'd2logistic_offset', 'logistic_offset_p',
           'logistic2_offset', 'dlogistic2_offset', 'd2logistic2_offset', 'logistic2_offset_p']

# -----------------------------------------------------------
# curvature of function
def curvature(x, dfunc, d2func, *args, **kwargs):
    """ Curvature of function f''/(1+f'^2)^3/2
          x         independent variable
          dfunc     first derivative of function f: f'
          d2func    second derivative of function f: f''
          args      arguments for dfunc and d2func
          kwargs    keyword arguments for dfunc and d2func
    """
    return d2func(x, *args, **kwargs)/(1.+dfunc(x, *args, **kwargs)**2)**1.5

# -----------------------------------------------------------
# a/(1+exp(-b(x-c))) - logistic function
def logistic(x, L, k, x0):
    """ logistic function L/(1+exp(-k(x-x0)))
          x         independent variable
          L         maximum
          k         steepness
          x0        inflection point
    """
    return L*sp.expit(k*(x-x0))

def logistic_p(x, p):
  """ logistic function p[0]/(1+exp(-p[1](x-p[2])))
        x        independent variable
        p        array of size 3, parameters
  """
  return logistic(x, p[0], p[1], p[2])

# -----------------------------------------------------------
# 1st derivative of logistic functions
def dlogistic(x, L, k, x0):
    """ First derivative of logistic function L/(1+exp(-k(x-x0)))
          x         independent variable
          L         maximum
          k         steepness
          x0        inflection point
    """
    return k*L/(2.*(np.cosh(k*(x-x0))+1.))

# -----------------------------------------------------------
# 2nd derivative of logistic functions
def d2logistic(x, L, k, x0):
    """ Second derivative of logistic function L/(1+exp(-k(x-x0)))
          x         independent variable
          L         maximum
          k         steepness
          x0        inflection point
    """
    return -k**2 * L * np.sinh(k*(x-x0))/(2.*(np.cosh(k*(x-x0))+1.)**2)

# -----------------------------------------------------------
# L/(1+exp(-k(x-x0))) + a - logistic function with offset
def logistic_offset(x, L, k, x0, a):
    """ logistic function with offset L/(1+exp(-k(x-x0))) + a
          x         independent variable
          L         maximum
          k         steepness
          x0        inflection point
          a         offset
    """
    return L*sp.expit(k*(x-x0)) + a

def logistic_offset_p(x, p):
  """ logistic function with offset p[0]/(1+exp(-p[1](x-p[2]))) + p[3]
        x    independent variable
        p    4D-array of parameters
  """
  return logistic_offset(x, p[0], p[1], p[2], p[3])

# -----------------------------------------------------------
# 1st derivative of logistic functions with offset
def dlogistic_offset(x, L, k, x0, a):
    """ First derivative of logistic function L/(1+exp(-k(x-x0))) + a
          x         independent variable
          L         maximum
          k         steepness
          x0        inflection point
          a         offset
    """
    return k*L/(2.*(np.cosh(k*(x-x0))+1.))

# -----------------------------------------------------------
# 2nd derivative of logistic functions with offset
def d2logistic_offset(x, L, k, x0, a):
    """ Second derivative of logistic function L/(1+exp(-k(x-x0))) + a
          x         independent variable
          L         maximum
          k         steepness
          x0        inflection point
          a         offset
    """
    return -k**2 * L * np.sinh(k*(x-x0))/(2.*(np.cosh(k*(x-x0))+1.)**2)

# -----------------------------------------------------------
# L/(1+exp(-k(x-x0))) + a - logistic function with offset
def logistic2_offset(x, L1, k1, x01, L2, k2, x02, a):
    """ double logistic function with offset L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a2
          x         independent variable
          L1        maximum 1st logistic function
          k1        steepness 1st logistic function
          x01       inflection point 1st logistic function
          L2        maximum 2nd logistic function
          k2        steepness 2nd logistic function
          x02       inflection point 2nd logistic function
          a         offset
    """
    return L1*sp.expit(k1*(x-x01)) - L2*sp.expit(k2*(x-x02)) + a

def logistic2_offset_p(x, p):
  """ double logistic function with offset p[0]/(1+exp(-p[1](x-p[2]))) - p[3]/(1+exp(-p[4](x-p[5]))) + p[6]
        x    independent variable
        p    4D-array of parameters
  """
  return logistic2_offset(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6])

# -----------------------------------------------------------
# 1st derivative of logistic functions with offset
def dlogistic2_offset(x, L1, k1, x01, L2, k2, x02, a):
    """ First derivative of double logistic function L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a2
          x         independent variable
          L1        maximum 1st logistic function
          k1        steepness 1st logistic function
          x01       inflection point 1st logistic function
          L2        maximum 2nd logistic function
          k2        steepness 2nd logistic function
          x02       inflection point 2nd logistic function
          a         offset
    """
    return ( k1*L1/(2.*(np.cosh(k1*(x-x01))+1.)) -
             k2*L2/(2.*(np.cosh(k2*(x-x02))+1.)) )

# -----------------------------------------------------------
# 2nd derivative of logistic functions with offset
def d2logistic2_offset(x, L1, k1, x01, L2, k2, x02, a):
    """ Second derivative of logistic function L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a2
          x         independent variable
          L1        maximum 1st logistic function
          k1        steepness 1st logistic function
          x01       inflection point 1st logistic function
          L2        maximum 2nd logistic function
          k2        steepness 2nd logistic function
          x02       inflection point 2nd logistic function
          a         offset
    """
    return ( -k1**2 * L1 * np.sinh(k1*(x-x01))/(2.*(np.cosh(k1*(x-x01))+1.)**2)
             +k2**2 * L2 * np.sinh(k2*(x-x02))/(2.*(np.cosh(k2*(x-x02))+1.)**2) )

# -----------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
