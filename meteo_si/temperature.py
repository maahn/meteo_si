# -*- coding: utf-8 -*-
# (c) Jan Schween 2005 (gnuplot) 
# (c) Mario Mech 2009 (python)
# (c) Maximilian Maahn 2011 (python)



from __future__ import absolute_import, division, print_function
import numpy as np

# from .due import due, Doi

from .constants import *
from . import humidity



__all__ = ["kelvin_2_celsius", "celsius_to_kelvin","T_virt_rh","T_virt_q"]


def kelvin_2_celsius(T):
	'''
	Calculate temperature in Celsius
	'''

	return T + Tnull

def celsius_to_kelvin(C):
	'''
	Calculate temperature in Kelvin
	'''
	return C - Tnull


def T_virt_rh(T,rh,p):
  '''
  Calculate the virtual temperature from air temperature,
  pressure, and relative humidity.

  Input:
  T is in K
  rh is in Pa/Pa
  p is in Pa
  
  Output:
  T_virt in K
  '''
  if np.any(rh > 5): raise TypeError("rh must not be in %")
  return T_virt_q(T,humidity.rh2q(rh,T,p))

def T_virt_q(T,q):
    '''
    Calculate the virtual temperature from air temperature and specific humidity.

    Input:
    T is in K
    q is in kg/kg
    
    Output:
    T_virt in K
    '''
    return T + T * (Rvapor/Rair-1) * q