# -*- coding: utf-8 -*-
# (c) Jan Schween 2005 (gnuplot) 
# (c) Mario Mech 2009 (python)
# (c) Maximilian Maahn 2011 (python)



from __future__ import absolute_import, division, print_function
import numpy as np

from .due import due, Doi

from .constants import *
from .temperature import *
from .density import *




def e2q(e,p):
	'''
	Calculate the specific humidity from water vapour pressure and air pressure. 

	Input:
	e is in Pa
	p is in Pa
	
	Output
	q in kg/kg
	'''
	q = Mwml*e/(p-(1-Mwml)*e)
	return q
  
def q2e(q,p):
	'''
	Calculate water vapour pressure from the specific humidity and air pressure. 

	Input:
	q in kg/kg
	p is in Pa
	
	Output
	e is in Pa
	'''
	e=p/((Mwml/q)+1-Mwml)
	return e

def rh2q(rh,T,p):
	'''
	Calculate the specific humidity from relative humidity, air temperature,
	and pressure. 

	Input:
	T is in K
	rh is in Pa/Pa
	p is in Pa
	
	Output
	q in kg/kg
	'''
	if np.any(rh > 5): raise TypeError("rh must not be in %")
	
	eStar = e_sat_gg_water(T)
	e = rh*eStar
	q = e2q(e,p)
	del e, eStar
	return q

	
def rh2a(rh,T):
	'''
	Calculate the absolute humidity from relative humidity, air temperature,
	and pressure.

	Input T is in K
	rh is in Pa/Pa
	p is in Pa
	Output 
	a in kg/m^3
	
	Source: Kraus: Chapter 8.1.2
	'''
	
	if np.any(rh > 5): raise TypeError("rh must not be in %")

	e = rh*e_sat_gg_water(T)
	a = e/(Rvapor*T)
	return a

def a2rh(a,T):
        '''
        Calculate the relative from absolute humidity and air temperature.

        Input
        T is in K
        a in kg/kg
        Output
        rh in Pa/Pa
        
        Source: Kraus: Chapter 8.1.2
        '''
        
        e = a*(Rvapor*T)
        rh = e/e_sat_gg_water(T)
        return rh
  
	
	
def q2rh(q,T,p):
	'''
	Calculate relative humidity from specific humidity
	
	Input:
	T is in K
	p is in Pa
	q in kg/kg
	
	Output:
	rh is in Pa/Pa
	'''
	
	
	if neAvail: e = ne.evaluate("p/(Mwml*((1/q)+(1/(Mwml)-1)))")
	else: e = p/(Mwml*((1/q)+(1/(Mwml)-1)))
	
	eStar = e_sat_gg_water(T)
	rh = e/eStar
	del e,eStar
	return rh
	
def e_sat_gg_water(T):
	'''
	Calculates the saturation pressure over water after Goff and Gratch (1946).
	It is the most accurate that you can get for a temperture range from -90°C to +80°C.
	Source: Smithsonian Tables 1984, after Goff and Gratch 1946
	http://cires.colorado.edu/~voemel/vp.html
	http://hurri.kean.edu/~yoh/calculations/satvap/satvap.html

	Input:
	T in Kelvin.
	Output:
	e_sat_gg_water in Pa.
	'''
	if neAvail: e_sat_gg_water = ne.evaluate("100 * 1013.246 * 10**( -7.90298*(373.16/T-1) + 5.02808*log10(373.16/T) - 1.3816e-7*(10**(11.344*(1-T/373.16))-1) + 8.1328e-3 * (10**(-3.49149*(373.16/T-1))-1) )")
	else: e_sat_gg_water = 100 * 1013.246 * 10**( -7.90298*(373.16/T-1) + 5.02808*np.log10(373.16/T) - 1.3816e-7*(10**(11.344*(1-T/373.16))-1) + 8.1328e-3 * (10**(-3.49149*(373.16/T-1))-1) )
	return e_sat_gg_water

def rh_to_iwv(relhum_lev,temp_lev,press_lev,hgt_lev):
	'''
	Calculate the integrated water vapour

	Input:
	T is in K
	rh is in Pa/Pa
	p is in Pa
	z is in m
	
	Output
	iwv in kg/m^2
	'''
	dz = np.diff(hgt_lev,axis=-1)
	relhum = (relhum_lev[...,0:-1] + relhum_lev[...,1:])/2.
	temp = (temp_lev[...,0:-1] + temp_lev[...,1:])/2.

	xp = -1.*np.log(press_lev[...,1:]/press_lev[...,0:-1])/dz
	press = -1.*press_lev[...,0:-1]/xp*(exp(-xp*dz)-1.)/dz

	q = meteoSI.rh2q(relhum,temp,press)
	rho_moist = meteoSI.moist_rho_q(press,temp,q)

	return np.sum(q*rho_moist*dz)
