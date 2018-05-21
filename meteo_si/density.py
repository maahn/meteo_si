# -*- coding: utf-8 -*-
# (c) Jan Schween 2005 (gnuplot) 
# (c) Mario Mech 2009 (python)
# (c) Maximilian Maahn 2011 (python)



from __future__ import absolute_import, division, print_function
import numpy as np

from .due import due, Doi

from .constants import *
from .humidity import rh2q
# from .temperature import *


__all__ = ["moist_rho_rh", "moist_rho_q"]



def moist_rho_rh(p,T,rh,*qm):
	'''
	Input:
	p is in Pa
	T is in K
	rh is in Pa/Pa
	Optional, several possible:
	qm is in kg/kg other species which contribute to the air mass! (ice, snow, cloud etc.)

	Output:
	density of moist air [kg/m^3]

	Example:
	moist_rho_rh(p,T,rh,q_ice,q_snow,q_rain,q_cloud,q_graupel,q_hail)

	
	'''
	if np.any(rh > 5): raise ValueError("rh must not be in %")
	
	q = rh2q(rh,T,p)
	
	return moist_rho_q(p,T,q,*qm)

def moist_rho_q(p,T,q,*qm):
	'''
	Input p is in Pa
	T is in K
	q is in kg/kg
	Optional, several possible:
	qm is in kg/kg other species which contribute to the air mass! (ice, snow, cloud etc.)
	
	Output:
	density of moist air [kg/m^3]
	
	Example:
	moist_rho_q(p,T,q,q_ice,q_snow,q_rain,q_cloud,q_graupel,q_hail)
	'''

	if len(qm)> 0: 
		#get rid of masked data!
		qm[qm<0] = 0
		qm = np.sum(qm,axis=0)
	else: 
		qm = 0
	
	moist_rho_q = p/(Rair*T*(1+(Rvapor/Rair-1)*q-qm))
	

	
	if np.any(moist_rho_q < 0):
		if np.any(moist_rho_q < -0.001): 
			raise ValueError("calculated negative densities!")
		else:
			try: moist_rho_q[moist_rho_q<0] = 0
			except: moist_rho_q = 0
			
	
	return moist_rho_q