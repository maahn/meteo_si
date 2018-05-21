# -*- coding: utf-8 -*-
# (c) Jan Schween 2005 (gnuplot)
# (c) Mario Mech 2009 (python)
# (c) Maximilian Maahn 2011 (python)

from __future__ import absolute_import, division, print_function
import numpy as np

# from .due import due, Doi

from .constants import Mwml, Rvapor
from .density import moist_rho_q
from . import temperature


__all__ = ["e2q", "q2e", "rh2q", "rh2a", "rh_to_iwv", "e_sat_gg_ice",
           "e_sat_gg_water", "q2rh", "a2rh"]


def e_sat_gg_water(T):
    '''
    Calculates the saturation pressure over water after "Guide to
    Meteorological Instruments and Methods of Observation" (CIMO Guide)
    (WMO, 2008).

    Input:
    T in Kelvin.
    Output:
    e_sat_gg_water in Pa.
    '''
    T = temperature.kelvin_2_celsius(T)
    e_sat_gg_water = 100 * 6.112 * np.exp(17.62 * T/(243.12 + T))
    return e_sat_gg_water


def e_sat_gg_ice(T):
    '''
    Calculates the saturation pressure over water after "Guide to
    Meteorological Instruments and Methods of Observation" (CIMO Guide)
    (WMO, 2008).

    Input:
    T in Kelvin.
    Output:
    e_sat_gg_ice in Pa.
    '''
    T = temperature.kelvin_2_celsius(T)
    e_sat_gg_ice = 100 * 6.112 * np.exp(22.46 * T/(272.62 + T))
    return e_sat_gg_ice


def e2q(e, p):
    '''
    Calculate the specific humidity from water vapour pressure and air
    pressure.

    Input:
    e is in Pa
    p is in Pa

    Output
    q in kg/kg
    '''
    q = Mwml*e/(p-(1-Mwml)*e)
    return q


def q2e(q, p):
    '''
    Calculate water vapour pressure from the specific humidity and air
    pressure.

    Input:
    q in kg/kg
    p is in Pa

    Output
    e is in Pa
    '''
    e = p/((Mwml/q)+1-Mwml)
    return e


def rh2q(rh, T, p, e_sat_func=e_sat_gg_water):
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
    with np.errstate(divide='ignore', invalid='ignore'):
        if np.any(rh > 5):
            raise TypeError("rh must not be in %")

    eStar = e_sat_func(T)
    e = rh*eStar
    q = e2q(e, p)
    del e, eStar
    return q


def rh2a(rh, T, e_sat_func=e_sat_gg_water):
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

    with np.errstate(divide='ignore', invalid='ignore'):
        if np.any(rh > 5):
            raise TypeError("rh must not be in %")

    e = rh*e_sat_func(T)
    a = e/(Rvapor*T)
    return a


def a2rh(a, T, e_sat_func=e_sat_gg_water):
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
        rh = e/e_sat_func(T)
        return rh


def q2rh(q, T, p, e_sat_func=e_sat_gg_water):
    '''
    Calculate relative humidity from specific humidity

    Input:
    T is in K
    p is in Pa
    q in kg/kg

    Output:
    rh is in Pa/Pa
    '''

    # if neAvail: e = ne.evaluate("p/(Mwml*((1/q)+(1/(Mwml)-1)))")
    e = p/(Mwml*((1/q)+(1/(Mwml)-1)))

    eStar = e_sat_func(T)
    rh = e/eStar
    return rh


def rh_to_iwv(relhum_lev, temp_lev, press_lev, hgt_lev):
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
    dz = np.diff(hgt_lev, axis=-1)
    relhum = (relhum_lev[..., 0:-1] + relhum_lev[..., 1:])/2.
    temp = (temp_lev[..., 0:-1] + temp_lev[..., 1:])/2.

    xp = -1.*np.log(press_lev[..., 1:]/press_lev[..., 0:-1])/dz
    press = -1.*press_lev[..., 0:-1]/xp*(np.exp(-xp*dz)-1.)/dz

    q = rh2q(relhum, temp, press)
    rho_moist = moist_rho_q(press, temp, q)

    return np.sum(q*rho_moist*dz)
