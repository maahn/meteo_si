# -*- coding: utf-8 -*-
# (c) Jan Schween 2005 (gnuplot)
# (c) Mario Mech 2009 (python)
# (c) Maximilian Maahn 2011 (python)

import numpy as np

# from .due import due, Doi

import meteo_si.constants
import meteo_si.temperature


__all__ = ['a2e', 'e2a', "e2q", "q2e", "rh2q", "rh2a", "rh_to_iwv",
           "e_sat_gg_ice", "e_sat_gg_water", "e_sat_goffgratch_water",
           "q2rh", "a2rh"]


def _reject_percent_rh(rh):
    """
    Guard against relative humidity accidentally passed in %, which would
    otherwise silently produce results off by a factor of 100. Values up to
    5 are allowed, i.e. this is a sanity check, not a physical bound.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        if np.any(rh > 5):
            raise TypeError("rh must not be in %")


def a2e(a, T):
    """
    Calculate water vapor pressure from the absolute humidity and air
    temperature.

    Parameters
    ----------
    a:
        absolute humidity [kg / m3]
    T:
        Temperature in K

    Returns
    -------

    float :
        vapor pressure [Pa]

    """
    e = a * T * meteo_si.constants.Rvapor

    return e


def e2a(e, T):
    """
    Calculate the absolute humidity from water vapor pressure and air
    temperature.

    Parameters
    ----------
    e:
        vapor pressure [Pa]
    T:
        Temperature in K

    Returns
    -------

    float :
        absolute humidity [kg / m3]

    """

    a = e / (T * meteo_si.constants.Rvapor)

    return a


def e_sat_gg_water(T):
    """
    Calculates the saturation pressure over water after "Guide to
    Meteorological Instruments and Methods of Observation" (CIMO Guide)
    (WMO, 2008).

    Parameters
    ----------
    T:
        Temperature in K

    Returns
    -------

    float :
        saturation pressure [Pa]

    """
    T = meteo_si.temperature.kelvin_2_celsius(T)
    return 100 * 6.112 * np.exp(17.62 * T / (243.12 + T))


def e_sat_gg_ice(T):
    """
    Calculates the saturation pressure over ice after "Guide to
    Meteorological Instruments and Methods of Observation" (CIMO Guide)
    (WMO, 2008).

    Parameters
    ----------
    T:
        Temperature in K

    Returns
    -------

    float :
        saturation pressure [Pa]

    """
    T = meteo_si.temperature.kelvin_2_celsius(T)
    return 100 * 6.112 * np.exp(22.46 * T / (272.62 + T))


def e_sat_goffgratch_water(T):
    """
    Calculates the saturation pressure over water after Goff and Gratch
    (1946). More accurate than :func:`e_sat_gg_water` over a wide
    temperature range (-90 degC to +80 degC), at the cost of a more
    expensive formula -- the two are not interchangeable results, just two
    different approximations of the same physical quantity.

    Source: Smithsonian Tables 1984, after Goff and Gratch 1946.

    Parameters
    ----------
    T:
        Temperature in K

    Returns
    -------

    float :
        saturation pressure [Pa]

    """
    return 100 * 1013.246 * 10**(
        -7.90298*(373.16/T-1)
        + 5.02808*np.log10(373.16/T)
        - 1.3816e-7*(10**(11.344*(1-T/373.16))-1)
        + 8.1328e-3*(10**(-3.49149*(373.16/T-1))-1)
    )


def e2q(e, p):
    """
    Calculate the specific humidity from vapor pressure and air
    pressure.

    Parameters
    ----------
    e:
        vapor pressure [Pa]
    p:
        pressure [Pa]

    Returns
    -------

    float :
        specific humidity [kg / kg]

    """
    q = meteo_si.constants.Mwml * e / (p - (1 - meteo_si.constants.Mwml) * e)
    return q


def q2e(q, p):
    """
    Calculate water vapor pressure from the specific humidity and air
    pressure.

    Parameters
    ----------
    q:
        specific humidity [kg / kg]
    p:
        pressure [Pa]

    Returns
    -------

    float :
        vapor pressure [Pa]

    """

    e = p / ((meteo_si.constants.Mwml / q)+1-meteo_si.constants.Mwml)
    return e


def rh2q(rh, T, p, e_sat_func=e_sat_gg_water):
    """
    Calculate the specific humidity from relative humidity, air temperature,
    and pressure.

    Parameters
    ----------
    rh:
        Relative humidity in Pa / Pa
    T:
        Temperature in K
    p:
        pressure [Pa]
    e_sat_func: func, optional
        Function to estimate the saturation pressure. E.g. e_sat_gg_water for
        water and e_sat_gg_ice for ice.

    Returns
    -------

    float :
        specific humidity [kg / kg]

    """
    _reject_percent_rh(rh)

    e = rh * e_sat_func(T)
    return e2q(e, p)


def rh2a(rh, T, e_sat_func=e_sat_gg_water):
    """
    Calculate the absolute humidity from relative humidity, air temperature,
    and pressure.

    Parameters
    ----------
    rh:
        Relative humidity in Pa / Pa
    T:
        Temperature in K
    e_sat_func: func, optional
        Function to estimate the saturation pressure. E.g. e_sat_gg_water for
        water and e_sat_gg_ice for ice.

    Returns
    -------

    float :
        absolute humidity [kg / m3]

    """

    _reject_percent_rh(rh)

    e = rh * e_sat_func(T)
    return e2a(e, T)


def a2rh(a, T, e_sat_func=e_sat_gg_water):
    """
    Calculate the relative humidity from absolute humidity and air
    temperature. Source: Kraus, 'Die Atmosphäre der Erde', Chapter 8.1.2

    Parameters
    ----------
    a:
        absolute humidity [kg / m3]
    T:
        Temperature in K
    e_sat_func: func, optional
        Function to estimate the saturation pressure. E.g. e_sat_gg_water for
        water and e_sat_gg_ice for ice.

    Returns
    -------

    float :
        relative humidity [Pa / Pa]

    """

    e = a2e(a, T)
    return e / e_sat_func(T)


def q2rh(q, T, p, e_sat_func=e_sat_gg_water):
    """
    Calculate relative humidity from specific humidity. Source: Kraus, 'Die
    Atmosphäre der Erde', Chapter 8.1.2

    Parameters
    ----------
    q:
        specific humidity [kg / kg]
    T:
        Temperature in K
    p:
        pressure [Pa]
    e_sat_func: func, optional
        Function to estimate the saturation pressure. E.g. e_sat_gg_water for
        water and e_sat_gg_ice for ice.

    Returns
    -------

    float :
        relative humidity [Pa / Pa]

    """

    e = q2e(q, p)
    return e / e_sat_func(T)


def rh_to_iwv(relhum_lev, temp_lev, press_lev, hgt_lev,
              e_sat_func=e_sat_gg_water):
    """
    Integrate relative humidity to obtain the integrated water vapor (IWV)
    column.

    Parameters
    ----------
    relhum_lev:
        relative humidity at levels humidity [Pa / Pa]
    temp_lev:
        Temperature at levels [K]
    press_lev:
        pressure at levels [Pa]
    hgt_lev:
        altitude of levels [m]
    e_sat_func: func, optional
        Function to estimate the saturation pressure. E.g. e_sat_gg_water for
        water and e_sat_gg_ice for ice.

    Returns
    -------

    float :
        IWV [kg / m^2]

    """

    dz = np.diff(hgt_lev, axis=-1)
    relhum = (relhum_lev[..., 0:-1] + relhum_lev[..., 1:]) / 2.
    temp = (temp_lev[..., 0:-1] + temp_lev[..., 1:]) / 2.

    xp = -1.*np.log(press_lev[..., 1:] / press_lev[..., 0:-1]) / dz
    press = -1.*press_lev[..., 0:-1] / xp*(np.exp(-xp*dz)-1.) / dz

    import meteo_si.density

    q = rh2q(relhum, temp, press, e_sat_func=e_sat_func)
    rho_moist = meteo_si.density.moist_rho_q(press, temp, q)

    return np.sum(q*rho_moist*dz)
