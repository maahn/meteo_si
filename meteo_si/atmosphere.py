# -*- coding: utf-8 -*-
# (c) Ralph Carmichael, Public Domain Aeronautical Software (original
#     algorithm, http://www.pdas.com/programs/atmos.f90)
# (c) Maximilian Maahn 2011 (Python translation)

from collections.abc import Iterable

import numpy as np

'''
Functions to compute properties of the 1976 US Standard Atmosphere.

'''

__all__ = ["usStandard"]


def _Atmosphere(alt):
    """
    Compute the properties of the 1976 standard atmosphere to 86 km.

    If ``alt`` > 86, the values returned will not be correct, but they
    will not be too far removed from the correct values for density. The
    reference document does not use the terms pressure and temperature
    above 86 km.

    Parameters
    ----------
    alt:
        geometric altitude [km]

    Returns
    -------

    sigma : float
        density ratio (density / sea-level standard density)
    delta : float
        pressure ratio (pressure / sea-level standard pressure)
    theta : float
        temperature ratio (temperature / sea-level standard temperature)

    """

    REARTH = 6369.0  # radius of the Earth (km)
    GMR = 34.163195  # hydrostatic constant
    NTAB = 8  # number of entries in the defining tables

    htab = np.array([0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852])
    ttab = np.array(
        [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946])
    ptab = np.array(
        [
            1.0,
            2.233611e-1,
            5.403295e-2,
            8.5666784e-3,
            1.0945601e-3,
            6.6063531e-4,
            3.9046834e-5,
            3.68501e-6,
        ]
    )
    gtab = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0])

    # convert geometric to geopotential altitude
    h = alt * REARTH / (alt + REARTH)

    i = 1
    j = NTAB  # setting up for binary search
    while True:
        k = (i + j) // 2  # integer division
        if h < htab[k - 1]:
            j = k
        else:
            i = k
        if j <= i + 1:
            break

    tgrad = gtab[i - 1]  # i will be in 1...NTAB-1
    tbase = ttab[i - 1]
    deltah = h - htab[i - 1]
    tlocal = tbase + tgrad * deltah
    theta = tlocal / ttab[0]  # temperature ratio

    if tgrad == 0.0:  # pressure ratio
        delta = ptab[i-1] * np.exp(-GMR * deltah / tbase)
    else:
        delta = ptab[i-1] * (tbase / tlocal) ** (GMR / tgrad)

    sigma = delta / theta  # density ratio

    return sigma, delta, theta


def usStandard(height):
    """
    Compute the US standard atmosphere (1976).

    Source: http://www.pdas.com/programs/atmos.f90

    Parameters
    ----------
    height:
        height [m], as a single value or an array

    Returns
    -------

    density : float
        air density [kg / m3]
    pressure : float
        air pressure [Pa]
    temperature : float
        air temperature [K]

    """

    # check whether height is float or array:
    if isinstance(height, Iterable):
        density = np.ones_like(height, dtype=float)
        pressure = np.ones_like(height, dtype=float)
        temperature = np.ones_like(height, dtype=float)
        for ii, hh in enumerate(height):
            # the reference algorithm processes only a single value at a
            # time, so make sure it is not itself a vector
            assert not isinstance(hh, Iterable)
            density[ii], pressure[ii], temperature[ii] = (
                _Atmosphere(hh / 1000.0)
            )
    else:
        density, pressure, temperature = _Atmosphere(height / 1000.0)

    # results are normed to standard conditions:
    density = density * 1.2250
    pressure = pressure * 101325
    temperature = temperature * 288.15

    return density, pressure, temperature
