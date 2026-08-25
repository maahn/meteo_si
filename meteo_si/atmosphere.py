# -*- coding: utf-8 -*-
# (c) Ralph Carmichael, Public Domain Aeronautical Software (original
#     algorithm, http://www.pdas.com/programs/atmos.f90)
# (c) Maximilian Maahn 2011 (Python translation)

'''
Functions to compute properties of the 1976 US Standard Atmosphere.

'''

from collections.abc import Iterable

import numpy as np

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

    # index of the layer containing h; altitudes below the first and above
    # the last table entry are extrapolated from the outermost layers
    i = int(np.clip(np.searchsorted(htab, h, side='right') - 1,
                    0, len(htab) - 2))

    tgrad = gtab[i]
    tbase = ttab[i]
    deltah = h - htab[i]
    tlocal = tbase + tgrad * deltah
    theta = tlocal / ttab[0]  # temperature ratio

    if tgrad == 0.0:  # pressure ratio
        delta = ptab[i] * np.exp(-GMR * deltah / tbase)
    else:
        delta = ptab[i] * (tbase / tlocal) ** (GMR / tgrad)

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

    if isinstance(height, Iterable):
        height = np.asarray(height, dtype=float)
        assert height.ndim == 1, (
            "height must be a scalar or a one-dimensional array")
        # the reference algorithm processes only a single value at a time
        ratios = np.array([_Atmosphere(hh / 1000.0) for hh in height])
        sigma, delta, theta = ratios.reshape(-1, 3).T
    else:
        sigma, delta, theta = _Atmosphere(height / 1000.0)

    # the ratios are relative to sea-level standard conditions
    density = sigma * 1.2250
    pressure = delta * 101325
    temperature = theta * 288.15

    return density, pressure, temperature
