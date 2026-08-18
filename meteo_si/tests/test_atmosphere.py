import numpy as np
import numpy.testing as npt

import meteo_si


def test_usStandard_sea_level():
    density, pressure, temperature = meteo_si.atmosphere.usStandard(0.0)
    npt.assert_allclose(density, 1.225, rtol=1e-6)
    npt.assert_allclose(pressure, 101325.0, rtol=1e-6)
    npt.assert_allclose(temperature, 288.15, rtol=1e-6)


def test_usStandard_decreases_with_height():
    heights = np.array([0.0, 5000.0, 11000.0, 20000.0])
    density, pressure, temperature = meteo_si.atmosphere.usStandard(heights)
    assert np.all(np.diff(density) < 0)
    assert np.all(np.diff(pressure) < 0)
    # Temperature drops through the troposphere, then holds ~constant
    # across the tropopause (11-20 km) in the 1976 standard atmosphere.
    assert temperature[1] < temperature[0]
    npt.assert_allclose(temperature[2], temperature[3], rtol=1e-3)


def test_usStandard_scalar_matches_array_input():
    heights = np.array([0.0, 11000.0, 20000.0])
    scalar_results = [meteo_si.atmosphere.usStandard(h) for h in heights]
    array_results = meteo_si.atmosphere.usStandard(heights)
    for i, (density, pressure, temperature) in enumerate(scalar_results):
        npt.assert_allclose(array_results[0][i], density)
        npt.assert_allclose(array_results[1][i], pressure)
        npt.assert_allclose(array_results[2][i], temperature)
