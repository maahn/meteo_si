import numpy as np
import numpy.testing as npt
import pytest

import meteo_si


def test_kelvin_celsius_roundtrip():
    T = np.array([200.0, 250.0, 273.15, 300.0, 320.0])
    C = meteo_si.temperature.kelvin_2_celsius(T)
    npt.assert_allclose(meteo_si.temperature.celsius_to_kelvin(C), T)


def test_kelvin_2_celsius_reference_point():
    # 0 degC is exactly 273.15 K by definition.
    npt.assert_allclose(meteo_si.temperature.kelvin_2_celsius(273.15), 0.0)


def test_T_virt_q_equals_T_when_dry():
    T = 290.0
    npt.assert_allclose(meteo_si.temperature.T_virt_q(T, 0.0), T)


def test_T_virt_q_exceeds_T_for_moist_air():
    # Water vapor is lighter than dry air, so virtual temperature of moist
    # air is always >= the actual temperature.
    T = 290.0
    T_virt = meteo_si.temperature.T_virt_q(T, 0.01)
    assert T_virt > T


def test_T_virt_rh_matches_T_virt_q_via_rh2q():
    T = 285.0
    rh = 0.6
    p = 95000.0
    q = meteo_si.humidity.rh2q(rh, T, p)
    npt.assert_allclose(
        meteo_si.temperature.T_virt_rh(T, rh, p),
        meteo_si.temperature.T_virt_q(T, q),
    )


def test_T_virt_rh_rejects_percent_input():
    # rh must be given as a fraction (Pa/Pa), not a percentage.
    with pytest.raises(TypeError):
        meteo_si.temperature.T_virt_rh(285.0, 60.0, 95000.0)
