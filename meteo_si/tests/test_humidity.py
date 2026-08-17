import numpy as np
import numpy.testing as npt
import pytest

import meteo_si


def test_a2e_e2a_roundtrip():
    a = np.array([0.001, 0.005, 0.01, 0.02])
    T = np.array([250.0, 270.0, 290.0, 310.0])
    e = meteo_si.humidity.a2e(a, T)
    npt.assert_allclose(meteo_si.humidity.e2a(e, T), a)


def test_e2q_q2e_roundtrip():
    e = np.array([500.0, 1000.0, 2000.0])
    p = 90000.0
    q = meteo_si.humidity.e2q(e, p)
    npt.assert_allclose(meteo_si.humidity.q2e(q, p), e)


def test_e_sat_gg_water_reference_point():
    # Saturation vapor pressure over water at 0 degC is ~611.2 Pa, the
    # standard reference value used to anchor the CIMO Guide formula.
    npt.assert_allclose(
        meteo_si.humidity.e_sat_gg_water(273.15), 611.2, rtol=1e-6)


def test_e_sat_gg_ice_reference_point():
    npt.assert_allclose(
        meteo_si.humidity.e_sat_gg_ice(273.15), 611.2, rtol=1e-6)


def test_e_sat_increases_with_temperature():
    T = np.array([250.0, 270.0, 290.0, 310.0])
    e_sat = meteo_si.humidity.e_sat_gg_water(T)
    assert np.all(np.diff(e_sat) > 0)


def test_rh2q_q2rh_roundtrip():
    rh = 0.5
    T = 280.0
    p = 95000.0
    q = meteo_si.humidity.rh2q(rh, T, p)
    npt.assert_allclose(meteo_si.humidity.q2rh(q, T, p), rh)


def test_rh2a_a2rh_roundtrip():
    rh = 0.7
    T = 295.0
    a = meteo_si.humidity.rh2a(rh, T)
    npt.assert_allclose(meteo_si.humidity.a2rh(a, T), rh)


def test_rh2q_rejects_percent_input():
    with pytest.raises(TypeError):
        meteo_si.humidity.rh2q(70.0, 280.0, 95000.0)


def test_rh2a_rejects_percent_input():
    with pytest.raises(TypeError):
        meteo_si.humidity.rh2a(70.0, 280.0)


def test_rh2q_with_e_sat_gg_ice():
    rh = 0.5
    T = 260.0
    p = 90000.0
    q = meteo_si.humidity.rh2q(
        rh, T, p, e_sat_func=meteo_si.humidity.e_sat_gg_ice)
    # Saturation pressure over ice is lower than over water at the same
    # sub-freezing temperature, so the resulting specific humidity is too.
    q_water = meteo_si.humidity.rh2q(
        rh, T, p, e_sat_func=meteo_si.humidity.e_sat_gg_water)
    assert q < q_water


def test_rh_to_iwv_constant_profile():
    # A well-mixed, constant relative-humidity profile should give a
    # strictly positive IWV, and doubling the layer thickness (height
    # spacing) should roughly double it.
    n = 10
    relhum_lev = np.full(n, 0.5)
    temp_lev = np.full(n, 280.0)
    press_lev = np.linspace(95000.0, 80000.0, n)
    hgt_lev = np.linspace(0.0, 2000.0, n)

    iwv = meteo_si.humidity.rh_to_iwv(
        relhum_lev, temp_lev, press_lev, hgt_lev)
    assert iwv > 0

    iwv_double_height = meteo_si.humidity.rh_to_iwv(
        relhum_lev, temp_lev, press_lev, hgt_lev * 2)
    npt.assert_allclose(iwv_double_height, iwv * 2, rtol=1e-6)
