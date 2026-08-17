import numpy as np
import numpy.testing as npt

import meteo_si


def test_circular_mean_deg_symmetric_around_north():
    # 350 and 10 degrees straddle north/0; the mean is 0 (equivalently
    # 360) modulo a full turn.
    result = meteo_si.wind.circular_mean_deg(np.array([350.0, 10.0]))
    npt.assert_allclose(result % 360.0, 0.0, atol=1e-9)


def test_circular_mean_deg_of_close_headings():
    npt.assert_allclose(
        meteo_si.wind.circular_mean_deg(np.array([80.0, 100.0])), 90.0)


def test_circular_mean_matches_arithmetic_mean_for_close_angles():
    # For angles that don't wrap around, the circular mean should agree
    # with the ordinary arithmetic mean.
    angles = np.deg2rad(np.array([10.0, 20.0, 30.0]))
    npt.assert_allclose(
        meteo_si.wind.circular_mean(angles), np.mean(angles), atol=1e-9)


def test_circular_mean_returns_nan_with_any_nan():
    angles = np.array([0.1, np.nan, 0.3])
    assert np.isnan(meteo_si.wind.circular_mean(angles))


def test_nan_circular_mean_ignores_nans():
    with_nan = np.array([0.1, np.nan, 0.3, 0.5])
    without_nan = np.array([0.1, 0.3, 0.5])
    npt.assert_allclose(
        meteo_si.wind.nan_circular_mean(with_nan),
        meteo_si.wind.nan_circular_mean(without_nan),
    )


def test_circular_mean_deg_returns_value_in_0_360():
    result = meteo_si.wind.circular_mean_deg(np.array([190.0, 200.0]))
    assert 0.0 <= result < 360.0 + 1e-9
