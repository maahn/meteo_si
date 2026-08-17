import numpy.testing as npt
import pytest

import meteo_si


def test_moist_rho_q_reduces_to_ideal_gas_law_when_dry():
    p = 90000.0
    T = 280.0
    npt.assert_allclose(
        meteo_si.density.moist_rho_q(p, T, q=0.0),
        p / (meteo_si.constants.Rair * T),
    )


def test_moist_air_is_less_dense_than_dry_air():
    # Water vapor is lighter than dry air, so adding humidity at fixed
    # pressure and temperature lowers the density.
    p = 90000.0
    T = 280.0
    rho_dry = meteo_si.density.moist_rho_q(p, T, q=0.0)
    rho_moist = meteo_si.density.moist_rho_q(p, T, q=0.01)
    assert rho_moist < rho_dry


def test_moist_rho_rh_matches_moist_rho_q_via_rh2q():
    p = 92000.0
    T = 285.0
    rh = 0.4
    q = meteo_si.humidity.rh2q(rh, T, p)
    npt.assert_allclose(
        meteo_si.density.moist_rho_rh(p, T, rh),
        meteo_si.density.moist_rho_q(p, T, q),
    )


def test_moist_rho_rh_rejects_percent_input():
    with pytest.raises(TypeError):
        meteo_si.density.moist_rho_rh(92000.0, 285.0, 40.0)


def test_moist_rho_q_hydrometeor_loading_increases_density_term():
    # qm represents the mass mixing ratio of condensed hydrometeors that
    # still contribute to the air mass; increasing it increases the
    # computed density for otherwise fixed p, T, q.
    p = 90000.0
    T = 280.0
    q = 0.005
    rho_no_hydrometeors = meteo_si.density.moist_rho_q(p, T, q, qm=0.0)
    rho_with_hydrometeors = meteo_si.density.moist_rho_q(p, T, q, qm=0.01)
    assert rho_with_hydrometeors > rho_no_hydrometeors


def test_moist_rho_q_raises_on_strongly_negative_density():
    with pytest.raises(ValueError):
        meteo_si.density.moist_rho_q(90000.0, 280.0, q=0.0, qm=10.0)
