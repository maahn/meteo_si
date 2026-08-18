# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`meteo_si` is a small Python library of atmospheric-science formulas (humidity, density, temperature, wind) implemented consistently in SI units unless otherwise noted in a function's docstring. Supports Python 3.10–3.14.

## Commands

Install dependencies and the package itself:
```bash
pip install -e .[test,lint]
```

Run the test suite (from repo root):
```bash
make test
# equivalent to:
py.test --pyargs meteo_si --cov-report term-missing --cov=meteo_si
```

Run a single test:
```bash
py.test meteo_si/tests/test_humidity.py::test_e_sat_gg_water_reference_point
```

Lint:
```bash
make flake8
# equivalent to:
flake8 --extend-ignore N802,N806 $(find . -name '*.py' | grep -v /doc/)
```
`N802`/`N806` (mixed-case naming) are intentionally ignored — variable/function names mirroring physics notation (e.g. `T`, `Rair`, `Mwml`) are fine. Use `--extend-ignore`, not `--ignore`: flake8's own default ignore list already excludes `E226`/`W504` (arithmetic-operator spacing), and `--ignore` replaces that default list instead of adding to it — that mistake previously made `make flake8` fail on pre-existing style in every formula-heavy module.

Build docs (Sphinx, config in `doc/conf.py`):
```bash
cd doc && make html
```

CI (`.github/workflows/tests.yml`) runs the test matrix on Python 3.10–3.14 across Linux/Windows/macOS, plus a separate flake8 job. `.github/workflows/publish.yml` builds and publishes to PyPI via trusted publishing (OIDC) on GitHub release.

## Architecture

- `meteo_si/__init__.py` imports each topic module as a **namespace**, not a flat function list:
  ```python
  from . import temperature, humidity, density, wind, constants, atmosphere
  ```
  So calling code uses `meteo_si.humidity.rh2q(...)`, `meteo_si.density.moist_rho_q(...)`, etc. — not `meteo_si.rh2q(...)`. Keep new functions inside the appropriate topic module and export them via that module's `__all__`.

- **Module responsibilities and dependency direction**: `constants.py` (no internal deps) → `temperature.py` / `humidity.py` (both depend on `constants`, and `temperature` also depends on `humidity` for `T_virt_rh`) → `density.py` (depends on `constants` and `humidity`). `wind.py` and `atmosphere.py` are both independent (circular-mean helpers, and the 1976 US Standard Atmosphere model, respectively — neither needs anything from the other modules). Respect this layering when adding functions — e.g. don't import `density` from `humidity`.

- `humidity.py` has two saturation-vapor-pressure-over-water formulas that are **not interchangeable**: `e_sat_gg_water` (WMO CIMO Guide 2008, Magnus-type) and `e_sat_goffgratch_water` (Goff & Gratch 1946, more accurate over a wider temperature range but pricier to evaluate). Despite the name, `e_sat_gg_water` is the CIMO one, not Goff-Gratch — that's legacy naming, kept for backward compatibility. Every humidity function that depends on saturation pressure (`rh2q`, `rh2a`, `a2rh`, `q2rh`, `rh_to_iwv`) takes an `e_sat_func` parameter (defaulting to `e_sat_gg_water`) for exactly this reason — pick explicitly rather than assuming the default is "the accurate one" or "the Goff-Gratch one".

- Physical quantities are passed as bare floats/arrays with units documented per-function in the docstring (numpy-style `Parameters`/`Returns` sections), not via unit-aware types. When adding a function, follow the existing docstring convention exactly (units in brackets, e.g. `[Pa]`, `[kg/kg]`, `[kg/m3]`) since that's the only unit documentation that exists.

- Relative-humidity functions (`rh2q`, `rh2a`, `T_virt_rh`, `moist_rho_rh`) accept `rh` as a fraction (Pa/Pa), not a percentage, and raise `TypeError` if any value exceeds 5 — this is a sanity guard against callers passing percentages by mistake, not a hard physical bound.

- Tests live one file per topic module (`meteo_si/tests/test_humidity.py`, `test_density.py`, `test_temperature.py`, `test_wind.py`), mirroring the module layout — add new tests to the matching file rather than a catch-all test module.

- `meteo_si/due.py` is an unmodified `duecredit` stub (for citation support) and isn't wired into any module — leave it as-is unless adding citation annotations.

- Version metadata lives in `meteo_si/version.py` (`__version__`, built from `_version_major`/`_version_minor`/`_version_micro`/`_version_extra`); `pyproject.toml` reads it dynamically via `tool.setuptools.dynamic`. All other package metadata (name, description, classifiers, dependencies) lives directly in `pyproject.toml`.

- `doc/sphinxext/` has two live custom Sphinx extensions: `math_dollar.py` (lets docstrings use `$...$` for inline math) and `github.py` (adds `:ghissue:`/`:ghpull:`/`:ghuser:`/`:ghcommit:` RST roles). `doc/conf.py` must use `autodoc_default_options` (a dict), not the removed `autodoc_default_flags`, or `.. automodule::` directives in `doc/index.rst` silently render with no members.
