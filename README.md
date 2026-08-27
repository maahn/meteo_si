## meteo_si

[![Tests](https://github.com/maahn/meteo_si/actions/workflows/tests.yml/badge.svg)](https://github.com/maahn/meteo_si/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/meteo-si/badge/?version=latest)](https://meteo-si.readthedocs.io/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/meteo_si.svg)](https://pypi.org/project/meteo_si/)
[![Downloads](https://static.pepy.tech/badge/meteo_si)](https://pepy.tech/project/meteo_si)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/meteo_si.svg)](https://anaconda.org/conda-forge/meteo_si)

Meteo SI is a collection of functions for atmospheric sciences following the
SI-convention (unless stated otherwise). See
https://meteo-si.readthedocs.io/en/latest/ for documentation.

### Installation

```
pip install meteo_si
```

or, via conda-forge:

```
conda install conda-forge::meteo_si
```

### Development

```
pip install -e .[test,lint]
make test     # run the test suite
make flake8   # run the linter
```

Supported Python versions: 3.10–3.14.
