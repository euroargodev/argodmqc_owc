# pyowc: OWC salinity calibration in Python

![build](https://github.com/euroargodev/argodmqc_owc/workflows/build/badge.svg)
[![codecov](https://codecov.io/gh/euroargodev/argodmqc_owc/branch/refactor-configuration/graph/badge.svg)](https://codecov.io/gh/euroargodev/argodmqc_owc)
[![Requirements Status](https://requires.io/github/euroargodev/argodmqc_owc/requirements.svg?branch=master)](https://requires.io/github/euroargodev/argodmqc_owc/requirements/?branch=refactor-configuration)
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)

[![Gitter](https://badges.gitter.im/Argo-floats/owc-python.svg)](https://gitter.im/Argo-floats/owc-python?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

This software is a python implementation of the "OWC" salinity calibration method used in Argo floats Delayed Mode Quality Control.

This software is in very active development and its usage may change any time. Please [post an issue to get involved if you're interested](https://github.com/euroargodev/argodmqc_owc/issues/new/choose).

# Installation

```bash
pip install pyowc
```

# Software usage

[![badge](https://img.shields.io/static/v1.svg?logo=Jupyter&label=Pangeo+Binder&message=Click+here+to+try+this+software+online+!&color=blue&style=for-the-badge)](https://binder.pangeo.io/v2/gh/euroargodev/argodmqc_owc/refactor-configuration?urlpath=lab/tree/docs/Tryit.ipynb)

Otherwise, in a python script, you can calibrate salinity data this way:

```python
import pyowc as owc

FLOAT_NAME = "3901960"
USER_CONFIG = owc.configuration.load() # fetch the default configuration and parameters

owc.configuration.set_calseries("/", FLOAT_NAME, USER_CONFIG)
owc.calibration.update_salinity_mapping("/", FLOAT_NAME, USER_CONFIG)
owc.calibration.calc_piecewisefit("/", FLOAT_NAME, USER_CONFIG)
```

### Parameters for your analysis

Parameters for the analysis are set in a configuration dictionary. 
The configuration has the same parameters as the Matlab software.

You can change parameters this way:
```python
USER_CONFIG['MAP_USE_PV'] = 0 # Possibly tune options
```

And you can see it content like this:
```python
print(owc.configuration.print_cfg(USER_CONFIG))
```

### Plots

```python
owc.plot.dashboard("/", FLOAT_NAME, USER_CONFIG)
```

# Software history

- Major refactoring of the software for performance optimisation and to fully embrace the Pythonic way of doing this !

- [UAT: Phase 1, 2, 3](https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool) 

- Migration of the code from from BODC/NOC git to euroargodev/argodmqc_owc. See https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues/10 for more details on the migration. Contribution from [G. Maze](https://github.com/gmaze)

- Alpha experts user testings with [feedbacks available here](https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues). Contributions from: [K. Walicka](https://github.com/kamwal), [C. Cabanes](https://github.com/cabanesc), [A. Wong](https://github.com/apswong)

- BODC created [the first version of the code](https://git.noc.ac.uk/bodc/owc-software-python), following the [Matlab implementation](https://github.com/ArgoDMQC/matlab_owc). Contributions from: [M. Donnelly](https://github.com/matdon17), [E. Small](https://github.com/edsmall-bodc), [K. Walicka](https://github.com/kamwal).


## New positioning of functions 
Note that functions name are not changed !

- **pyowc/core**
  - **stats.py**: brk_pt_fit, build_cov, covarxy_pv, covar_xyt_pv, noise_variance, signal_variance, fit_cond, nlbpfun
  - **finders.py**: find_10thetas, find_25boxes, find_besthit, find_ellipse, nearest_neighbour

- **pyowc/data**
  - **fetchers.py**: get_region_data, get_region_hist_locations, get_data, get_topo_grid
  - **wrangling.py**: interp_climatology, map_data_grid 

- **pyowc/plot**
  - **dashboard.py**: plot_diagnostics
  - **plots.py**: cal_sal_curve_plot, sal_var_plot, t_s_profile_plot, theta_sal_plot, trajectory_plot
  - **utils.py**: create_dataframe

- **pyowc/calibration.py**: update_salinity_mapping, calc_piecewisefit

- **pyowc/configuration.py**: load_configuration, set_calseries, print_cfg

- **pyowc/tests**  # Contain all the unit tests !

- **pyowc/utilities.py**: change_dates, cal2dec, potential_vorticity, wrap_longitudes, sorter, spatial_correlation