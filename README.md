|<img src="https://raw.githubusercontent.com/euroargodev/argodmqc_owc/master/docs/_static/logo_pyowc.png" alt="pyowc logo" width="100"/>|``pyowc`` is a python library for OWC salinity calibration in Python <br><br> [![Status](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable) ![Python](https://img.shields.io/badge/python-3.6%20%7C%203.7%20%7C%203.8-blue.svg) [![Gitter](https://badges.gitter.im/Argo-floats/owc-python.svg)](https://gitter.im/Argo-floats/owc-python?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge) <br> [![codecov](https://codecov.io/gh/euroargodev/argodmqc_owc/branch/refactor-configuration/graph/badge.svg)](https://codecov.io/gh/euroargodev/argodmqc_owc) [![build](https://github.com/euroargodev/argodmqc_owc/actions/workflows/main.yml/badge.svg)](https://github.com/euroargodev/argodmqc_owc/actions/workflows/main.yml)|
|:---------:|:---------:|
This software is a python implementation of the "OWC" salinity calibration method used in Argo floats Delayed Mode Quality Control.

[Post an issue to get involved if you're interested](https://github.com/euroargodev/argodmqc_owc/issues/new/choose).

# Installation

If you are using Linux, Windows or macOS, you can simply install this package using `pip` in a virtual environment.
Assuming your virtual environment is activated:

```bash
pip install git+https://github.com/euroargodev/argodmqc_owc
```

# Software usage

A short tutorial is available on the [argopy documentation here](https://argopy.readthedocs.io/en/latest/data_quality_control.html#running-the-calibration).

For Python beginners, you can run the pyowc in this way:

In start_with_pycharm.py code, you can specify the WMO float number that you want to do analysis.
You can also add more float numbers, then the calculations of all floats will be done at the
same time.

```python
import pyowc as owc

warnings.filterwarnings("ignore", category=RuntimeWarning)

if __name__ == '__main__':

    FLOAT_NAMES = ["3901960"]  # add float names here e.g. ["3901960","3901961","3901962"]
    USER_CONFIG = owc.configuration.load()  # fetch the default configuration and parameters
    print(owc.configuration.print_cfg(USER_CONFIG))
```

## Parameters for your analysis

Parameters for the analysis are set in a configuration.py python code. 
The configuration has the same parameters as the Matlab software (https://github.com/ArgoDMQC/matlab_owc).

- You can change the default directories to locations of your historical data.
```python
        #    Climatology Data Input Paths
        'HISTORICAL_DIRECTORY': "data/climatology/"
        'HISTORICAL_CTD_PREFIX': "/historical_ctd/ctd_"
        'HISTORICAL_BOTTLE_PREFIX': "/historical_bot/bot_"
        'HISTORICAL_ARGO_PREFIX': "/historical_argo/argo_"
```
- To run the analysis,you need to have the float source file in .mat format. 
```python
        #    Float Input Path
        'FLOAT_SOURCE_DIRECTORY': "data/float_source/"
        'FLOAT_SOURCE_POSTFIX': ".mat"
```
- The output from the analysis will be saved in default directory of the code.You can change 
the default directories to locations of your constants.
```python
        #    Constants File Path
        'CONFIG_DIRECTORY': "data/constants/"
        'CONFIG_COASTLINES': "coastdat.mat"
        'CONFIG_WMO_BOXES': "wmo_boxes.mat"
        'CONFIG_SAF': "TypicalProfileAroundSAF.mat"
```
- Final step is to set your objective mapping parameters, e.g.
```python
        'MAP_USE_PV': 0
        'MAP_USE_SAF': 0

        'MAPSCALE_LONGITUDE_LARGE': 8
        'MAPSCALE_LONGITUDE_SMALL': 4
        'MAPSCALE_LATITUDE_LARGE': 4
        'MAPSCALE_LATITUDE_SMALL': 2
 ```
- Additionally, you can set a specific ranges of theta bounds for salinity anomaly plot.
The code will crete two separate plots with set ranges.
```python 
     #    Plotting Parameters
        # Theta bounds for salinity anomaly plot
        'THETA_BOUNDS': [[0, 5], [5, 20]]
```

## Plots
The plots are automatically generated. It is worth to note that only one plot will be 
displayed at one time in the PyCharm. The next plot will be displayed after closing
the window of the current plot. 

The number of generated plots at specific theta levels (from 1 to 10 theta levels) can be
currently changed in the dashboard.py code. The default is set to 2. The plots will be 
generated separately for each theta level.

```python
def plot_diagnostics(float_dir, float_name, config, levels=2):
```

# Building the documentation

If you wish to build the documentation locally, you will need a virtual environment.
Assuming your virtual environment is activated, follow these steps:

1. Install the required documentation packages
    ```bash
    pip install -r requirements-docs.txt
    ```
2. Change directory to the `docs` directory, for example on Linux:
    ```bash
    cd docs
    ```
3. Run the `sphinx-build` command:
    ```bash
    sphinx-build -M html source build -W
    ```

This will build the HTML documentation under the `docs/build/html` directory and can be viewed
using your normal web browser.

```{admonition} Note
:class: note

If you make modifications to the code or documentation configuration, you may need to delete
the `docs/source/generated` directory for the documentation to build correctly.
```

# Software history

- Major refactoring of the software for performance optimisation and to fully embrace the Pythonic way of doing this !

- [UAT: Phase 1, 2, 3](https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool) 

- Migration of the code from from BODC/NOC git to euroargodev/argodmqc_owc. See https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues/10 for more details on the migration. Contribution from [G. Maze](https://github.com/gmaze)

- Alpha experts user testings with [feedbacks available here](https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues). Contributions from: [K. Walicka](https://github.com/kamwal), [C. Cabanes](https://github.com/cabanesc), [A. Wong](https://github.com/apswong)

- BODC created [the first version of the code](https://git.noc.ac.uk/bodc/owc-software-python), following the [Matlab implementation](https://github.com/ArgoDMQC/matlab_owc).
  Contributions from: [M. Donnelly](https://github.com/matdon17), [E. Small](https://github.com/edsmall-bodc),
   [K. Walicka](https://github.com/kamwal), [A. Hale](https://github.com/halebodc), [T. Gardner](https://github.com/thogar-computer).


## New positioning of functions 
Note that functions name are not changed !

- **pyowc/core**
  - **stats.py**: brk_pt_fit, build_cov, covarxy_pv, covar_xyt_pv, noise_variance, signal_variance, fit_cond, nlbpfun
  - **finders.py**: find_10thetas, find_25boxes, find_besthit, find_ellipse, nearest_neighbour

- **pyowc/data**
  - **fetchers.py**: get_region_data, get_region_hist_locations, get_data, get_topo_grid, frontal_constraint_saf
  - **wrangling.py**: interp_climatology, map_data_grid 

- **pyowc/plot**
  - **dashboard.py**: plot_diagnostics
  - **plots.py**: cal_sal_curve_plot, sal_var_plot, t_s_profile_plot, theta_sal_plot, trajectory_plot
  - **utils.py**: create_dataframe

- **pyowc/calibration.py**: update_salinity_mapping, calc_piecewisefit

- **pyowc/configuration.py**: load_configuration, set_calseries, print_cfg

- **pyowc/tests**  # Contain all the unit tests !

- **pyowc/utilities.py**: change_dates, cal2dec, potential_vorticity, wrap_longitudes, sorter, spatial_correlation
