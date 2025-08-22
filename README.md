|<img src="https://raw.githubusercontent.com/euroargodev/argodmqc_owc/master/docs/_static/pyowc-logo.png" alt="pyowc logo" width="100"/>|``pyowc`` is a python library for OWC salinity calibration in Python <br><br> [![Status](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable) ![Python](https://img.shields.io/badge/python-3.6%20%7C%203.7%20%7C%203.8-blue.svg) [![Gitter](https://badges.gitter.im/Argo-floats/owc-python.svg)](https://gitter.im/Argo-floats/owc-python?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge) <br> [![codecov](https://codecov.io/gh/euroargodev/argodmqc_owc/branch/refactor-configuration/graph/badge.svg)](https://codecov.io/gh/euroargodev/argodmqc_owc) [![build](https://github.com/euroargodev/argodmqc_owc/actions/workflows/main.yml/badge.svg)](https://github.com/euroargodev/argodmqc_owc/actions/workflows/main.yml)|
|:---------:|:---------:|

This software is a python implementation of the "OWC" salinity calibration method used in Argo floats Delayed Mode Quality Control.

[Post an issue to get involved if you're interested](https://github.com/euroargodev/argodmqc_owc/issues/new/choose).

# General Guidance

To use this software, you'll need Python, and ideally a virtual environment with the package installed.  
A virtual environment is not absolutely essential — you can install the package globally — but it is recommended to avoid issues.

There are two ways of working with this software:

## 1. Installation and General Usage via PyPI

This method is intended for general usage, without modifying the codebase.  
If you intend to use the software this way, follow the documentation from the **General Usage** section.

## 2. Installation and Development Work via GitHub

This method is intended for development work and gives access to the codebase, it is intended for those wanting to develop the code. If this is what you intend to do, follow the documentation from the **Developer Usage** section.

---

## Overview

To use the app installed via PyPi there are 3 steps to follow.

1. Pip install it
  Run `pip install argodmqc-owc`


2. Setup folder structure
  The app requires a specific folder structure, and these folders & files are referenced in the config JSON file.
  The example structure here is reflected in the example config JSON file.



## General Usage

To use the app installed via PyPi there are 3 steps to follow.

1. Pip install it

  Run `pip install argodmqc-owc`


2. Setup folder structure

  The app requires a specific folder structure, and these folders & files are referenced in the config JSON file.
  The example structure here is reflected in the example config JSON file. Please note that the `data` folder is available [here](https://github.com/euroargodev/argodmqc_owc/tree/feature/pypi-preparation/data)
  - data
    - climatology
      - historical_argo
      - historical_bot
      - historical_ctd
    - constants
      - bathymetry
      - coastline
      - reefs
    - float_calib
    - float_mapped
    - float_plots
    - float_source

3. Create the config JSON file

  See the example config JSON file [here](https://github.com/euroargodev/argodmqc_owc/blob/feature/pypi-preparation/owc_config.json)


4. Run the software

  Use the example script provided here [here](https://github.com/euroargodev/argodmqc_owc/blob/feature/pypi-preparation/start.py) to run the DMQC software.

---

## Developer Usage

TODO

---


3. Create the config JSON file
  See the example config JSON file here: 



4. Create the script to import and then use the DMQC software
  See the example python file here:

---

## Developer Usage

TODO

---

## Virtual Environments


To create a virtual environment:
- **Mac/Linux**

  `python3 -m venv .venv`

  `source .venv/bin/activate`

- **Windows**

  `python -m venv .venv`

  `.\.venv\Scripts\Activate`



## Installation and usage with Poetry

- **Step 1**  
  Make sure you have Python installed, along with `virtualenv`.

- **Step 2**  
  Clone the repository and open it in your code editor.

- **Step 3**  
  Create a new virtual environment.

- **Step 4**  
  Install Poetry: `pip install poetry`

- **Step 5**
  Install the dependencies: `poetry install --no-root.

If any messages appear with 'poetry not found', try prefixing your command with `python` or `python -m`



## Running the Linter with Poetry

`poetry install --no-root --with lint`

`poetry run ruff check`


## Running the Docs Builder with Poetry

`poetry install --with docs`

`cd docs`

`poetry run sphinx-build -M html source build -W`

## Running the Tests with Poetry

`poetry install --with tests`

`poetry run pytest`



## Software usage

- **Running with Poetry**
  

Run the code (start.py): `poetry run run-floats`.


A short tutorial is available on the [argopy documentation here](https://argopy.readthedocs.io/en/latest/data_quality_control.html#running-the-calibration).

For Python beginners, you can run the pyowc in this way:

In start.py, you can specify the WMO float number that you want to do analysis.
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




## Software history

- Major refactoring of the software for performance optimisation and to fully embrace the Pythonic way of doing this !

- [UAT: Phase 1, 2, 3](https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool) 

- Migration of the code from from BODC/NOC git to euroargodev/argodmqc_owc. See https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues/10 for more details on the migration. Contribution from [G. Maze](https://github.com/gmaze)

- Alpha experts user testings with [feedbacks available here](https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues). Contributions from: [K. Walicka](https://github.com/kamwal), [C. Cabanes](https://github.com/cabanesc), [A. Wong](https://github.com/apswong)

- BODC created [the first version of the code](https://git.noc.ac.uk/bodc/owc-software-python), following the [Matlab implementation](https://github.com/ArgoDMQC/matlab_owc).
  Contributions from: [M. Donnelly](https://github.com/matdon17), [E. Small](https://github.com/edsmall-bodc),
   [K. Walicka](https://github.com/kamwal), [A. Hale](https://github.com/halebodc), [T. Gardner](https://github.com/thogar-computer).
