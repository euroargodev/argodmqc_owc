|<img src="https://raw.githubusercontent.com/euroargodev/argodmqc_owc/master/docs/_static/pyowc-logo.png" alt="pyowc logo" width="100"/>|``pyowc`` is a python library for OWC salinity calibration in Python <br><br> [![Status](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable) [![Pypi][pip-badge]][pip-link] [![Chat][chat-badge]][chat-link] <br> ![codecov][cov-badge] ![CI][ci-badge]|
|:---------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|

[ci-badge]: https://github.com/euroargodev/argodmqc_owc/actions/workflows/main.yml/badge.svg
[pip-badge]: https://img.shields.io/pypi/v/argopy
[pip-link]: https://pypi.org/project/argopy/
[chat-badge]: https://badges.gitter.im/Argo-floats/owc-python.svg
[chat-link]: https://gitter.im/Argo-floats/owc-python?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
[cov-badge]: https://codecov.io/gh/euroargodev/argodmqc_owc/branch/master/graph/badge.svg
[cov-link]: https://codecov.io/gh/euroargodev/argodmqc_owc

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
  The JSON file includes all data directory and float configurations and settings. 
  
  See the example config JSON file [here](https://github.com/euroargodev/argodmqc_owc/blob/feature/pypi-preparation/owc_config.json)


4. Run the software

  Use the example script provided here [here](https://github.com/euroargodev/argodmqc_owc/blob/feature/pypi-preparation/start_owc_python.py) to run the DMQC software.
  
  Please note that the config JSON is checked before any processing runs, so any errors in the config will be reported back.

  The script contains a line which has the floats to be processed: `FLOAT_NAMES = ["3901960"]`.
  This can be one or more floats, and when the code is ran they will be processed in turn.

---

## Developer Usage

Using the app as a developer requires git-cloning rather than pip-installing so direct access to the code is possible, and modifying the code is easy. To use the app as a developer it is recommended you follow the sections below in the order prescribed

1. Virtual Environments
2. Installing Poetry
3. Cloning the repository
4. Installing the dependencies
5. Running the linting & tests and docs builder
6. Executing the DMQC code

---

## Virtual Environments

A virtual environment is recommended to work in as the dependencies wont conflict with any globally installed packages.

To create a virtual environment:

- **Mac/Linux**

  `python3 -m venv .venv`

  `source .venv/bin/activate`

- **Windows**

  `python -m venv .venv`

  `.\.venv\Scripts\Activate`



## Installing Poetry

Poetry is a dependency management tool and the software uses a `pyproject.toml` file to handle the dependencies.


To Install Poetry, run: `pip install poetry`

If any messages appear with 'poetry not found', try prefixing your command with `python` or `python -m`


## Cloning the Repository

The repository can be cloned by clicking the green `<> Code` button near the top of the main page on Github. Follow the prompts to either clone it via the command line, or open with Github Desktop.

It is recommended to clone the repository to a new folder. Make sure you are in this folder with your virtual environment activated and the repository cloned before moving to the next step. You need to make sure you are at the same level as the `pyproject.toml` file.


## Installing the dependencies

Install the dependencies with: `poetry install --no-root`

This will take a few seconds, and you should see a list of the installed packages in the terminal window.


## Running the linting & tests and docs builder

The dependencies for running these utiltiies are also packaged up with Poetry, and they can be ran as follows:

### Running the Linter with Poetry

`poetry install --no-root --with lint`

`poetry run ruff check`


### Running the Docs Builder with Poetry

`poetry install --with docs`

`cd docs`

`poetry run sphinx-build -M html source build -W`

### Running the Tests with Poetry

`poetry install --with tests`

`poetry run pytest`




## Executing the DMQC code
  

Open the file `start_owc_python.py`

Look at lines 15 & 16, these are to be changed if different floats are to be processed, or a different configuration is to be used.


Run the code (start_owc_python.py): `poetry run run-floats`.


A short tutorial is available on the [argopy documentation here](https://argopy.readthedocs.io/en/latest/data_quality_control.html#running-the-calibration).

For Python beginners, you can run the pyowc in this way:

In start_owc_python.py, you can specify the WMO float number that you want to run the analysis on.
You can also add more float numbers, then the calculations of all floats will be done at the
same time.

```python

    FLOAT_NAMES = ["3901960"]  # add float names here e.g. ["3901960","3901961","3901962"]
```


## Parameters for your analysis

Parameters for the analysis are set in a owc_config.json python code. 
The configuration has the same parameters as the Matlab software (https://github.com/ArgoDMQC/matlab_owc).

- You can change the default directories to locations of your historical data.
 
      #    Climatology Data Input Paths
      'HISTORICAL_DIRECTORY': "data/climatology/"
      'HISTORICAL_CTD_PREFIX': "/historical_ctd/ctd_"
      'HISTORICAL_BOTTLE_PREFIX': "/historical_bot/bot_"
      'HISTORICAL_ARGO_PREFIX': "/historical_argo/argo_"

- To run the analysis,you need to have the float source file in .mat format. 

      #    Float Input Path
      'FLOAT_SOURCE_DIRECTORY': "data/float_source/"
      'FLOAT_SOURCE_POSTFIX': ".mat"

- The output from the analysis will be saved in default directory of the code.You can change 
the default directories to locations of your constants.

      #    Constants File Path
      'CONFIG_DIRECTORY': "data/constants/"
      'CONFIG_COASTLINES': "coastdat.mat"
      'CONFIG_WMO_BOXES': "wmo_boxes.mat"
      'CONFIG_SAF': "TypicalProfileAroundSAF.mat"

- To set your objective mapping parameters update the following, e.g.

      'MAP_USE_PV': 0
      'MAP_USE_SAF': 0

      'MAPSCALE_LONGITUDE_LARGE': 8
      'MAPSCALE_LONGITUDE_SMALL': 4
      'MAPSCALE_LATITUDE_LARGE': 4
      'MAPSCALE_LATITUDE_SMALL': 2

- To manually set the calibration values update the following entires below.
The calibration contains two additional fields `splits` and `exclusions`. The splits field refers to the splitting of profiles, and if left empty will process all the data. If the data is to be split up, then the value needs to be changed to a list of profiles to split into. 
For example if splits is [5, 10] then splits will be made at profiles 1-4, 5-10, and then 10 until the end.
If the exclusions list is empty then all data is processed, but profiles can be added to the exclusions list if they are to be omitted from the program.

      'CALIB_PROFILE_NO': [],
      'BREAKS': [],
      'MAX_BREAKS': 4,
      'SPLITS': [],
      'EXCLUSIONS': [],
      'USE_THETA_LT': [],
      'USE_THETA_GT': [],
      'USE_PRES_LT': [],
      'USE_PRES_GT': [],
      'USE_PERCENT_GT': 0.5,

- Additionally, you can set a specific ranges of theta bounds for salinity anomaly plot.
The code will crete two separate plots with set ranges 'THETA_BOUNDS': [[0, 5], [5, 20]]. To create only one plot with a specific range use the 'THETA_BOUNDS': [[0, 20]].

      #    Plotting Parameters
      'THETA_BOUNDS': [[0, 5], [5, 20]]



## Software history

Check out software history at: https://github.com/euroargodev/argodmqc_owc/network
