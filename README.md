# py_owc: OWC salinity calibration in Python

![build](https://github.com/euroargodev/argodmqc_owc/workflows/build/badge.svg)
[![codecov](https://codecov.io/gh/euroargodev/argodmqc_owc/branch/master/graph/badge.svg)](https://codecov.io/gh/euroargodev/argodmqc_owc)
[![Requirements Status](https://requires.io/github/euroargodev/argodmqc_owc/requirements.svg?branch=master)](https://requires.io/github/euroargodev/argodmqc_owc/requirements/?branch=master)
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)

[![Gitter](https://badges.gitter.im/Argo-floats/owc-python.svg)](https://gitter.im/Argo-floats/owc-python?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

This software is a python implementation of the "OWC" salinity calibration method used in Argo floats Delayed Mode Quality Control.

This software is in very active development and its usage may change any time. Please [post an issue to get involved if you're interested](https://github.com/euroargodev/argodmqc_owc/issues/new/choose).

# Installation

Download the code from this repository, using a git software or the command line:
```bash
git clone git@github.com:euroargodev/argodmqc_owc.git
```

Go to your local copy and then install dependencies by running the following command in terminal:

```bash
pip install -r requirements.txt
```

To test your installation before changing anything, run the code in your terminal by entering the ``argodmqc_owc`` directory and using the command:

```bash
python -W ignore owc_calibration.py
```

The -W ignore flags removes warnings.

The warnings are mostly due to comparing by NaN. Python complains about this, but it can be ignored in this case.

The code should run, but you should get some custom error messages that some WMO box data is missing (as the whole data set is not included in the download).

The software is currently set up just to run and test update salinity mapping.

# Software usage


You can click here to [![badge](https://img.shields.io/badge/launch-Pangeo%20binder-579ACA.svg?logo=data:image/png;base64,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)](https://binder.pangeo.io/v2/gh/euroargodev/argodmqc_owc/master?urlpath=lab/tree/docs/Usage.ipynb) and try this software before you even install it (thanks [Pangeo](pangeo.io)).

Otherwise, in a python script, you can calibrate salinity data this way:

```python
from ow_calibration.load_configuration.load_configuration import load_configuration
from ow_calibration.update_salinity_mapping.update_salinity_mapping import update_salinity_mapping

# store name for profile
FLOAT_NAME = "3901960"

# fetch the configuration and parameters set by the user
USER_CONFIG = load_configuration()

# get historical salinity data for comparison
update_salinity_mapping("/", FLOAT_NAME, USER_CONFIG)
```

### Parameters for your analysis

Most parameters for the analysis are set in the load configuration file.
Go to the load configuration file to set the paths of your wmo_boxes, and your historical data (climatology). You also may need to change the path to your float profile.

Outputs will be saved, by default, to the data folder in this project, but feel free to change this. If you keep the same output folder selected, the data will be saved over if you run the analysis again.

The update salinity mapping is simplified for testing purposes. The SAF code is removed, as well as potential vorticity code. So set these parameters to 0 until this code is pulled into the repo.

To run this code, set up your parameters in load_config, enter a float name, and enter a list of "/".

Eventually this will all be put into a `for` loop so multiple floats can be processed all at once

### Plots

The trajectory plotting routine is included, but commented out. This is because some of the dependencies of geopandas require libraries writen in C. If you have anaconda, geopandas can be installed using:

```bash
conda install geopandas
```

If you refuse to use anaconda, you can install using binary files. Check the documentation for geopandas online.

Then, to run the classic suite of OWC plots, simply type:

```python
from ow_calibration.plot_diagnostics.plot_diagnostics import plot_diagnostics
plot_diagnostics("/", FLOAT_NAME, USER_CONFIG)
```

# Software development plan

[...]

# Software history

- Migration of the code from from BODC/NOC git to euroargodev/argodmqc_owc. See https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues/10 for more details on the migration. Contribution from [G. Maze](https://github.com/gmaze)

- Alpha experts user testings with [feedbacks available here](https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues). Contributions from: [K. Walicka](https://github.com/kamwal), [C. Cabanes](https://github.com/cabanesc), [A. Wong](https://github.com/apswong)

- BODC created [the first version of the code](https://git.noc.ac.uk/bodc/owc-software-python), following the [Matlab implementation](https://github.com/ArgoDMQC/matlab_owc). Contributions from: [M. Donnelly](https://github.com/matdon17), [E. Small](https://github.com/edsmall-bodc), [K. Walicka](https://github.com/kamwal).