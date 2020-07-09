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

owc.calibration.update_salinity_mapping("/", FLOAT_NAME, USER_CONFIG)
```

### Parameters for your analysis

Parameters for the analysis are set in a configuration file. 
The configuration file has the same format as the Matlab software.

### Plots

```python
owc.plot.dashboard("/", FLOAT_NAME, USER_CONFIG)
```

# Software history

- Major refactoring of the software for performance optimisation and to fully embrace the Pythonic way of doing this ! 

- Phase 3

- Phase 2

- Phase 1

- Migration of the code from from BODC/NOC git to euroargodev/argodmqc_owc. See https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues/10 for more details on the migration. Contribution from [G. Maze](https://github.com/gmaze)

- Alpha experts user testings with [feedbacks available here](https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues). Contributions from: [K. Walicka](https://github.com/kamwal), [C. Cabanes](https://github.com/cabanesc), [A. Wong](https://github.com/apswong)

- BODC created [the first version of the code](https://git.noc.ac.uk/bodc/owc-software-python), following the [Matlab implementation](https://github.com/ArgoDMQC/matlab_owc). Contributions from: [M. Donnelly](https://github.com/matdon17), [E. Small](https://github.com/edsmall-bodc), [K. Walicka](https://github.com/kamwal).