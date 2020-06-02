# Argo float salinity calibration software
## OWC in Python

![build](https://github.com/euroargodev/argodmqc_owc/workflows/build/badge.svg)
[![codecov](https://codecov.io/gh/euroargodev/argodmqc_owc/branch/master/graph/badge.svg)](https://codecov.io/gh/euroargodev/argodmqc_owc)
[![Requirements Status](https://requires.io/github/euroargodev/argodmqc_owc/requirements.svg?branch=master)](https://requires.io/github/euroargodev/argodmqc_owc/requirements/?branch=master)

[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)


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

[...]

# Software development plan

[...]

# Software history

- Migration of the code from from BODC/NOC git to euroargodev/argodmqc_owc. See https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues/10 for more details on the migration. Contribution from @gmaze.

- Alpha experts user testings with [feedbacks available here](https://github.com/euroargodev/User-Acceptance-Test-Python-version-of-the-OWC-tool/issues). Contributions from: @kamwal, @cabanesc, @apswong.

- BODC created [the first version of the code](https://git.noc.ac.uk/bodc/owc-software-python), following the [Matlab implementation](https://github.com/ArgoDMQC/matlab_owc). Contributions from: @matdon17, @edsmall-bodc, @kamwal.