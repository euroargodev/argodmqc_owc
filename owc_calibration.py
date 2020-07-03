"""
-----ow calibration-----

Written by: Edward Small
When 01/05/2020

Python file to run OWC analysis on float profiles

Install dependencies by running the following command in terminal:

pip install -r requirements.txt

Before changing anything, run the code in your terminal by entering the
owc_software_python directory and using the command:

python -W ignore owc_calibration.py

The -W ignore flags removes warnings

The warnings are mostly due to comparing by NaN. Python complains about this,
but it can be ignored in this case

The code should run, but you should get some custom error messages that some wmo
box data is missing (as the whole data set is not included in the download)

Currently set up just to run and test update salinity mapping

Most parameters for the analysis are set in the load configuration file.
Go to the load configuration file to set the paths of your wmo_boxes, and
your historical data (climatology). You also may need to change the path to
your float profile.

Outputs will be saved, by default, to the data folder in this project,
but feel free to change this. If you keep the same output folder selected, the
data will be saved over if you run the analysis again.

The update salinity mapping is simplified for testing purposes. The SAF code is
removed, as well as potential vorticity code. So set these parameters to 0
until this code is pulled into the repo.

To run this code, set up your parameters in load_config, enter a float name,
and enter a list of "/".

Eventually this will all be put into a `for` loop so multiple floats can be
processed all at once

The trajectory plotting routine is included, but commented out. This is because some of the
dependencies of geopandas require libraries wirrten in C. If you have anaconda, geopandas can be
installed using:

conda install geopandas.

If you refuse to use anaconda, you can install using binary files. Check the documentation for
geopandas online

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

from ow_calibration.calc_piecewisefit.calc_piecewisefit import calc_piecewisefit
from ow_calibration.load_configuration.load_configuration import load_configuration
from ow_calibration.set_calseries.set_calseries import set_calseries
from ow_calibration.update_salinity_mapping.update_salinity_mapping import update_salinity_mapping
from ow_calibration.plot_diagnostics.plot_diagnostics import plot_diagnostics

# store name for profile
FLOAT_NAME = "3901960"

# fetch the configuration and parameters set by the user
USER_CONFIG = load_configuration()

# get historical salinity data for comparison
update_salinity_mapping("/", FLOAT_NAME, USER_CONFIG)

# get the calibration setting parameters
set_calseries("/", FLOAT_NAME, USER_CONFIG)

# calibrate the salinities
calc_piecewisefit("/", FLOAT_NAME, USER_CONFIG)

# create the diagnostic plots
plot_diagnostics("/", FLOAT_NAME, USER_CONFIG)
