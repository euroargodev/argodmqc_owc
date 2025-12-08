"""Prepare configuration for unit testing"""

import os

import pyowc
from pyowc.utilities import load_configuration_from_json_file

print(os.getcwd(), os.listdir())
pyowc_root = os.path.split(os.path.abspath(pyowc.__file__))[0]
data_root = pyowc_root.replace("pyowc", "data/")

TESTS_CONFIG = load_configuration_from_json_file("pyowc/tests/test_config.json")

TESTS_CONFIG["TEST_DIRECTORY"] = "data/test_data"
TESTS_CONFIG["FLOAT_PLOTS_DIRECTORY"] = "data/test_data/float_plots/"
TESTS_CONFIG["FLOAT_PLOTS_FORMAT"] = "eps"

TESTS_CONFIG["FLOAT_SOURCE_DIRECTORY"] = "data/test_data/float_source/"

# Define specific values used in tests that depends on a specific float processing:
TESTS_CONFIG["TEST_FLOAT_SOURCE"] = "3901960"
TESTS_CONFIG["TEST_FLOAT_WMO_BOXES"] = [3505, 3506]
TESTS_CONFIG["TEST_FLOAT_WMO_BOXES_Nhist"] = [10, 33, 787]  # Nb of CTD, Bottle and Argo measurements in 3505

TESTS_CONFIG["TEST_FLOAT_IN_HIST"] = "1900193"
TESTS_CONFIG["TEST_FLOAT_N_TO_REMOVE"] = 30
TESTS_CONFIG["TEST_FLOAT_N_DATA"] = 830

# Fix paths for tests:
for path in [k for k in TESTS_CONFIG.keys() if "DIRECTORY" in k]:
    TESTS_CONFIG[path] = os.path.abspath(TESTS_CONFIG[path].replace("data/", data_root, 1))
