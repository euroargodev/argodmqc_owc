""" Prepare configuration for unit testing """
import os
from pyowc.configuration import load as load_configuration

import pyowc
pyowc_root = os.path.split(os.path.abspath(pyowc.__file__))[0]
data_root = pyowc_root.replace("pyowc", "data/")

TESTS_CONFIG = load_configuration()
TESTS_CONFIG['TEST_DIRECTORY'] = "data/test_data"
TESTS_CONFIG['TEST_FLOAT_SOURCE'] = "3901960"

for path in [k for k in TESTS_CONFIG.keys() if "DIRECTORY" in k]:
    TESTS_CONFIG[path] = os.path.abspath(TESTS_CONFIG[path].replace("data/", data_root))
