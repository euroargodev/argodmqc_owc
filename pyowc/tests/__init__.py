""" Prepare configuration for unit testing """
import os
from pyowc.configuration import load as load_configuration


TESTS_CONFIG = load_configuration()
TESTS_CONFIG['TEST_DATA'] = "data/test_data"
TESTS_CONFIG['TEST_FLOAT_SOURCE'] = "3901960"

for path in ['TEST_DATA', 'CONFIG_DIRECTORY', 'FLOAT_CALIB_DIRECTORY',
             'FLOAT_SOURCE_DIRECTORY', 'FLOAT_MAPPED_DIRECTORY']:
    TESTS_CONFIG[path] = os.path.abspath(TESTS_CONFIG[path].replace("data/", "../../data/"))
