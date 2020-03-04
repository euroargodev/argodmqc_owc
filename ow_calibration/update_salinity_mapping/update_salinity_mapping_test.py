import unittest
from ow_calibration.load_configuration.load_configuration import load_configuration
from ow_calibration.update_salinity_mapping.update_salinity_mapping import update_salinity_mapping


class MyTestCase(unittest.TestCase):
    def test_something(self):
        update_salinity_mapping("/", "3901960", load_configuration())


if __name__ == '__main__':
    unittest.main()
