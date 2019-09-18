import unittest
from ow_calibration.lo_system_configuration.ow_config import lo_system_configuration


class lo_system_configuration_test(unittest.TestCase):

    def test_all_keys_have_value(self):
        for i in lo_system_configuration:
            self.assertNotEqual(lo_system_configuration[i], '')


if __name__ == '__main__':
    unittest.main()
