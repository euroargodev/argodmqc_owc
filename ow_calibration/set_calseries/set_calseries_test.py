import unittest
from ow_calibration.set_calseries.set_calseries import set_calseries


class MyTestCase(unittest.TestCase):
    def test_something(self):
        float_dir = "/"
        float_name = "3901960"
        system_config = {'FLOAT_SOURCE_POSTFIX': ".mat",
                         'FLOAT_CALSERIES_PREFIX': "calseries_",
                         'FLOAT_CALIB_POSTFIX': ".mat"}

        set_calseries(float_dir, float_name, system_config)



if __name__ == '__main__':
    unittest.main()
