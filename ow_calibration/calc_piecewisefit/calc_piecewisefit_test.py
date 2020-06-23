import unittest
from ow_calibration.calc_piecewisefit.calc_piecewisefit import calc_piecewisefit
from ow_calibration.load_configuration.load_configuration import load_configuration

class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.float_source = "3901960"
        calc_piecewisefit("/", self.float_source, load_configuration())


if __name__ == '__main__':
    unittest.main()
