import unittest
import numpy as np
from ow_calibration.build_cov.covarxy_pv.covarxy_pv import covarxy_pv


class MyTestCase(unittest.TestCase):
    def test_something(self):
        input_coords = np.array([0.0572, -0.0592, 5.1083]) * 1.0e+03
        coords = np.array([[0.0572, -0.0592, 5.1083],
                           [0.0578, -0.0591, 5.0993],
                           [0.0586, -0.0585, 5.0861]]) * 1.0e+03
        long = 4
        lat = 2
        phi = 0.1
        pv = 0

        covarxy_pv(input_coords, coords, long, lat, phi, pv)


if __name__ == '__main__':
    unittest.main()
