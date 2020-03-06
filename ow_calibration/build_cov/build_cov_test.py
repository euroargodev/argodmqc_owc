import unittest
import numpy as np
from ow_calibration.build_cov.build_cov import build_cov


class MyTestCase(unittest.TestCase):
    def test_something(self):
        config = {"MAPSCALE_LONGITUDE_SMALL": 4,
                  "MAPSCALE_LATITUDE_SMALL": 2,
                  "MAPSCALE_PHU_SMALL": 0.1,
                  "MAP_USE_PV": 0}
        ptmp = np.array([[0.7058, 0.7039, 0.8285],
                         [0.6713, 0.6664, 0.7432],
                         [0.8257, 0.8224, 0.7804],
                         [0.7452, 0.7411, 1.1980],
                         [0.7836, 0.7802, 1.1504],
                         [1.2008, 1.2010, 1.2497],
                         [1.1496, 1.1481, 1.3036],
                         [1.2520, 1.2553, 1.0921],
                         [1.3039, 1.3046, np.nan],
                         [1.0947, 1.0962, np.nan]])
        coord_float = np.array([[0.0572, -0.0592, 5.1083],
                                [0.0578, -0.0591, 5.0993],
                                [0.0586, -0.0585, 5.0861]]) * 1.0e+03

        build_cov_vec = np.vectorize(build_cov)

        build_cov(ptmp, coord_float, config)


if __name__ == '__main__':
    unittest.main()
