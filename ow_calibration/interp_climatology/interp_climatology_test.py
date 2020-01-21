import unittest
import numpy as np
from ow_calibration.interp_climatology.interp_climatology import interp_climatology


class MyTestCase(unittest.TestCase):
    def test_something(self):
        grid_sal = np.array([[33.6620, 33.7630, 33.9060, 33.8380, 33.8650],
                             [np.inf, 33.7615, 33.9060, 33.8380, 33.8650],
                             [np.inf, 33.7630, 33.9060, 33.8380, 33.8660],
                             [np.inf, 33.7605, 33.9070, 33.8380, 33.8660],
                             [33.6623, 33.7620, 33.9065, 33.8390, 33.8660]])
        #grid_sal = np.full((5, 5), np.inf)
        grid_theta = np.array([[np.inf, 5.4031, 3.8174, 5.0465, 4.7226],
                               [1.2131, 5.3989, 3.8168, 5.0455, 4.7224],
                               [1.2124, 5.4023, 3.8163, 5.0434, 4.7223],
                               [1.2143, 5.3962, 3.8150, 5.0433, 4.7221],
                               [1.2120, 5.3851, 3.8165, 5.0442, 4.7230]])
        grid_pres = np.array([[12, 7, 9, 6, 6],
                             [13, 8, 10, 7, 8],
                             [14, 9, 11, 8, 10],
                             [15, 10, 15, 9, 12],
                             [np.inf, 11, 16, 10, 14]])
        float_sal = np.array([33.8030,
                              33.8050,
                              33.8030,
                              33.8050,
                              33.8270])
        float_theta = np.array([1.0769,
                                1.0758,
                                1.0794,
                                1.0609,
                                0.7645])
        float_pres = np.array([3.0000,
                               5.0000,
                               15.1000,
                               25.1000,
                               36.0000])

        interp_climatology(grid_sal, grid_theta, grid_pres, float_sal, float_theta, float_pres)

if __name__ == '__main__':
    unittest.main()
