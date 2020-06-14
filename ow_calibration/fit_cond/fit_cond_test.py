from ow_calibration.fit_cond.fit_cond import fit_cond
import numpy as np
import scipy.io as scipy
import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        input = scipy.loadmat("data/test_data/fit_cond/fit_cond_input.mat")

        fit_cond(input['x'], input['y'], input['n_err'],
                 input['lvcov'], 'max_no_breaks', np.array([4]))




if __name__ == '__main__':
    unittest.main()
