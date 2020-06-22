from ow_calibration.fit_cond.fit_cond import fit_cond, nlbpfun
import numpy as np
import scipy.io as scipy
import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        fit_input = scipy.loadmat("data/test_data/fit_cond/fit_cond_input.mat")

        fit_cond(fit_input['x'], fit_input['y'], fit_input['n_err'],
                 fit_input['lvcov'], 'max_no_breaks', np.array([4]))


if __name__ == '__main__':
    unittest.main()
