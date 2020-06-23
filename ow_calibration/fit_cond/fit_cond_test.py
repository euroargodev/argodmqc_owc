"""
-----fit condition Test File-----

Written by: Edward Small
When: 31/10/2019

Contains unit tests to check the functionality of the `fit_cond` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python

"""

import unittest
import numpy as np
import scipy.io as scipy
from ow_calibration.fit_cond.fit_cond import fit_cond


class FitCondTestCase(unittest.TestCase):
    """
    Test cases for 'fit_cond' function
    """
    def setUp(self):
        fit_input = scipy.loadmat("data/test_data/fit_cond/fit_cond_input.mat")
        fit_out = scipy.loadmat("data/test_data/fit_cond/fit_cond_output.mat")

        self.in_x = fit_input['x']
        self.in_y = fit_input['y']
        self.in_err = fit_input['n_err']
        self.in_cov = fit_input['lvcov']

        self.xfit = fit_out['xfit']
        self.condslope = fit_out['condslope']
        self.condslope_err = fit_out['condslope_err']
        self.time_derivself = fit_out['time_deriv']
        self.time_deriv_err = fit_out['time_deriv_err']
        self.sta_mean = fit_out['sta_mean']
        self.sta_rms = fit_out['sta_rms']
        self.ndf = fit_out['NDF']
        self.fit_coef = fit_out['fit_coef']
        self.fit_breaks = fit_out['fit_breaks']

    def test_max_no_breaks(self):
        """
        Check return values if we specify a maximum number of breaks
        :return: nothing
        """
        print("Testing fit_cond for max_no_breaks")

        python_test = fit_cond(self.in_x, self.in_y, self.in_err,
                               self.in_cov, 'max_no_breaks', 4)

        for i in range(python_test[0].__len__()):
            self.assertEqual(python_test[0][i], self.xfit[i],
                             "xfit is incorrect")

        for i in range(python_test[2][0].__len__()):
            self.assertAlmostEqual(python_test[2][0, i], 0.000122830414419824, 12,
                                   "slope error is incorrect")

        for i in range(python_test[5].shape[0]):
            for j in range(python_test[5].shape[1]):
                if ~np.isnan(python_test[5][i, j]):
                    self.assertAlmostEqual(python_test[5][i, j], self.sta_mean[i, j], 12,
                                           "mean is incorrect")

        for i in range(python_test[6].shape[0]):
            for j in range(python_test[6].shape[1]):
                if ~np.isnan(python_test[6][i, j]):
                    self.assertAlmostEqual(python_test[6][i, j], self.sta_rms[i, j], 12,
                                           "rms is incorrect")

        self.assertEqual(python_test[7], self.ndf, "degrees of freeodm is incorrect")


if __name__ == '__main__':
    unittest.main()
