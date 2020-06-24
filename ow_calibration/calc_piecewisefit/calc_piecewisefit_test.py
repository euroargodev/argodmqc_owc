"""
-----calculate piecewise fit Test File-----

Written by: Edward Small
When: 14/06/2020

Contains unit tests to check the functionality of the `fit_cond` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python

"""

import unittest
import numpy as np
import scipy.io as scipy
from ow_calibration.calc_piecewisefit.calc_piecewisefit import calc_piecewisefit
from ow_calibration.load_configuration.load_configuration import load_configuration


class MyTestCase(unittest.TestCase):
    """
    Test cases for 'calc_piecewisefit' function
    """

    def test_custom(self):
        """
        Change variables in this test to use different mapped outputs
        :return: nothing
        """
        float_source = "3901960"
        calc_piecewisefit("/", float_source, load_configuration())

        test = scipy.loadmat("data/float_calib/cal_" + float_source + ".mat")
        matlab = scipy.loadmat("data/test_data/float_calib_test/cal_" +
                               float_source + ".mat")

        python_sal = test['cal_SAL']
        matlab_sal = matlab['cal_SAL']

        self.assertEqual(python_sal.shape, matlab_sal.shape)

        for i in range(python_sal.shape[0]):
            for j in range(python_sal.shape[1]):
                if ~np.isnan(python_sal[i, j]):
                    self.assertAlmostEqual(python_sal[i, j], matlab_sal[i, j], 3)

        python_sal_err = test['cal_SAL_err']
        matlab_sal_err = matlab['cal_SAL_err']

        for i in range(python_sal_err.shape[0]):
            for j in range(python_sal_err.shape[1]):
                if ~np.isnan(python_sal_err[i, j]):
                    self.assertAlmostEqual(python_sal_err[i, j], matlab_sal_err[i, j], 3)


if __name__ == '__main__':
    unittest.main()
