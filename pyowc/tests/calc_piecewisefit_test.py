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
from scipy.io import loadmat
import pyowc
from pyowc import core

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
        config = pyowc.configuration.load()
        for path in ['CONFIG_DIRECTORY', 'FLOAT_CALIB_DIRECTORY', 'FLOAT_SOURCE_DIRECTORY', 'FLOAT_MAPPED_DIRECTORY']:
            config[path] = config[path].replace("data/", "../../data/")

        core.stats.calc_piecewisefit("/", float_source, config)

        test = loadmat("../../data/float_calib/cal_" + float_source + ".mat")
        matlab = loadmat("../../data/test_data/float_calib_test/cal_" +
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
