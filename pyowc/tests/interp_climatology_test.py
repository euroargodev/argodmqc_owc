"""
-----Interpolate Climatology Test File-----

Written by: Edward Small
When: 27/01/2020

Contains unit tests to check the functionality of the `interp_climatology` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from scipy.io import loadmat
from pyowc.data.wrangling import interp_climatology

#pylint: disable=too-many-instance-attributes
class InterpClimatologyTestCase(unittest.TestCase):
    """
    Test cases for interp_climatology function
    """

    def setUp(self):
        """
        Set up some values for testing, pulling from matlab matrices.
        The matlab matrices were used in the matlab version of this code, so answers should match
        :return: Nothing
        """
        # load in the data for testing
        test = loadmat('../../data/test_data/interp_climatology/testfile.mat')
        results = loadmat('../../data/test_data/interp_climatology/results.mat')

        # set the test variables from the loaded .mat file
        self.grid_sal = test['S']
        self.grid_theta = test['Theta']
        self.grid_pres = test['P']
        self.float_sal = test['S_f']
        self.float_theta = test['Theta_f']
        self.float_pres = test['P_f']

        self.expected_interp_pres = results['P_h']
        self.expected_interp_sal = results['S_h']

    def test_throws_error(self):
        """
        Test that an error is thrown if we have no good climatology data
        :return: Nothing
        """
        print("Testing that interp_climatology throws an error if all data is bad")

        bad_grid_sal = np.full((5, 5), np.inf)

        with self.assertRaises(ValueError) as bad_climatology:
            interp_climatology(bad_grid_sal, self.grid_theta, self.grid_pres,
                               self.float_sal, self.float_theta, self.float_pres)

        self.assertTrue("No good climatological data has been found"
                        in str(bad_climatology.exception))

    def test_returns_correct_shape(self):
        """
        Test that the returned matrix is the shape we expect it to be
        :return: Nothing
        """
        print("Testing that interp_climatology returns a matrix of the correct shape")

        sal, pres = interp_climatology(self.grid_sal, self.grid_theta, self.grid_pres,
                                       self.float_sal, self.float_theta, self.float_pres)

        self.assertTrue(sal.shape == self.expected_interp_sal.shape,
                        "salinity matrix shape is incorrect")
        self.assertTrue(pres.shape == self.expected_interp_pres.shape,
                        "pressure matrix shape is incorrect")

    def test_returns_same_shape(self):
        """
        Test that the salinity and pressures matrices are equal shapes
        :return: Nothing
        """
        print("Testing that interp_climatology returns a matrices of the same shape")

        sal, pres = interp_climatology(self.grid_sal, self.grid_theta, self.grid_pres,
                                       self.float_sal, self.float_theta, self.float_pres)

        self.assertTrue(sal.shape == pres.shape,
                        "salinity and pressure matrices shaped differently")

    def test_returns_correct_values(self):
        """
        Test that the output of this function matches the matlab version value by value.
        We cannot compare NaNs directly (as it will always give false, so only do a direct
        comparison if they are both not NaN
        :return: Nothing
        """
        print("Testing that interp_climatology returns matrices with correct values")

        sal, pres = interp_climatology(self.grid_sal, self.grid_theta, self.grid_pres,
                                       self.float_sal, self.float_theta, self.float_pres)

        for i in range(0, sal.shape[0]):
            for j in range(0, sal.shape[1]):
                if not (np.isnan(sal[i, j]) and np.isnan(self.expected_interp_sal[i, j])):
                    self.assertTrue(sal[i, j] == self.expected_interp_sal[i, j],
                                    ("Values at ", i, " and ", j, " do not match for salinity"))
                    self.assertTrue(pres[i, j] == self.expected_interp_pres[i, j],
                                    ("Values at ", i, " and ", j, " do not match for pressure"))


if __name__ == '__main__':
    unittest.main()
