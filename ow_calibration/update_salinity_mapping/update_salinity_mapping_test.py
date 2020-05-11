"""
-----Update Salinity Mapping Test File-----

Written by: Edward Small
When: 21/11/2019

Contains unit tests to check the functionality of the `update_salinity_mapping` function

Changes can be made to this file to test all different kinds of floats, as long
as you have data for comparison

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import os
import numpy as np
import scipy.io as scipy
from ow_calibration.load_configuration.load_configuration import load_configuration
from ow_calibration.update_salinity_mapping.update_salinity_mapping import update_salinity_mapping


class MyTestCase(unittest.TestCase):
    """
    Test cases for update_salinity_mapping function
    """

    def setUp(self):
        """
        Only run if we are missing our test file
        :return: Nothing
        """

        float_source = "3901960"
        python_output_path = "data/float_mapped/map_" + float_source + ".mat"
        matlab_output_path = "data/test_data/float_mapped_test/map_" + float_source + ".mat"

        if not os.path.exists(python_output_path):
            print("Getting mapped data for testing...")
            update_salinity_mapping("/", float_source, load_configuration())

        self.matlab_output = scipy.loadmat(matlab_output_path)
        self.python_output = scipy.loadmat(python_output_path)
        self.acceptable_diff = 3

    def test_ptmp_output(self):
        """
        check that ptmp matrices match across version
        :return: Nothing
        """

        print("Testing update_salinity_mapping gives correct potential temperature")

        test = self.matlab_output['la_ptmp']
        result = self.python_output['la_ptmp']

        test_mean = np.nanmean(test)
        result_mean = np.nanmean(result)

        self.assertAlmostEqual(test_mean, result_mean, self.acceptable_diff,
                               "error: mean of potential temperatures differ between python andd matlab")

    def test_mapped_sal_output(self):
        """
        check that ptmp matrices match across version
        :return: Nothing
        """

        print("Testing update_salinity_mapping gives correct mapped salinity")

        test = self.matlab_output['la_mapped_sal']
        result = self.python_output['la_mapped_sal']
        test_mean = np.nanmean(test)
        result_mean = np.nanmean(result)

        self.assertAlmostEqual(test_mean, result_mean, self.acceptable_diff,
                               "error: mean of mapped salinities differ between python andd matlab")

    def test_mapped_salerrors_output(self):
        """
        check that ptmp matrices match across version
        :return: Nothing
        """

        print("Testing update_salinity_mapping gives correct mapped salinity")

        test = self.matlab_output['la_mapsalerrors']
        result = self.python_output['la_mapsalerrors']
        test_mean = np.nanmean(test)
        result_mean = np.nanmean(result)

        self.assertAlmostEqual(test_mean, result_mean, self.acceptable_diff,
                               "error: mean of mapped salinities differ between python andd matlab")

    def test_noise_sal_output(self):
        """
        check that ptmp matrices match across version
        :return: Nothing
        """

        print("Testing update_salinity_mapping gives correct mapped salinity")

        test = self.matlab_output['la_noise_sal']
        result = self.python_output['la_noise_sal']
        test_mean = np.nanmean(test)
        result_mean = np.nanmean(result)

        self.assertAlmostEqual(test_mean, result_mean, self.acceptable_diff,
                               "error: mean of mapped salinities differ between python andd matlab")

    def test_signal_sal_output(self):
        """
        check that ptmp matrices match across version
        :return: Nothing
        """

        print("Testing update_salinity_mapping gives correct mapped salinity")

        test = self.matlab_output['la_signal_sal']
        result = self.python_output['la_signal_sal']
        test_mean = np.nanmean(test)
        result_mean = np.nanmean(result)

        self.assertAlmostEqual(test_mean, result_mean, self.acceptable_diff,
                               "error: mean of mapped salinities differ between python andd matlab")


if __name__ == '__main__':
    unittest.main()
