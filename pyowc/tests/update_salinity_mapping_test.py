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
import io
import sys
import numpy as np
from scipy.io import loadmat

from pyowc.configuration import load as load_configuration
from pyowc.calibration import update_salinity_mapping

#pylint: disable=bare-except
class MyTestCase(unittest.TestCase):
    """
    Test cases for update_salinity_mapping function
    """

    def setUp(self):
        """
        Only run if we are missing our test file
        :return: Nothing
        """
        self.float_source = "3901960"
        self.python_output_path = "../../data/float_mapped/map_" + self.float_source + ".mat"
        matlab_output_path = "../../data/test_data/float_mapped_test/map_" + self.float_source + ".mat"
        DEFAULT_CONFIG = load_configuration()
        for path in ['CONFIG_DIRECTORY', 'FLOAT_CALIB_DIRECTORY', 'FLOAT_SOURCE_DIRECTORY',
                     'FLOAT_MAPPED_DIRECTORY', 'HISTORICAL_DIRECTORY']:
            DEFAULT_CONFIG[path] = DEFAULT_CONFIG[path].replace("data/", "../../data/")
        self.config = DEFAULT_CONFIG

        if not os.path.exists(self.python_output_path):
            print("Getting mapped data for testing...")
            update_salinity_mapping("/", self.float_source, self.config)

        self.matlab_output = loadmat(matlab_output_path)
        self.python_output = loadmat(self.python_output_path)
        self.acceptable_diff = 3

    def test_salinity_mapping(self):
        """
        Check that the salinity mapping protocol runs
        :return: Nothing
        """

        print("testing that update salinity mapping runs through")

        if os.path.exists(self.python_output_path):
            os.remove(self.python_output_path)

        try:
            update_salinity_mapping("/", self.float_source, self.config)

        except:
            self.fail("Update salinity mapping encountered an unexpected error")

        # Should use precalculated data
        captured_output = io.StringIO()
        sys.stdout = captured_output
        update_salinity_mapping("/", self.float_source, self.config)
        sys.stdout = sys.__stdout__

        self.assertTrue("Using precalculated data" in captured_output.getvalue(),
                        "Should use precalculated data !")


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
                               "error: mean of potential temperatures differ "
                               "between python and matlab")

    def test_mapped_sal_output(self):
        """
        check that mapped salinity matrices match across version
        :return: Nothing
        """

        print("Testing update_salinity_mapping gives correct mapped salinity")

        test = self.matlab_output['la_mapped_sal']
        result = self.python_output['la_mapped_sal']
        test_mean = np.nanmean(test)
        result_mean = np.nanmean(result)

        self.assertAlmostEqual(test_mean, result_mean, self.acceptable_diff,
                               "error: mean of mapped salinities "
                               "differ between python andd matlab")

    def test_mapped_salerrors_output(self):
        """
        check that salinity errors matrices match across version
        :return: Nothing
        """

        print("Testing update_salinity_mapping gives correct salinity errors")

        test = self.matlab_output['la_mapsalerrors']
        result = self.python_output['la_mapsalerrors']
        test_mean = np.nanmean(test)
        result_mean = np.nanmean(result)

        self.assertAlmostEqual(test_mean, result_mean, self.acceptable_diff,
                               "error: mean of salinity errors "
                               "differ between python andd matlab")

    def test_noise_sal_output(self):
        """
        check that noise salinity matrices match across version
        :return: Nothing
        """

        print("Testing update_salinity_mapping gives correct noise salinity")

        test = self.matlab_output['la_noise_sal']
        result = self.python_output['la_noise_sal']
        test_mean = np.nanmean(test)
        result_mean = np.nanmean(result)

        self.assertAlmostEqual(test_mean, result_mean, self.acceptable_diff,
                               "error: mean of noise salinity "
                               "differ between python andd matlab")

    def test_signal_sal_output(self):
        """
        check that signal salinity matrices match across version
        :return: Nothing
        """

        print("Testing update_salinity_mapping gives correct signal salinity")

        test = self.matlab_output['la_signal_sal']
        result = self.python_output['la_signal_sal']
        test_mean = np.nanmean(test)
        result_mean = np.nanmean(result)

        self.assertAlmostEqual(test_mean, result_mean, self.acceptable_diff,
                               "error: mean of signal salinity "
                               "differ between python andd matlab")


if __name__ == '__main__':
    unittest.main()
