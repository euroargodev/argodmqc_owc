
import os
import sys
import io
import unittest
import tempfile
import math
import numpy as np
from scipy.io import loadmat

from pyowc.configuration import load as load_configuration
from pyowc.calibration import update_salinity_mapping, calc_piecewisefit
from . import TESTS_CONFIG


#pylint: disable=bare-except
class UpdateSalinityMapping(unittest.TestCase):
    """
    Test cases for update_salinity_mapping function
    """

    def setUp(self):
        """
        Only run if we are missing our test file
        :return: Nothing
        """
        # self.float_source = "3901960"
        # self.python_output_path = "../../data/float_mapped/map_" + self.float_source + ".mat"
        # matlab_output_path =     "../../data/test_data/float_mapped_test/map_" + self.float_source + ".mat"
        # DEFAULT_CONFIG = load_configuration()
        # for path in ['CONFIG_DIRECTORY', 'FLOAT_CALIB_DIRECTORY', 'FLOAT_SOURCE_DIRECTORY',
        #              'FLOAT_MAPPED_DIRECTORY', 'HISTORICAL_DIRECTORY']:
        #     DEFAULT_CONFIG[path] = DEFAULT_CONFIG[path].replace("data/", "../../data/")
        # self.config = DEFAULT_CONFIG
        self.float_source = TESTS_CONFIG['TEST_FLOAT_SOURCE']
        self.python_output_path = os.path.sep.join([TESTS_CONFIG['FLOAT_MAPPED_DIRECTORY'],
                                                    TESTS_CONFIG['FLOAT_MAPPED_PREFIX'] +
                                                    TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                                    TESTS_CONFIG['FLOAT_MAPPED_POSTFIX']])
        matlab_output_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "float_mapped_test",
                                                    TESTS_CONFIG['FLOAT_MAPPED_PREFIX'] +
                                                    TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                                    TESTS_CONFIG['FLOAT_MAPPED_POSTFIX']])
        self.config = TESTS_CONFIG

        if not os.path.exists(self.python_output_path):
            print("Getting mapped data for testing...")
            update_salinity_mapping("/", self.config, self.float_source)

        self.matlab_output = loadmat(matlab_output_path)
        self.python_output = loadmat(self.python_output_path)
        self.acceptable_diff = 3

    def test_salinity_mapping(self):
        """
        Check that the salinity mapping protocol runs
        :return: Nothing
        """

        print("testing that update salinity mapping runs through")

        # use temporary directory for output to ensure the file doesn't exist
        with tempfile.TemporaryDirectory() as tmpdir:
            # override the output directory for this call
            tmp_config = {**self.config, "FLOAT_MAPPED_DIRECTORY": tmpdir}

            try:
                update_salinity_mapping("/", tmp_config, self.float_source)
            except:
                self.fail("Update salinity mapping encountered an unexpected error")

        # Should use precalculated data
        captured_output = io.StringIO()
        sys.stdout = captured_output
        update_salinity_mapping("/", self.config, self.float_source)
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


class CalcPiecewiseFit(unittest.TestCase):
    """
    Test cases for 'calc_piecewisefit' function
    """

    def test_custom(self):
        """
        Change variables in this test to use different mapped outputs
        :return: nothing
        """
        calc_piecewisefit("/", TESTS_CONFIG['TEST_FLOAT_SOURCE'], TESTS_CONFIG)

        test_data_path = os.path.sep.join([TESTS_CONFIG['FLOAT_CALIB_DIRECTORY'],
                                      TESTS_CONFIG['FLOAT_CALIB_PREFIX'] +
                                      TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                      TESTS_CONFIG['FLOAT_CALIB_POSTFIX']])
        test = loadmat(test_data_path)

        matlab_data_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], 'float_calib_test',
                                             TESTS_CONFIG['FLOAT_CALIB_PREFIX'] +
                                             TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                             TESTS_CONFIG['FLOAT_CALIB_POSTFIX']])
        matlab = loadmat(matlab_data_path)

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
