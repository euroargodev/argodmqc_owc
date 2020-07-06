"""
-----Set Calseries Test File-----

Written by: Edward Small
When: 21/11/2019

Contains unit tests to check the functionality of the `set_calseries` function

Changes can be made to this file to test all different kinds of floats, as long
as you have data for comparison

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import os
import scipy.io as scipy
from ow_calibration.set_calseries.set_calseries import set_calseries


class SetCalSeriesTestCase(unittest.TestCase):
    """
    Test cases forset_calseries function
    """

    def setUp(self):
        """
        Only run if we are missing our test file
        :return: Nothing
        """

        self.float_source = "3901960"
        self.float_dir = "/"
        self.system_config = {'FLOAT_SOURCE_DIRECTORY': "data/float_source",
                              'FLOAT_SOURCE_POSTFIX': ".mat",
                              'FLOAT_CALIB_DIRECTORY': "data/float_calib",
                              'FLOAT_CALSERIES_PREFIX': "calseries_",
                              'FLOAT_CALIB_POSTFIX': ".mat"}
        self.python_output_path = "data/float_calib/calseries_" + \
                                  self.float_source + ".mat"
        matlab_output_path = "data/test_data/float_calib_test/calseries_" + \
                             self.float_source + ".mat"

        if not os.path.exists(self.python_output_path):
            print("Getting calibrated data for testing...")
            set_calseries(self.float_dir, self.float_source, self.system_config)

        self.matlab_output = scipy.loadmat(matlab_output_path)
        self.python_output = scipy.loadmat(self.python_output_path)

    def test_set_calseries(self):
        """
        Check that the function runs and saves
        :return: Nothing
        """
        print("Testing that set_calseries saves file")

        if os.path.exists(self.python_output_path):
            os.remove(self.python_output_path)

        try:
            set_calseries(self.float_dir, self.float_source, self.system_config)
            self.assertTrue(os.path.exists(self.python_output_path),
                            "calseries file was not saved")

        except KeyError:
            self.fail("Update salinity mapping encountered an unexpected error")

        # should use old file
        set_calseries(self.float_dir, self.float_source, self.system_config)

    def test_bad_boundaries(self):
        """
        Check that, even if we pass  in bad boundaries, we still get values
        :return: nothing
        """
        print("Testing that set_calseries sets parameters even if boundaries are bad")

        if os.path.exists(self.python_output_path):
            os.remove(self.python_output_path)

        set_calseries(self.float_dir, self.float_source, self.system_config)

        self.assertTrue(os.path.exists(self.python_output_path))


if __name__ == '__main__':
    unittest.main()
