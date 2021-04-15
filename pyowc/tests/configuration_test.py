import os
import unittest
from scipy.io import loadmat

from pyowc.configuration import load as load_configuration, set_calseries, print_cfg
from . import TESTS_CONFIG


class PrintCFG(unittest.TestCase):
    """ Test print_cfg """

    def setUp(self):
        self.config = TESTS_CONFIG

    def test_is_str(self):
        self.assertTrue(isinstance(print_cfg(self.config), str))



class LoadConfiguration(unittest.TestCase):
    """
    Test cases for load_configuration function
    """

    def setUp(self):
        self.lo_system_configuration = load_configuration()

    def test_is_dict(self):
        """
        Check return type is dictionary
        :return: Nothing
        """
        print("Testing that load_configuration() returns a dictionary...\n")
        self.assertTrue(isinstance(self.lo_system_configuration, dict))

    def test_contains_correct_number_of_keys(self):
        """
        Check returned dictionary contains the correct number of
        pairs
        :return: Nothing
        """
        print("Testing length of dictionary...\n")
        self.assertEqual(len(self.lo_system_configuration), 33)

    def test_all_keys_have_value(self):
        """
        Check that each pair is non-empty
        :return: Nothing
        """
        print("Testing that the keys in the dictionary all have a value...\n")
        for i in self.lo_system_configuration:
            #print(i, ": ", lo_system_configuration[i]) #uncomment this line to find error location
            self.assertNotEqual(self.lo_system_configuration[i], '')


class SetCalSeries(unittest.TestCase):
    """
    Test cases forset_calseries function
    """

    def setUp(self):
        """
        Only run if we are missing our test file
        :return: Nothing
        """

        self.float_source = TESTS_CONFIG['TEST_FLOAT_SOURCE']
        self.float_dir = "/"
        self.python_output_path = os.path.sep.join([TESTS_CONFIG['FLOAT_CALIB_DIRECTORY'],
                                                    TESTS_CONFIG['FLOAT_CALSERIES_PREFIX'] +
                                                    TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                                    TESTS_CONFIG['FLOAT_CALIB_POSTFIX']])
        matlab_output_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "float_calib_test",
                                                    TESTS_CONFIG['FLOAT_CALSERIES_PREFIX'] +
                                                    TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                                    TESTS_CONFIG['FLOAT_CALIB_POSTFIX']])
        self.system_config = TESTS_CONFIG

        if not os.path.exists(self.python_output_path):
            print("Getting calibrated data for testing...")
            set_calseries(self.float_dir, self.float_source, self.system_config)

        self.matlab_output = loadmat(matlab_output_path)
        self.python_output = loadmat(self.python_output_path)

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
