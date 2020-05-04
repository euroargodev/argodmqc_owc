import unittest
import os.path as path
import numpy as np
import scipy.io as scipy
from ow_calibration.load_configuration.load_configuration import load_configuration
from ow_calibration.update_salinity_mapping.update_salinity_mapping import update_salinity_mapping


class MyTestCase(unittest.TestCase):

    def setUp(self):
        """
        Only run if we are missing our test file
        :return: Nothing
        """

        if not path.exists("data/float_mapped/map_3901960.mat"):
            print("Getting mapped data for testing...")
            update_salinity_mapping("/", "3901960", load_configuration())


    def test_run_salinity_mapping(self):
        """
        We need to run update salinity mapping to check it runs through.
        Save the output so we can test it
        :return: Nothing
        """

        print("Testing that salinity mapping runs smoothly...")
        update_salinity_mapping("/", "3901960", load_configuration())

    def test_ptmp_output(self):
        """
        check that ptmp matrices match across version
        :return: Nothing
        """

        test = scipy.loadmat("data/test_data/float_mapped_test/map_3901960.mat")['la_ptmp']
        result = scipy.loadmat("data/float_mapped/map_3901960.mat")['la_ptmp']

        indices = test.shape
        for i in range(indices[0]):
            for j in range(indices[1]):
                if not (np.isnan(test[i, j]) and np.isnan(result[i, j])):
                    self.assertEqual(test[i, j], result[i, j])

    def test_ptmp_output(self):
        """
        check that ptmp matrices match across version
        :return: Nothing
        """

        test = scipy.loadmat("data/test_data/float_mapped_test/map_3901960.mat")['la_ptmp']
        result = scipy.loadmat("data/float_mapped/map_3901960.mat")['la_ptmp']

        indices = test.shape
        for i in range(indices[0]):
            for j in range(indices[1]):
                if not (np.isnan(test[i, j]) and np.isnan(result[i, j])):
                    self.assertEqual(test[i, j], result[i, j])


if __name__ == '__main__':
    unittest.main()
