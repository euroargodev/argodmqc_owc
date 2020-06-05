"""
-----Update Salinity Mapping Test File-----

Written by: Edward Small
When: 05/06/2020

Contains unit tests to check the functionality of the `sorter` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import math
import numpy as np
from ow_calibration.sorter.sorter import sorter


class SorterTestCase(unittest.TestCase):
    """
    Sorter test case
    """

    def setUp(self):
        """
        Set up repeated test values
        :return: Nothing
        """
        self.sites = A = np.array([-1, 0])
        self.msites = np.arange(-1, 1.01, 0.01)

    def test_returns_array(self):
        """
        Check that the sorter returns a numpy array
        :return: Nothing
        """

        print("Testing that sorter returns numpy array")

        sorted_test = sorter(self.sites, self.msites)

        self.assertEqual(type(sorted_test), np.ndarray,
                         "Return type for sorter is incorrect")

    def test_returns_length(self):
        """
        Check that the sorter returns a numpy array of the correct length
        :return: Nothing
        """

        print("Testing that sorter returns numpy array of correct length")

        sorted_test = sorter(self.sites, self.msites)

        self.assertEqual(sorted_test.__len__(), self.msites.__len__(),
                         "Return for sorter is incorrect length")

    def test_returns_ones(self):
        """
        Check that the sorter returns a numpy array of the correct length
        :return: Nothing
        """

        print("Testing that sorter returns numpy array of ones "
              "if boundaries are bad")

        sorted_test = sorter(np.array([-1, -1]), self.msites)

        for i in range(sorted_test.__len__()):
            self.assertEqual(sorted_test[i], 1, "Expected all ones")

    def test_returns_in_bounds(self):
        """
        Check that the sorter returns 0's for in bounds, 1's for out bounds
        :return: Nothing
        """

        print("Testing sorter boundaries")

        sorted_test = sorter(self.sites, self.msites)

        for i in range(math.floor(sorted_test.__len__()/2)):
            self.assertEqual(sorted_test[i], 0, "Expected 0's")

        for i in range(math.floor(sorted_test.__len__()/2), sorted_test.__len__()):
            self.assertEqual(sorted_test[i], 1, "Expected 1's")


if __name__ == '__main__':
    unittest.main()
