"""
-----Load Configuration Test File-----

Written by: Edward Small
When: 24/09/2019

Contains unit tests to check the functionality of the `find_25boxes` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import scipy.io as scipy
import numpy as np
from pyowc import core

class Find25BoxesTestCase(unittest.TestCase):
    """
    Test cases for find_25boxes function"
    """

    def setUp(self):
        """
        Setting up variables to be used in test cases
        Is run before each test case
        :return: Nothing
        """
        self.pa_wmo_boxes = scipy.loadmat('../../data/constants/wmo_boxes.mat')
        self.pn_float_long = 57.1794
        self.pn_float_lat = -59.1868

    def test_can_load_matrix(self):
        """
        Check that we can open a matlab matrix file (.mat)
        :return: Nothing
        """
        print("Testing that a matlab matrix can be loaded in...")
        self.assertNotEqual(self.pa_wmo_boxes, '', "Failed to load matlab matrix")

    def test_nearest_neighbour_nan(self):
        """
        Check that the nearest neighbour algorithm returns grid point
        :return: Nothing
        """
        print("Testing nearest_neighbour returns correct grid point")
        x_axis = np.arange(0, 3, 1, int)
        y_axis = np.arange(0, 3, 1, int)
        data = np.full((1, 9), np.arange(3, 12)).reshape(3, 3)
        self.assertEqual(
            core.finders.nearest_neighbour(
                x_axis, y_axis, data, 0, 0
            ), 3, "Nearest neighbour interpolation is wrong"
        )

    def test_nearest_neighbour_returns_point_in_grid(self):
        """
        Check that, even if the nearest neighbour function is given a number outside
        the grid, it still returns a value that exists within the grid
        :return: Nothing
        """
        print("Testing nearest_neighbour returns data from the grid")
        x_axis = np.arange(0, 10000, 5, int)
        y_axis = np.arange(0, 10000, 5, int)
        data = np.full((1, 10000 * 10000), np.arange(1, 10000 * 10000 + 1)).reshape(10000, 10000)
        self.assertIn(
            core.finders.nearest_neighbour(
                x_axis, y_axis, data, 93820495.3958, -9283499184.439
            ), data, "Answer is not in data grid"
        )

    def test_returned_array_is_correct_size(self):
        """
        Check that we get back 25 boxes, each with 4 values
        :return: Nothing
        """
        print("Testing that the returned value is 25x4 matrix")
        pa_wmo_numbers = core.finders.find_25boxes(self.pn_float_long, self.pn_float_lat, self.pa_wmo_boxes)
        self.assertEqual(pa_wmo_numbers.shape, (25, 4), 'returned matrix shape is incorrect')

    def test_returns_array_of_nan_for_nan_input(self):
        """
        Check that even if we get empty inputs, we still receive boxes
        :return: Nothing
        """
        print("Testing nan input")
        pa_wmo_numbers = core.finders.find_25boxes(np.nan, np.nan, self.pa_wmo_boxes)
        self.assertEqual(
            np.allclose(
                pa_wmo_numbers, np.full((25, 4), np.nan), equal_nan=True
            ), True, "nan input doesn't lead to nan output"
        )

    def test_returned_value_exists_in_wmo_boxes(self):
        """
        Check that each of the 25 boxes we receive definitely exists in the original
        collection of wmo boxes
        :return:
        """
        print("Testing returned value exists in wmo_boxes.mat")
        pa_wmo_numbers = core.finders.find_25boxes(self.pn_float_long, self.pn_float_lat, self.pa_wmo_boxes)
        result = np.full((25, 1), False)

        count = 0
        for i in np.array(pa_wmo_numbers):
            for j in np.array(self.pa_wmo_boxes.get('la_wmo_boxes')):
                if i[0] == j[0] and i[1] == j[1] and i[2] == j[2] and i[3] == j[3]:
                    result[count] = True
            count += 1

        self.assertEqual(
            np.all([result, np.full((25, 1), True)]), True, "Some values don't match wmo_boxes.mat"
        )


if __name__ == '__main__':
    unittest.main()
