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
from .find_25boxes import nearest_neighbour, find_25boxes

class Find25BoxesTestCase(unittest.TestCase):

    def setUp(self):
        self.pa_wmo_boxes = scipy.loadmat('../data/constants/wmo_boxes.mat')
        self.pn_float_long = 57.1794
        self.pn_float_lat = -59.1868

    def test_can_load_matrix(self):
        print("Testing that a matlab matrix can be loaded in...")
        self.assertNotEqual(self.pa_wmo_boxes, '', "Failed to load matlab matrix")

    def test_nearest_neighbour_nan(self):
        print("Testing nearest_neighbour returns correct grid point")
        x = np.arange(0,3,1, int)
        y = np.arange(0, 3, 1, int)
        data = np.full((1,9), np.arange(3, 12)).reshape(3, 3)
        self.assertEqual(nearest_neighbour(x,y,data,0,0), 3, "Nearest neighbour interpolation is wrong")

    def test_nearest_neighbour_returns_point_in_grid(self):
        print("Testing nearest_neighbour returns data from the grid")
        x = np.arange(0, 10000, 5, int)
        y = np.arange(0, 10000, 5, int)
        data = np.full((1, 10000*10000), np.arange(1, 10000*10000+1)).reshape(10000, 10000)
        self.assertIn(nearest_neighbour(x,y,data, 93820495.3958, -9283499184.439), data, "Answer is not in data grid")

    def test_returned_array_is_correct_size(self):
        print("Testing that the returned value is 25x4 matrix")
        pa_wmo_numbers = find_25boxes(self.pn_float_long, self.pn_float_lat, self.pa_wmo_boxes)
        self.assertEqual(pa_wmo_numbers.shape, (25, 4), 'returned matrix shape is incorrect')

    def test_returns_array_of_nan_for_nan_input(self):
        print("Testing nan input")
        pa_wmo_numbers = find_25boxes(np.nan, np.nan, self.pa_wmo_boxes)
        self.assertEqual(np.allclose(pa_wmo_numbers, np.full((25, 4), np.nan), equal_nan=True), True,
                         "nan input doesn't lead to nan output")

    def test_returned_value_exists_in_wmo_boxes(self):
        print("Testing returned value exists in wmo_boxes.mat")
        pa_wmo_numbers = find_25boxes(self.pn_float_long, self.pn_float_lat, self.pa_wmo_boxes)
        result = np.full((25,1), False)

        count = 0
        for i in np.array(pa_wmo_numbers):
            for j in np.array(self.pa_wmo_boxes.get('la_wmo_boxes')):
                if i[0] == j[0] and i[1] == j[1] and i[2] == j[2] and i[3] == j[3]:
                    result[count] = True
            count += 1

        self.assertEqual(np.all([result, np.full((25, 1), True)]), True, "Some values don't match wmo_boxes.mat")


if __name__ == '__main__':
    unittest.main()
