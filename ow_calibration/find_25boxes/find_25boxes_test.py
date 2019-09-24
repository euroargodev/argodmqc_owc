"""
-----Load Configuration Test File-----

Written by: Edward Small
When: 24/09/2019

Contains unit tests to check the functionaility of the `find_25boxes` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import scipy.io as scipy

class Find25BoxesTestCase(unittest.TestCase):

    def setUp(self):
        self.pa_wmo_boxes = scipy.loadmat('../../data/constants/wmo_boxes.mat')
        self.pn_float_long = 57.1794
        self.pn_float_lat = -59.1868

    def test_can_load_matrix(self):
        print("Testing that a matlab matrix can be loaded in...")
        self.assertNotEqual(self.pa_wmo_boxes, '', 'Failed to load matlab matrix')

    def test_returned_array_is_correct_size(self):
        print("Testing that the returned value is 25x4 matrix")
        pa_wmo_numbers = find_25boxes(self.pn_float_long, self.pn_float_lat, self.pa_wmo_boxes)
        self.assertEqual(pa_wmo_numbers.shape(), (25, 4, 0), 'returned matrix shape is incorrect')


if __name__ == '__main__':
    unittest.main()
