"""
-----Wrap Longitude Test File-----

Written by: Edward Small
When: 05/12/2019

Contains unit tests to check the functionality of the `wrap_longitude` function.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from pyowc.core.utils import wrap_longitude

class GetDataTestCase(unittest.TestCase):
    """
    Test cases for get_data function
    """

    def test_does_not_wrap(self):
        """
        Shouldn't wrap if longitude values are > 40 and < 320
        :return: Nothing
        """
        print("Testing that wrap_longitude does not wrap values unnecessarily")

        long = np.array([50, 100, 150, 68, 300])
        new_long = wrap_longitude(long)

        for i in range(long.__len__()):
            self.assertTrue(long[i] == new_long[i],
                            "Shouldn't wrap values if they are > 40 and < 320")

    def test_wraps_values(self):
        """
        Should wrap values of the longitudes are around 0 or 360
        :return: Nothing
        """
        print("Testing that wrap_longitude wraps values of they are near 0 or 360")

        long = np.array([-10, 20, 30, 350, 340, 300, 359, 100, 35])
        expected_long = np.array([350, 380, 390, 350, 340, 300, 359, 100, 395])
        new_long = wrap_longitude(long)

        for i in range(long.__len__()):
            self.assertTrue(new_long[i] == expected_long[i],
                            "Should wrap values if they are < 40 and > 320")

    def test_length_equal(self):
        """
        This function shouldn't throw away any spatial values. Output length should equal input
        :return: Nothing
        """
        print("Testing that wrap_longitude returns array of equal length to input")

        long = np.array([-100, 200, 35, 352, 34, 397, 359, 100, 35])
        new_long = wrap_longitude(long)

        self.assertTrue(long.__len__() == new_long.__len__(),
                        "input and output should be equal in length")


if __name__ == '__main__':
    unittest.main()
