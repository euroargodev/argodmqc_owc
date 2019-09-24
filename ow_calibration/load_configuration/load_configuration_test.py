"""
-----Load Configuration Test File-----

Written by: Edward Small
When: 18/09/2019

Contains unit tests to check the functionality of the `load_configuration` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
from .load_configuration import load_configuration


class LoadConfigurationTestCase(unittest.TestCase):

    def setUp(self):
        self.lo_system_configuration = load_configuration()

    def test_is_dict(self):
        print("Testing that load_configuration() returns a dictionary...\n")
        self.assertTrue(isinstance(self.lo_system_configuration, dict))

    def test_contains_correct_number_of_keys(self):
        print("Testing length of dictionary...\n")
        self.assertEqual(self.lo_system_configuration.__len__(), 31)

    def test_all_keys_have_value(self):
        print("Testing that the keys in the dictionary all have a value...\n")
        for i in self.lo_system_configuration:
            #print(i, ": ", lo_system_configuration[i]) #uncomment this line to find error location
            self.assertNotEqual(self.lo_system_configuration[i], '')


if __name__ == '__main__':
    unittest.main()
