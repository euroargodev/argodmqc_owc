import unittest
from load_configuration import load_configuration


class LoSystemConfigurationTestCase(unittest.TestCase):

    def test_is_dict(self):
        print("Testing that load_configuration() returns a dictionary...\n")
        lo_system_configuration = load_configuration()
        self.assertTrue(isinstance(lo_system_configuration, dict))

    def test_all_keys_have_value(self):
        print("Testing that the keys in the dictionary all have a value...\n")
        lo_system_configuration = load_configuration()
        for i in lo_system_configuration:
            #print(i, ": ", lo_system_configuration[i]) #uncomment this line to find error location
            self.assertNotEqual(lo_system_configuration[i], '')


if __name__ == '__main__':
    unittest.main()
