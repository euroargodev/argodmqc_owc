import unittest
from HelloWorld import hello_world_string


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(hello_world_string(), "Hello World")


if __name__ == '__main__':
    unittest.main()
