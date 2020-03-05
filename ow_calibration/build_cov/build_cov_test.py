import unittest
from ow_calibration.build_cov.build_cov import build_cov

class MyTestCase(unittest.TestCase):
    def test_something(self):
        build_cov(1, 2, 3)


if __name__ == '__main__':
    unittest.main()
