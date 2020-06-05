import unittest
import numpy as np
from ow_calibration.sorter.sorter import sorter

class MyTestCase(unittest.TestCase):
    def test_something(self):

        A = np.array([-1, 0])
        B = np.arange(-1, 1.01, 0.01)

        sorter(A, B)



if __name__ == '__main__':
    unittest.main()
