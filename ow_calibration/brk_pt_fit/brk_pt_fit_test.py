"""
-----break point fit Test File-----
Written by: Edward Small
When: 01/06/2020
Contains unit tests to check the functionality of the `brk_pt_fit` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

from ow_calibration.brk_pt_fit.brk_pt_fit import brk_pt_fit
import numpy as np
import unittest


class MyTestCase(unittest.TestCase):
    """
    Test case for 'brk_pt_fit' function
    """
    def test_something(self):
        x = np.arange(-1, 1 + 0.1, 0.1)
        x = np.flip(np.delete(x, 10))
        y = np.array([1.94417105917954,
                      2.21737581122686,
                      -0.119039551119482,
                      2.24012756841706,
                      1.39707773867623,
                      -0.207378785001771,
                      0.335494656601145,
                      1.14064455761495,
                      2.37252050630289,
                      2.39466560559783,
                      -0.0271607549673552,
                      2.41177834528185,
                      2.37150084472884,
                      0.956126946168524,
                      1.90084140666640,
                      -0.0743409841183540,
                      0.765283847878825,
                      2.24720657556720,
                      1.87662198867866,
                      2.37847727917871])

        weights = np.array([1.48361104873488,
                            0.553567517861284,
                            1.77369395880317,
                            1.90098987163633,
                            1.51810273228666,
                            1.63661019586750,
                            1.61469870218737,
                            1.08834052930125,
                            1.48321683526633,
                            0.756780031717343,
                            1.55906913202941,
                            0.547749269566131,
                            0.915384477441335,
                            0.569257085946731,
                            0.645697671853771,
                            1.73518674249094,
                            1.54224293446373,
                            0.975649220091291,
                            1.92533307325753,
                            0.551669120754363])

        w_i = np.zeros((20, 20))
        no_weights = w_i.shape[0]

        for i in range(no_weights):
            w_i[i, i,] = weights[i]

        brk_pt_fit(x, y, w_i)


if __name__ == '__main__':
    unittest.main()
