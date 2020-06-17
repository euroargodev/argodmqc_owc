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

    def setUp(self):
        self.x = np.arange(-1, 1 + 0.1, 0.1)
        self.x = np.flip(np.delete(self.x, 10))
        self.y = np.array([1.94417105917954,
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

        self.w_i = np.zeros((20, 20))
        no_weights = self.w_i.shape[0]

        for i in range(no_weights):
            self.w_i[i, i] = weights[i]

    def test_bad_inputs_gives_outputs(self):
        """
        Check that, if observations are the wrong shape, we get a bad output
        :return: Nothing
        """
        print("Testing brk_pt_fit returns values for bad arrays")

        b = [1, 2, 3]

        a, residual = brk_pt_fit([1, 2, 3], [1, 2], [1, 2, 3], b)

        self.assertEqual(residual, 999, "residual should be 999 if bad inputs occur")
        self.assertEqual(a.shape, (b.__len__() + 2, 1), "matrix shape should be (b + 2, 1)")

        for i in a:
            self.assertEqual(i, 0, "A should be all zeroes")

    def test_stuff(self):
        """

        :return:
        """
        brk_pt_fit(self.x, self.y, self.w_i)

if __name__ == '__main__':
    unittest.main()
