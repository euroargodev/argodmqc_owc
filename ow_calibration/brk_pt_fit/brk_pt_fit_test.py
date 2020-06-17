"""
-----break point fit Test File-----
Written by: Edward Small
When: 01/06/2020
Contains unit tests to check the functionality of the `brk_pt_fit` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
from ow_calibration.brk_pt_fit.brk_pt_fit import brk_pt_fit


class MyTestCase(unittest.TestCase):
    """
    Test case for 'brk_pt_fit' function
    """

    def setUp(self):
        self.x_obvs = np.arange(-1, 1 + 0.1, 0.1)
        self.x_obvs = np.flip(np.delete(self.x_obvs, 10))
        self.y_obvs = np.array([1.94417105917954,
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

        breaks = [1, 2, 3]
        y_obvs = [1, 2]

        fit_param, residual = brk_pt_fit([1, 2, 3], y_obvs, [1, 2, 3], breaks)

        for i in range(y_obvs.__len__()):
            self.assertEqual(residual[i], y_obvs[i],
                             "residuals should equal independent observations for bad inputs")

        self.assertEqual(fit_param.shape, (breaks.__len__() + 2, 1),
                         "matrix shape should be (b + 2, 1)")

        for i in fit_param:
            self.assertEqual(i, 0, "A should be all zeroes")

    def test_needs_no_weights(self):
        """
        Check that, even without weights, we get an estimate
        :return: Nothing
        """
        print("testing brk_pt_fit runs without weights")
        fit_param_ans = np.array([1.37597459912525,
                                  -0.501251340026727])

        residual_ans = np.array([0.568196460054282,
                                 0.841401212101603,
                                 -1.49501415024474,
                                 0.864152969291804,
                                 0.0211031395509740,
                                 -1.58335338412703,
                                 -1.04047994252411,
                                 -0.235330041510303,
                                 0.996545907177638,
                                 1.01869100647258,
                                 -1.40313535409261,
                                 1.03580374615659,
                                 0.995526245603582,
                                 -0.419847652956731,
                                 0.524866807541146,
                                 -1.45031558324361,
                                 -0.610690751246430,
                                 0.871231976441947,
                                 0.500647389553409,
                                 -4.44089209850063e-16])
        fit_param, residual = brk_pt_fit(self.x_obvs, self.y_obvs, [])

        for i in range(fit_param.__len__()):
            self.assertAlmostEqual(fit_param[i], fit_param_ans[i], 12,
                                   "matlab parameters should match python ones")

        for i in range(residual.__len__()):
            self.assertAlmostEqual(residual[i], residual_ans[i], 12,
                                   "matlab parameters should match python ones")


if __name__ == '__main__':
    unittest.main()
