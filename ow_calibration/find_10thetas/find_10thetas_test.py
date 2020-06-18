"""
-----Find 10 Theta Levels Test File-----

Written by: Edward Small
When: 01/06/2020

Contains unit tests to check the functionality of the `find_10thetas` function

Changes can be made to this file to test all different kinds of floats, as long
as you have data for comparison

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import unittest
import numpy as np
import scipy.io as scipy
from ow_calibration.find_10thetas.find_10thetas import find_10thetas


# pylint: disable=fixme
# pylint: disable=too-many-instance-attributes
class Find10ThetasTestCase(unittest.TestCase):
    """
    Test cases for find_10thetas function
    """

    def setUp(self):
        # set up data to run find_10thetas
        float_source = "3901960"
        mapped_values_path = "data/test_data/float_mapped_test/map_" + float_source + ".mat"
        source_values_path = "data/float_source/" + float_source + ".mat"

        mapped_values = scipy.loadmat(mapped_values_path)
        source_values = scipy.loadmat(source_values_path)

        self.sal = source_values['SAL']
        self.ptmp = source_values['PTMP']
        self.pres = source_values['PRES']

        self.la_ptmp = mapped_values['la_ptmp']

        # matlab answers

        self.tlevels = np.array([1.30830169059140,
                                 1.36299325725812,
                                 1.25487242273343,
                                 1.41236755027244,
                                 1.52565232734614,
                                 1.20028176089126,
                                 0.823372822547086,
                                 1.57739641634854,
                                 1.67653434735605,
                                 1.14931045824918]).reshape(-1, 1)

        self.plevels = np.array([1405.10000610352,
                                 1355.10000610352,
                                 1455.10000610352,
                                 1305.10000610352,
                                 1205.10000610352,
                                 1505.10000610352,
                                 1905.10000610352,
                                 1155.10000610352,
                                 1055.10000610352,
                                 1555.10000610352]).reshape(-1, 1)

        self.var_s_thetalevels = np.array([np.nan,
                                           0.00147259240080181,
                                           0.00690521308242283,
                                           0.00599714167308804,
                                           0.00598545268753434,
                                           0.00802480124455128,
                                           0.00463859096757298,
                                           0.00686978997540036,
                                           0.00190670919307289,
                                           0.00717484312550677,
                                           0.00745579840256897,
                                           0.00154219768993262,
                                           0.00475705221676092,
                                           0.00107598836150616,
                                           0.00167595421613088,
                                           0.00201337210437592,
                                           0.00604797739660125,
                                           2.21183340374061e-05,
                                           0.00546642420107933,
                                           1.58194626443039e-05,
                                           8.94793528356452e-06,
                                           0.00762037347540679,
                                           7.90075020336859e-06,
                                           5.47403808324673e-06,
                                           4.92215091502965e-06,
                                           6.55854597121638e-06,
                                           1.23629324722680e-05,
                                           5.46616782275005e-05,
                                           9.31582204233218e-05,
                                           9.74824249506302e-05,
                                           7.63924400054304e-05,
                                           5.73756473973588e-05,
                                           0.0126875498106624,
                                           0.0128534629746730,
                                           1.54955553212721e-05,
                                           np.nan,
                                           np.nan]).reshape(-1, 1)

        self.thetalevels = np.array([np.nan,
                                     1.93961370997963,
                                     1.99500305452812,
                                     2.02025501639417,
                                     2.02694204816290,
                                     2.00848644455079,
                                     1.99059381638721,
                                     1.97075994671766,
                                     1.94046894653132,
                                     1.92705528342800,
                                     1.89375691941123,
                                     1.86741374137747,
                                     1.83113346775466,
                                     1.79263813584680,
                                     1.74856506989514,
                                     1.71174508509146,
                                     1.68506761991558,
                                     1.67653434735605,
                                     1.62733096351643,
                                     1.57739641634854,
                                     1.52565232734614,
                                     1.47091410199299,
                                     1.41236755027244,
                                     1.36299325725812,
                                     1.30830169059140,
                                     1.25487242273343,
                                     1.20028176089126,
                                     1.14931045824918,
                                     1.09558301978575,
                                     1.04691796354743,
                                     1.00240760069712,
                                     0.954500067130170,
                                     0.910028340866476,
                                     0.866025309391749,
                                     0.823372822547086,
                                     0.781478873355349,
                                     0.741563516899611]).reshape(-1, 1)

    def test_theta_levels_shape(self):
        """
        Check that we get 10 levels
        :return: Nothing
        """
        print("Testing that we only get 10 theta levels")
        t_levels, p_levels, index, var_sal_theta, theta_levels = find_10thetas(self.sal,
                                                                               self.ptmp,
                                                                               self.pres,
                                                                               self.la_ptmp,
                                                                               0, 0, 0, 0, 0.5)

        self.assertEqual(t_levels.shape, self.tlevels.shape,
                         "Got incorrect number of theta levels")
        self.assertEqual(p_levels.shape, self.plevels.shape,
                         "Got incorrect number of pressure levels")
        self.assertEqual(index.shape, (10, self.sal.shape[1]),
                         "have incorrect number of indices")
        self.assertEqual(var_sal_theta.shape, self.var_s_thetalevels.shape,
                         "Got incorrect number of theta levels")
        self.assertEqual(theta_levels.shape, self.thetalevels.shape,
                         "Got incorrect number of theta levels")

    def test_theta_levels_values(self):
        """
        Check that we get 10 levels
        :return: Nothing
        """
        print("Testing that we only get 10 theta levels")
        t_levels, p_levels, index, var_sal_theta, theta_levels = find_10thetas(self.sal,
                                                                               self.ptmp,
                                                                               self.pres,
                                                                               self.la_ptmp,
                                                                               0, 0, 0, 0, 0.5)

        for i in range(t_levels.__len__()):
            self.assertAlmostEqual(t_levels[i, 0], self.tlevels[i, 0], 10,
                                   "Got incorrect theta level")
            self.assertAlmostEqual(p_levels[i, 0], self.plevels[i, 0], 10,
                                   "Got incorrect pressure level")

        for i in range(self.var_s_thetalevels.shape[0]):
            for j in range(self.var_s_thetalevels.shape[1]):
                if not np.isnan(var_sal_theta[i, j]):
                    self.assertAlmostEqual(var_sal_theta[i, j], self.var_s_thetalevels[i, j], 10,
                                           "salinity variance is incorrect")

        for i in range(self.thetalevels.shape[0]):
            for j in range(self.thetalevels.shape[1]):
                if not np.isnan(theta_levels[i, j]):
                    self.assertAlmostEqual(theta_levels[i, j], self.thetalevels[i, j], 10,
                                           "salinity variance is incorrect")

        self.assertGreater(index.__len__(), 0, "Should have some indices")

        # TODO: add further tests for bounded values


if __name__ == '__main__':
    unittest.main()
