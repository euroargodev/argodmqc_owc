""" Tests for core.stats module functions """
import os
import unittest
import numpy as np
from scipy.io import loadmat

from pyowc import core
from pyowc.core.stats import fit_cond, noise_variance, signal_variance
from . import TESTS_CONFIG


class BrkPtFit(unittest.TestCase):
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

        fit_param, residual = core.stats.brk_pt_fit([1, 2, 3], y_obvs, [1, 2, 3], breaks)

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
        fit_param, residual = core.stats.brk_pt_fit(self.x_obvs, self.y_obvs, [])

        for i in range(fit_param.__len__()):
            self.assertAlmostEqual(fit_param[i], fit_param_ans[i], 12,
                                   "matlab parameters should match python ones")

        for i in range(residual.__len__()):
            self.assertAlmostEqual(residual[i], residual_ans[i], 12,
                                   "matlab parameters should match python ones")

    def test_runs_with_weights(self):
        """
        Check we get the right answer when running with weights
        :return: nothing
        """
        print("Testing brk_pt_ft runs using weighted values")

        fit_param_ans = np.array([1.20259478507721,
                                  -0.587941247050749])
        residual_ans = np.array([0.741576274102326,
                                 1.01478102614965,
                                 -1.32163433619669,
                                 1.03753278333985,
                                 0.194482953599018,
                                 -1.40997357007898,
                                 -0.867100128476065,
                                 -0.0619502274622590,
                                 1.16992572122568,
                                 1.19207082052062,
                                 -1.22975554004457,
                                 1.20918356020464,
                                 1.16890605965163,
                                 -0.246467838908687,
                                 0.698246621589190,
                                 -1.27693576919556,
                                 -0.437310937198385,
                                 1.04461179048999,
                                 0.674027203601453,
                                 0])

        fit_param, residual = core.stats.brk_pt_fit(self.x_obvs, self.y_obvs, self.w_i)

        for i in range(fit_param.__len__()):
            self.assertAlmostEqual(fit_param[i], fit_param_ans[i], 12,
                                   "matlab parameters should match python ones")

        for i in range(residual.__len__()):
            self.assertAlmostEqual(residual[i], residual_ans[i], 12,
                                   "matlab parameters should match python ones")

    def test_det_zero(self):
        """
        Check that we still get an answer if determinant is o
        :return: nothing
        """
        print("Testing brk_pt_fit returns if det(A) == 0")

        y_obvs = np.array([0, 0, 0])

        fit_param, residual = core.stats.brk_pt_fit(np.array([0, 0, 0]), y_obvs, [])

        for i in range(fit_param.__len__()):
            self.assertEqual(fit_param[i], 0, "should be all zeroes when det(a) == 0")

        for i in range(residual.__len__()):
            self.assertEqual(residual[i], y_obvs[i], "should be all zeroes when det(a) == 0")


class BuildCov(unittest.TestCase):
    """
    Test cases for build_cov function
    """

    def setUp(self):
        """
        Set up for test
        :return: Nothing
        """

        self.config = {"MAPSCALE_LONGITUDE_SMALL": 4,
                       "MAPSCALE_LATITUDE_SMALL": 2,
                       "MAPSCALE_PHI_SMALL": 0.1,
                       "MAP_USE_PV": 0}
        self.ptmp = np.array([[0.7058, 0.7039, 0.8285],
                              [0.6713, 0.6664, 0.7432],
                              [0.8257, 0.8224, 0.7804],
                              [0.7452, 0.7411, 1.1980],
                              [0.7836, 0.7802, 1.1504],
                              [1.2008, 1.2010, 1.2497],
                              [1.1496, 1.1481, 1.3036],
                              [1.2520, 1.2553, 1.0921],
                              [1.3039, 1.3046, np.nan],
                              [1.0947, 1.0962, np.nan]])
        self.coord_float = np.array([[0.0572, -0.0592, 5.1083],
                                     [0.0578, -0.0591, 5.0993],
                                     [0.0586, -0.0585, 5.0861]]) * 1.0e+03

    def test_returns_numpy_array(self):
        """
        Check that build_cov returns a numpy array
        :return: Nothing
        """

        print("Testing that build_cov returns a numpy array")

        test = core.stats.build_cov(self.ptmp, self.coord_float, self.config)

        self.assertEqual(type(test), np.ndarray, "build_cov did not return a numpy array")

    def test_returns_correct_size(self):
        """
        Check that build_cov returns a matrix that is the correct size
        :return: Nothing
        """

        print("Testing that build_cov returns correct shape matrix")

        test = core.stats.build_cov(self.ptmp, self.coord_float, self.config)

        self.assertEqual(test.shape,
                         (self.coord_float.shape[0] * self.ptmp.shape[0],
                          self.coord_float.shape[0] * self.ptmp.shape[0]),
                         "build_cov returned a matrix of incorrect size")

    def test_returns_correct_elements(self):
        """
        Check that build_cov returns a matrix that is the correct size
        :return: Nothing
        """

        print("Testing that build_cov returns correct shape matrix")

        # expected = loadmat("../../data/test_data/build_cov/cov.mat")['test_cov_1']
        expected_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "build_cov", "cov.mat"])
        expected = loadmat(expected_path)['test_cov_1']

        expected_size = expected.shape

        test = core.stats.build_cov(self.ptmp, self.coord_float, self.config)

        for i in range(0, expected_size[0]):
            for j in range(0, expected_size[1]):
                self.assertAlmostEqual(test[i, j], expected[i, j], 15,
                                       "covariance matrix is incorrect")


class CalcPiecewiseFit(unittest.TestCase):
    """
    Test cases for 'calc_piecewisefit' function
    """

    def test_custom(self):
        """
        Change variables in this test to use different mapped outputs
        :return: nothing
        """
        core.stats.calc_piecewisefit("/", TESTS_CONFIG['TEST_FLOAT_SOURCE'], TESTS_CONFIG)

        test_data_path = os.path.sep.join([TESTS_CONFIG['FLOAT_CALIB_DIRECTORY'],
                                      TESTS_CONFIG['FLOAT_CALIB_PREFIX'] +
                                      TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                      TESTS_CONFIG['FLOAT_CALIB_POSTFIX']])
        test = loadmat(test_data_path)

        matlab_data_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], 'float_calib_test',
                                             TESTS_CONFIG['FLOAT_CALIB_PREFIX'] +
                                             TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                             TESTS_CONFIG['FLOAT_CALIB_POSTFIX']])
        matlab = loadmat(matlab_data_path)

        python_sal = test['cal_SAL']
        matlab_sal = matlab['cal_SAL']

        self.assertEqual(python_sal.shape, matlab_sal.shape)

        for i in range(python_sal.shape[0]):
            for j in range(python_sal.shape[1]):
                if ~np.isnan(python_sal[i, j]):
                    self.assertAlmostEqual(python_sal[i, j], matlab_sal[i, j], 3)

        python_sal_err = test['cal_SAL_err']
        matlab_sal_err = matlab['cal_SAL_err']

        for i in range(python_sal_err.shape[0]):
            for j in range(python_sal_err.shape[1]):
                if ~np.isnan(python_sal_err[i, j]):
                    self.assertAlmostEqual(python_sal_err[i, j], matlab_sal_err[i, j], 3)


class Covarxytpv(unittest.TestCase):
    """
    Test cases for covar_xyt_pv function
    """

    def setUp(self):
        """
        set up variables to use for testing
        :return: Nothing
        """
        self.points1 = np.array([[-0.057996, 0.053195, 1.9740875, 5.229838],
                                 [-0.0564902, 0.0631170, 1.9870367, 4.6300392],
                                 [-0.05208, 0.0619770, 1.9941118, 4.6536932]]) * 10 ** 3
        self.points2 = self.points1
        self.lat = 4
        self.long = 8
        self.age = 20
        self.phi = 0.5
        self.p_v = 0

    def test_returns_array(self):
        """
        Check that the function returns an array if given an array
        :return: Nothing
        """
        print("Testing that covar_xyt_pv returns an array")
        covar = core.stats.covar_xyt_pv(self.points1, self.points2, self.lat, self.long,
                             self.age, self.phi, self.p_v)

        self.assertTrue(isinstance(covar, np.ndarray), "Covariance isn't numpy array")

    def test_returns_ones_for_same_value(self):
        """
        Check that entering the same value for points 1 and 2 gives a matrix of ones
        :return: Nothing
        """
        print("Testing that covar_xyt_pv return 1's if input points are identical")

        covar = core.stats.covar_xyt_pv(np.array([[1, 2, 3, 4], [1, 2, 3, 4]]),
                             np.array([[1, 2, 3, 4], [1, 2, 3, 4]]),
                             self.lat, self.long, self.age, self.phi, self.p_v)
        expected = np.array([[1, 1], [1, 1]])

        for i in range(0, covar.__len__()):
            for j in range(0, covar[i].__len__()):
                self.assertEqual(covar[i][j], expected[i][j], "covariance isn't just 1's")

    def test_returns_zeroes_for_extremely_different_value(self):
        """
        Check that entering extreme value for points 1 and 2 gives a matrix of zeroes
        :return: Nothing
        """
        print("Testing that covar_xyt_pv return 0's if input points are extremely different")

        covar = core.stats.covar_xyt_pv(np.array([[1, 2, 3, 4], [1000, 2000, 3000, 4000]]),
                             np.array([[-1000, -2000, -3000, 4000],
                                       [-1 * 9999, 2 * 9999, -3 * 9999, -4 * 999]]),
                             self.lat, self.long, self.age, self.phi, self.p_v)
        expected = np.array([[0, 0], [0, 0]])

        for i in range(0, covar.__len__()):
            for j in range(0, covar[i].__len__()):
                self.assertAlmostEqual(covar[i][j] * 1000000, expected[i][j] * 1000000,
                                       "covariance isn't just 1's")

    def test_return_matrix_shape_correct_multidimensional(self):
        """
        Check that, if given an m*4 points 1 matrix and a n n*5 points 2 matrix
        the function returns an m*n matrix (m > 1)
        :return: Nothing
        """
        print("Testing that covar_xyt_pv returns matrix of the correct shape for multidemensional")

        covar = core.stats.covar_xyt_pv(self.points1, self.points2, self.lat, self.long,
                             self.age, self.phi, self.p_v)

        expected = (3, 3)

        for i in range(0, covar.shape.__len__()):
            self.assertEqual(covar.shape[i], expected[i], "covariance matrix is wrong shape")

    def test_return_matrix_shape_correct_one_demensional(self):
        """
        Check that, if given an 1*4 points 1 matrix and a n n*5 points 2 matrix
        the function returns an 1*n matrix (m > 1)
        :return: Nothing
        """
        print("Testing that covar_xyt_pv returns matrix of the correct shape one dimensional")

        covar = core.stats.covar_xyt_pv(np.array([-0.057996, 0.053195, 1.9740875, 5.229838]) * 10 ** 3,
                             self.points2, self.lat, self.long, self.age, self.phi, self.p_v)
        expected_shape = (1, 3)

        self.assertEqual(covar.shape, expected_shape, "covariance matrix is wrong shape")

    def test_returns_expected_answers(self):
        """
        Check that we get the answers we expect
        :return: Nothing
        """
        print("Testing that covar_xyt_pv returns the expected result")

        expected = np.array([[1, 0.001350185414046, 0.001712995003897],
                             [0.001350185414046, 1, 0.600332430927527],
                             [0.001712995003897, 0.600332430927527, 1]])
        covar = core.stats.covar_xyt_pv(self.points1, self.points2, self.lat, self.long,
                             self.age, self.phi, self.p_v)

        for i in range(0, covar.__len__()):
            for j in range(0, covar[i].__len__()):
                self.assertAlmostEqual(covar[i][j], expected[i][j], 15,
                                       "covariances is not as expected")

    def test_allows_1_dimensional_data(self):
        """
        Check that using 1 dimensional data sets does not through an error
        :return: Nothing
        """
        print("Testing that covar_xyt_pv can use one dimensional data")

        covar = core.stats.covar_xyt_pv(np.array([-0.057996, 0.053195, 1.9740875, 5.229838]) * 10 ** 3,
                             self.points2, self.lat, self.long, self.age, self.phi, self.p_v)

        expected_ans = [1, 0.122561894349782, 0.012339774448825]

        for i in range(covar.__len__()):
            self.assertAlmostEqual(covar[0][i], expected_ans[i], 15,
                                   "unexpected answer with one dimensional data")


class Covarxypv(unittest.TestCase):
    """
    Test cases for covarxy_pv function
    """

    def setUp(self):
        """
        Set up for test
        :return: Nothing
        """

        self.input_coords = np.array([0.0572, -0.0592, 5.1083]) * 1.0e+03
        self.coords = np.array([[0.0572, -0.0592, 5.1083],
                                [0.0578, -0.0591, 5.0993],
                                [0.0586, -0.0585, 5.0861]]) * 1.0e+03
        self.long = 4
        self.lat = 2
        self.phi = 0.1
        self.no_pv = 0
        self.yes_pv = 1

    def test_return_shape(self):
        """
        Check that the returned covariance matrix is the correct size
        :return: Nothing
        """
        print("Testing that covarxy_pv returns 1 dimensional matrix")

        cov = core.stats.covarxy_pv(self.input_coords, self.coords, self.long, self.lat, self.phi, self.no_pv)

        self.assertEqual(cov.shape, (3,), "covarxy_pv matrix is incorrect shape")

    def test_returns_correct(self):
        """
        Check that the returned covriance matrix contains the correct values
        :return: Nothing
        """

        print("Testing that covarxy_pv returns correct values")

        ans = np.array([1, 0.9134, 0.5941])

        cov = core.stats.covarxy_pv(self.input_coords, self.coords, self.long, self.lat, self.phi, self.no_pv)

        for i in range(0, ans.size):
            self.assertAlmostEqual(ans[i], cov[i],
                                   4, "elements in covarxy_pv matrix are not correct")

    def test_returns_correct_pv(self):
        """
        Check that the returned covriance matrix contains the correct values with
        potential vorticity
        :return: Nothing
        """

        print("Testing that covarxy_pv returns correct values with potential vorticity")

        ans = np.array([1, 0.9101, 0.5828])

        cov = core.stats.covarxy_pv(self.input_coords, self.coords, self.long, self.lat,
                         self.phi, self.yes_pv)

        for i in range(0, ans.size):
            self.assertAlmostEqual(ans[i], cov[i],
                                   4, "elements in covarxy_pv matrix are not correct")

    def test_returns_ones(self):
        """
        Check that we get a matrix of almost ones if data is very close (according to scale)
        :return: nothing
        """

        print("Testing that covarxy_pv returns almost ones if data is close")

        cov = core.stats.covarxy_pv(self.input_coords, self.coords, 99999999, 99999999, self.phi, self.no_pv)
        for i in cov:
            self.assertAlmostEqual(i, 1, 15, "elements in covarxy_pv matrix are not correct")


#pylint: disable=too-many-instance-attributes
class FitCond(unittest.TestCase):
    """
    Test cases for 'fit_cond' function
    """
    def setUp(self):
        # fit_input = loadmat("../../data/test_data/fit_cond/fit_cond_input.mat")
        # fit_out = loadmat("../../data/test_data/fit_cond/fit_cond_output.mat")
        fit_input = loadmat(os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "fit_cond", "fit_cond_input.mat"]))
        fit_out = loadmat(os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "fit_cond", "fit_cond_output.mat"]))

        self.in_x = fit_input['x']
        self.in_y = fit_input['y']
        self.in_err = fit_input['n_err']
        self.in_cov = fit_input['lvcov']

        self.xfit = fit_out['xfit']
        self.condslope = fit_out['condslope']
        self.condslope_err = fit_out['condslope_err']
        self.time_derivself = fit_out['time_deriv']
        self.time_deriv_err = fit_out['time_deriv_err']
        self.sta_mean = fit_out['sta_mean']
        self.sta_rms = fit_out['sta_rms']
        self.ndf = fit_out['NDF']
        self.fit_coef = fit_out['fit_coef']
        self.fit_breaks = fit_out['fit_breaks']

    def test_max_no_breaks(self):
        """
        Check return values if we specify a maximum number of breaks
        :return: nothing
        """
        print("Testing fit_cond for max_no_breaks")

        python_test = fit_cond(self.in_x, self.in_y, self.in_err,
                               self.in_cov, 'max_no_breaks', 4)

        for i in range(python_test[0].__len__()):
            self.assertEqual(python_test[0][i], self.xfit[i],
                             "xfit is incorrect")

        for i in range(python_test[2][0].__len__()):
            self.assertAlmostEqual(python_test[2][0, i], 0.000122830414419824, 12,
                                   "slope error is incorrect")

        for i in range(python_test[5].shape[0]):
            for j in range(python_test[5].shape[1]):
                if ~np.isnan(python_test[5][i, j]):
                    self.assertAlmostEqual(python_test[5][i, j], self.sta_mean[i, j], 12,
                                           "mean is incorrect")

        for i in range(python_test[6].shape[0]):
            for j in range(python_test[6].shape[1]):
                if ~np.isnan(python_test[6][i, j]):
                    self.assertAlmostEqual(python_test[6][i, j], self.sta_rms[i, j], 12,
                                           "rms is incorrect")

        self.assertEqual(python_test[7], self.ndf, "degrees of freedom is incorrect")

    def test_fixed_breaks(self):
        """
        Check that we can run this function with set break points
        :return: nothing
        """
        print("Testing that fit_cond returns values when using fixed breaks")

        python_test = fit_cond(self.in_x, self.in_y, self.in_err,
                               self.in_cov, 'breaks', np.array([0.3, 0.7]), 'max_no_breaks', 4)

        self.assertEqual(python_test.__len__(), 10, "should return 10 outputs")


class NoiseVariance(unittest.TestCase):
    """
    Test cases for noise_variance function
    """

    def setUp(self):
        self.sal = np.array([34.4988, 34.3267, 34.0346])
        self.lat = np.array([-57.9960, -56.4902, -52.0800])
        self.long = np.array([53.1950, 63.1170, 61.9770])

    def test_returns_float(self):
        """
        Check that noise_variance returns a float
        :return: Nothing
        """
        print("Testing that noise_variance returns a float")

        noise_var1 = noise_variance(self.sal, self.lat, self.long)
        noise_var2 = noise_variance(np.array([0, 2]), np.array([1, -1]), np.array([-1, 1]))

        self.assertTrue(isinstance(noise_var1, float), "noise variance is not a float")
        self.assertTrue(isinstance(noise_var2, float), "noise variance is not a float")

    def test_returns_0_if_no_unique_points(self):
        """
        Check that noise_variance returns 0 if it cannot find any unique points
        :return: Nothing
        """
        print("Testing that noise_variance returns 0 for no unique points")

        noise_var = noise_variance(np.array([0, 0]), np.array([-1, -1]), np.array([1, 1]))
        self.assertEqual(noise_var, 0, "Variance is not 0 for equal points")

    def test_returns_expected(self):
        """
        Check that noise_variance returns the expected answer
        :return: Nothing
        """
        print("Testing that noise_variance returns the expected value")

        expected = 0.033377205000001
        noise_var = noise_variance(self.sal, self.lat, self.long)

        self.assertAlmostEqual(noise_var, expected, 16, "Did not receive expected answer")


class SignalVariance(unittest.TestCase):
    """
    Test cases for signal_variance function
    """

    def test_returns_float(self):
        """
        Check that we return a float if given some data
        :return: Nothing
        """
        print("Testing that signal_variance returns a float")

        var = signal_variance([1, 2, 3, 4, 5])
        self.assertTrue(isinstance(var, float), "signal variance is not a float")

    def test_throws_exception(self):
        """
        Check that we thrown an exception if no valid salinities are given
        :return: Nothing
        """
        print("Testing that signal_variance throws an exception for no valid salinities")

        with self.assertRaises(Exception) as no_valid_sal:
            signal_variance([0, 0, float('nan'), float('nan')])

        self.assertTrue('Received no valid salinity values when calculating signal variance'
                        in str(no_valid_sal.exception))

    def test_nans_are_ignored(self):
        """
        Check that we ignore 0's and nan values
        :return: Nothing
        """
        print("Testing that signal_variance ignores 0's/NaNs correctly")

        expected = signal_variance([1, 2, 3, 4, 5])
        zeroes = signal_variance([1, 0, 2, 0, 3, 0, 4, 0, 5, 0])
        nans = signal_variance([1, 2, 3, 4, 5, float('nan'), float('nan')])

        self.assertEqual(expected, zeroes, "signal variance is not ignoring 0's")
        self.assertEqual(expected, nans, "sign_variance is not ignoring NaN's")

    def test_negtive_inputs_against_positive(self):
        """
        Check that giving a set of negative inputs gives the same result as positive inputs
        :return: Nothing
        """
        print("Testing that signal_variance returns the same result for input and -input")

        positive = signal_variance([30, 20, 10, 40, 50])
        negative = signal_variance([-30, -20, -10, -40, -50])

        self.assertEqual(positive, negative,
                         "signal_variance is not returning the same result for input and -input")

    def test_returns_correct_result(self):
        """
        Check that we get the expected result from some inputs
        :return: Nothing
        """
        print("Testing that signal_variance gives the expected results")

        expected = 0.185000000000001
        ans = signal_variance([-35, -35.5, -35.7, -36.2])

        self.assertEqual(ans, expected, "signal_variance did not give the expected result")


if __name__ == '__main__':
    unittest.main()
