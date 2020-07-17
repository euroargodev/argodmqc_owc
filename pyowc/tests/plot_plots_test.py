import os
import unittest
import copy
from unittest.mock import patch
from scipy.io import loadmat
import numpy as np

from pyowc.core.finders import find_10thetas
from pyowc.plot import plots, utils
from . import TESTS_CONFIG


# pylint: disable=bare-except
# pylint: disable=unused-argument
class CalSalCurve(unittest.TestCase):
    """
    Test cases for cal_sal_curve_plot function
    """

    @patch("pyowc.plot.plots.plt.show")
    def test_plot_runs(self, mockshow):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that calibrated salinity curve plot throws no errors")

        # get the data we need
        cal_data_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "float_calib_test",
                                          TESTS_CONFIG['FLOAT_CALIB_PREFIX'] +
                                          TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                          TESTS_CONFIG['FLOAT_CALIB_POSTFIX']])
        cal_data = loadmat(cal_data_path)
        float_data_path = os.path.sep.join([TESTS_CONFIG['FLOAT_SOURCE_DIRECTORY'],
                                            TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                            TESTS_CONFIG['FLOAT_SOURCE_POSTFIX']])
        float_data = loadmat(float_data_path)

        sal = float_data['SAL']
        cal_sal = cal_data['cal_SAL']
        cal_sal_err = cal_data['cal_SAL_err']
        sta_sal = cal_data['sta_SAL']
        sta_sal_err = cal_data['sta_SAL_err']
        sta_mean = cal_data['sta_mean']
        pcond_factor = np.array(cal_data['pcond_factor'])
        pcond_factor_err = np.array(cal_data['pcond_factor_err'])
        float_name = "3901960"
        profile_no = float_data['PROFILE_NO']

        config = TESTS_CONFIG
        config['FLOAT_PLOTS_DIRECTORY'] = "data/test_data/float_plots/"
        config['FLOAT_PLOTS_FORMAT'] = "eps"

        self.assertEqual(plots.cal_sal_curve_plot(sal, cal_sal, cal_sal_err, sta_sal,
                                                  sta_sal_err, sta_mean, pcond_factor,
                                                  pcond_factor_err, profile_no,
                                                  float_name, config),
                         None)


# pylint: disable=bare-except
# pylint: disable=unused-argument
# pylint: disable=too-many-locals
class SalVar(unittest.TestCase):
    """
    Test cases for sal_var_plot function
    """

    @patch("pyowc.plot.plots.plt.show")
    def test_plot_runs(self, mockshow):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that salinity variance plot throws no errors")

        # grid_data = loadmat("../../data/test_data/float_mapped_test/map_3901960.mat")
        # float_data = loadmat("../../data/float_source/3901960.mat")
        # cal_data = loadmat("../../data/test_data/float_calib_test/cal_3901960.mat")

        cal_data_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "float_calib_test",
                                          TESTS_CONFIG['FLOAT_CALIB_PREFIX'] +
                                          TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                          TESTS_CONFIG['FLOAT_CALIB_POSTFIX']])
        cal_data = loadmat(cal_data_path)
        float_data_path = os.path.sep.join([TESTS_CONFIG['FLOAT_SOURCE_DIRECTORY'],
                                            TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                            TESTS_CONFIG['FLOAT_SOURCE_POSTFIX']])
        float_data = loadmat(float_data_path)

        grid_data_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "float_mapped_test",
                                           TESTS_CONFIG['FLOAT_MAPPED_PREFIX'] +
                                           TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                           TESTS_CONFIG['FLOAT_MAPPED_POSTFIX']])
        grid_data = loadmat(grid_data_path)

        sal = float_data['SAL']
        theta = float_data['PTMP']
        pres = float_data['PRES']
        grid_sal = grid_data['la_mapped_sal']
        grid_ptmp = grid_data['la_ptmp']
        grid_errors = grid_data['la_mapsalerrors']
        cal_sal = cal_data['cal_SAL']
        cal_sal_errors = cal_data['cal_SAL_err']
        float_name = "3901960"

        no_boundaries = [[], [], [], [], 0.5]
        low_bound = [[250], [], [250], [], 0.5]
        up_bound = [[], [1], [], [1], 0.5]
        mid_bound = [[10], [-20], [1750], [500], 0.5]
        band_bound = [[-0.1], [0.1], [1000], [1100], 0.5]

        profile_no = float_data['PROFILE_NO']

        config = TESTS_CONFIG
        config['FLOAT_PLOTS_DIRECTORY'] = "data/test_data/float_plots/"
        config['FLOAT_PLOTS_FORMAT'] = "eps"

        # Check various types run

        self.assertEqual(plots.sal_var_plot(2, copy.deepcopy(sal), copy.deepcopy(pres),
                                            copy.deepcopy(theta), copy.deepcopy(grid_sal),
                                            copy.deepcopy(grid_errors), copy.deepcopy(grid_ptmp),
                                            copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_errors),
                                            no_boundaries,
                                            profile_no, float_name, config), None)
        self.assertEqual(plots.sal_var_plot(2, copy.deepcopy(sal), copy.deepcopy(pres),
                                            copy.deepcopy(theta), copy.deepcopy(grid_sal),
                                            copy.deepcopy(grid_errors), copy.deepcopy(grid_ptmp),
                                            copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_errors),
                                            low_bound,
                                            profile_no, float_name, config), None)
        self.assertEqual(plots.sal_var_plot(2, copy.deepcopy(sal), copy.deepcopy(pres),
                                            copy.deepcopy(theta), copy.deepcopy(grid_sal),
                                            copy.deepcopy(grid_errors), copy.deepcopy(grid_ptmp),
                                            copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_errors),
                                            up_bound,
                                            profile_no, float_name, config), None)
        self.assertEqual(plots.sal_var_plot(2, copy.deepcopy(sal), copy.deepcopy(pres),
                                            copy.deepcopy(theta), copy.deepcopy(grid_sal),
                                            copy.deepcopy(grid_errors), copy.deepcopy(grid_ptmp),
                                            copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_errors),
                                            band_bound,
                                            profile_no, float_name, config), None)
        self.assertEqual(plots.sal_var_plot(2, copy.deepcopy(sal), copy.deepcopy(pres),
                                            copy.deepcopy(theta), copy.deepcopy(grid_sal),
                                            copy.deepcopy(grid_errors), copy.deepcopy(grid_ptmp),
                                            copy.deepcopy(cal_sal), copy.deepcopy(cal_sal_errors),
                                            mid_bound,
                                            profile_no, float_name, config), None)


# pylint: disable=bare-except
# pylint: disable=unused-argument
# pylint: disable=too-many-locals
class TS(unittest.TestCase):
    """
    Test cases for t_s_plot function
    """

    @patch("pyowc.plot.plots.plt.show")
    def test_plot_runs(self, mockshow):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """
        print("Check t_s_profile_plot runs")
        # float_data = scipy.loadmat("../../data/float_source/3901960.mat")
        # grid_data = scipy.loadmat("../../data/test_data/float_mapped_test/map_3901960.mat")

        float_data_path = os.path.sep.join([TESTS_CONFIG['FLOAT_SOURCE_DIRECTORY'],
                                            TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                            TESTS_CONFIG['FLOAT_SOURCE_POSTFIX']])
        float_data = loadmat(float_data_path)

        grid_data_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "float_mapped_test",
                                           TESTS_CONFIG['FLOAT_MAPPED_PREFIX'] +
                                           TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                           TESTS_CONFIG['FLOAT_MAPPED_POSTFIX']])
        grid_data = loadmat(grid_data_path)

        sal = np.array(float_data['SAL'])
        ptmp = np.array(float_data['PTMP'])
        pres = float_data['PRES']
        grid_ptmp = grid_data['la_ptmp']

        thetas = find_10thetas(copy.deepcopy(sal),
                               copy.deepcopy(ptmp),
                               copy.deepcopy(pres),
                               copy.deepcopy(grid_ptmp),
                               [], [],
                               [], [],
                               0.5)

        sal_var = thetas[3]
        theta_levels = thetas[4]
        tlevels = thetas[0]
        plevels = thetas[1]

        config = TESTS_CONFIG
        config['FLOAT_PLOTS_DIRECTORY'] = "data/test_data/float_plots/"
        config['FLOAT_PLOTS_FORMAT'] = "eps"

        self.assertEqual(plots.t_s_profile_plot(sal, ptmp, pres, sal_var,
                                                theta_levels, tlevels, plevels,
                                                "3901960", config), None)


# pylint: disable=bare-except
# pylint: disable=unused-argument
class ThetaSal(unittest.TestCase):
    """
    Test cases for theta_sal_plot function
    """

    @patch("pyowc.plot.plots.plt.show")
    def test_plot_runs(self, mockshow):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        print("Test that theta salinity plot throws no errors")

        # grid_data = loadmat("../../data/test_data/float_mapped_test/map_3901960.mat")
        # float_data = loadmat("../../data/float_source/3901960.mat")

        float_data_path = os.path.sep.join([TESTS_CONFIG['FLOAT_SOURCE_DIRECTORY'],
                                            TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                            TESTS_CONFIG['FLOAT_SOURCE_POSTFIX']])
        float_data = loadmat(float_data_path)

        grid_data_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "float_mapped_test",
                                           TESTS_CONFIG['FLOAT_MAPPED_PREFIX'] +
                                           TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                           TESTS_CONFIG['FLOAT_MAPPED_POSTFIX']])
        grid_data = loadmat(grid_data_path)

        sal = np.array(float_data['SAL'])
        theta = np.array(float_data['PTMP'])
        pres = float_data['PRES']
        map_sal = grid_data['la_mapped_sal']
        map_ptmp = grid_data['la_ptmp']
        map_errors = grid_data['la_mapsalerrors']
        profiles = float_data['PROFILE_NO'][0]

        thetas = find_10thetas(copy.deepcopy(sal),
                               copy.deepcopy(theta),
                               copy.deepcopy(pres),
                               copy.deepcopy(map_ptmp),
                               [], [],
                               [], [],
                               0.5)

        index = thetas[2]

        config = TESTS_CONFIG
        config['FLOAT_PLOTS_DIRECTORY'] = "data/test_data/float_plots/"
        config['FLOAT_PLOTS_FORMAT'] = "eps"

        # Check various types run

        self.assertEqual(plots.theta_sal_plot(sal.transpose(), theta.transpose(),
                                              map_sal, map_ptmp, map_errors, index,
                                              profiles, config, "3901960"),
                         None)


# pylint: disable=bare-except
class Trajectory(unittest.TestCase):
    """
    Test cases for trajectory_plot function
    """

    def test_plot_runs(self):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """
        print("Test that trajectory plot throws no errors")

        config = TESTS_CONFIG
        config['FLOAT_PLOTS_DIRECTORY'] = "data/test_data/float_plots/"
        config['FLOAT_PLOTS_FORMAT'] = "eps"

        float_data_path = os.path.sep.join([TESTS_CONFIG['FLOAT_SOURCE_DIRECTORY'],
                                            TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                            TESTS_CONFIG['FLOAT_SOURCE_POSTFIX']])
        float_data = loadmat(float_data_path)

        grid_data_path = os.path.sep.join([TESTS_CONFIG['TEST_DIRECTORY'], "float_mapped_test",
                                           TESTS_CONFIG['FLOAT_MAPPED_PREFIX'] +
                                           TESTS_CONFIG['TEST_FLOAT_SOURCE'] +
                                           TESTS_CONFIG['FLOAT_MAPPED_POSTFIX']])
        grid_data = loadmat(grid_data_path)

        grid, floats = utils.create_dataframe(grid_data, float_data)

        try:
            plots.trajectory_plot(1, 1, floats, grid, "3901960", TESTS_CONFIG)

        except:
            self.fail("Trajectory plotting routine failed unexpectedly")


if __name__ == '__main__':
    unittest.main()
