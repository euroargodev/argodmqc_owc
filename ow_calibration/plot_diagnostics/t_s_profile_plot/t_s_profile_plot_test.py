"""
-----profile plot test file-----

Written by: Edward Small
When: 11/05/2019

Contains unit tests to check the functionality of the `t_s_profile_plot` function

To run this test specifically, look at the documentation at:
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import unittest
from unittest.mock import patch
import scipy.io as scipy
import numpy as np
from ow_calibration.plot_diagnostics.t_s_profile_plot.t_s_profile_plot import t_s_profile_plot


# pylint: disable=bare-except
# pylint: disable=unused-argument
# pylint: disable=too-many-locals
class MyTestCase(unittest.TestCase):
    """
    Test cases for t_s_plot function
    """

    #@patch("ow_calibration.plot_diagnostics.sal_var_plot.sal_var_plot.plt.show")
    #def test_plot_runs(self, mockshow):
    def test_plot_runs(self):
        """
        Check we get no errors during the plotting routine
        :return: nothing
        """

        t_s_profile_plot()


if __name__ == '__main__':
    unittest.main()