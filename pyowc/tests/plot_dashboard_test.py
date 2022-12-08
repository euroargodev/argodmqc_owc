import os
import unittest
from unittest.mock import patch

from pyowc.plot import dashboard
from . import TESTS_CONFIG


class Dashboard(unittest.TestCase):
    """
    Test dashboard function
    plot_diagnostics(float_dir, float_name, config, levels=2)

    """

    @patch("pyowc.plot.plots.plt.show")
    def test_dashboard(self, mockshow):
        print("Test that dashboard plot throws no errors")
        float_dir = "/"
        float_name = TESTS_CONFIG['TEST_FLOAT_SOURCE']
        config = TESTS_CONFIG
        self.assertEqual(dashboard(float_dir, float_name, config, levels=2), None)
