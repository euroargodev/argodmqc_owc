import unittest
from unittest.mock import patch

import numpy as np

from pyowc.helper import create_la_wmo_boxes_file

from . import TESTS_CONFIG


class CreateWMOBoxesTest(unittest.TestCase):
    """Tests for create_la_wmo_boxes_file."""

    @patch("pyowc.helper.savemat")
    def test_create_wmo_boxes_file(self, save_mat):
        create_la_wmo_boxes_file(TESTS_CONFIG)
        save_mat.assert_called_once()
