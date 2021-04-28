""" Global import controler
"""

#pylint: disable=broad-except
try:
    import pkg_resources
    __version__ = pkg_resources.get_distribution("pyowc").version
except Exception:
    # Local copy, not installed with setuptools, or setuptools is not available.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"

from . import core
from . import data
from . import plot
from . import utilities
from . import helper
from .plot import dashboard

from . import calibration
from . import configuration

#
__all__ = (
    # Sub-packages,
    "configuration",
    "calibration",
    "core",
    "data",
    "plot",
    "dashboard",
    "utilities",
    # Constants
    "__version__"
)
