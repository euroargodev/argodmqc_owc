"""
Global import controler
"""

try:
    import pkg_resources
    __version__ = pkg_resources.get_distribution("pyowc").version
except ImportError:
    # Local copy, not installed with setuptools, or setuptools is not available.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"
