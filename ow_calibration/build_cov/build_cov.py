"""
-----Build Covariance Matrix-----

Written by: Annie Wong
When: xx/05/2001
Converted to python by: Edward Small
When: 05/03/2020

Function builds a square covariance matrix that has n*n tiles, and each
tile is of size m*m

The vertical covariance is the building tile. It contains a 1 down its diagonal,
which decreases exponentially in the off-diagonals (representing the vertical covariance
between water masses).

We assume that each profile is independent from the other profiles

N.B. Change to code on the xx/11/2007: Undefined - Breck Owens

N.B. Change to code on the xx/06/2013: Take into account the horizontal covariance between
different mapped profiles. This lateral covariance takes into account the fact that a mapped
profile on an Argo position is built from a set of historical profiles that are not very
different from the set used to build a mapped profile for the next or previous Argo position.
This lateral covariance between two mapped profiles is constructed using a Guassian function
and the large spatial scales - Cecile Cabanes

N.B. Change to code on the xx/xx/2017: Use the small spatial scales instead of the large spatial
scales to build the lateral covariance. This has been found to be the best comprimise between the
informative errors and large enough NDF for AIC criterion, at least for the scales defined for the
North Atlantic basin - Cecile Cabanes

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np


def build_cov(ptmp, coord_float, config):
    """
    Build the covariance matrix
    :param ptmp: matrix of potential temperatures
    :param coord_float_config: the x, y, z position of the float
    :return: covariance matrix
    """

    # Set up theta boundaries for water masses

    ptboundaries = np.array([30, 24, 18, 12, 8, 4, 2.5, 1, -2])