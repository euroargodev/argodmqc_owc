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
import scipy.interpolate as interpolate
from ow_calibration.build_cov.covarxy_pv.covarxy_pv import covarxy_pv


def build_cov(ptmp, coord_float, config):
    """
    Build the covariance matrix
    :param ptmp: matrix of potential temperatures
    :param coord_float_config: the x, y, z position of the float
    :return: covariance matrix
    """

    # Set up theta boundaries for water masses

    ptboundaries = np.array([30, 24, 18, 12, 8, 4, 2.5, 1, -2])
    ptscale_down = np.array([6, 6, 6, 4, 4, 1.5, 1.5, 1, 1])
    ptscale_up = np.array([6, 6, 6, 6, 4, 4, 1.5, 1.5, 1])

    # Set up the building tile = vertical covariance matrix

    # The upper triangle of the matrix (top right values) are the covariance of each ptlevel with
    # every ptlevel below it, looking down the water column from the diagonal
    # The lower triangle of the matrix (bottom left values) are the covariance of each ptlevel with
    # every ptlevel above it, looking up the water column from the diagonal

    ptmp_rows = ptmp.shape[0]
    ptmp_columns = ptmp.shape[1]

    # set up the covariance matrix
    cov = np.zeros((ptmp_rows * ptmp_columns, ptmp_rows))

    # set up interpolation grids
    upper_interp = interpolate.interp1d(ptboundaries, ptscale_down)
    lower_interp = interpolate.interp1d(ptboundaries, ptscale_up)

    # go through each profile
    for profile in range(0, ptmp_columns):

        profile_1 = profile * ptmp_rows

        # go through each level
        for i in range(0, ptmp_rows):
            for j in range(0, ptmp_rows):

                # belongs in the upper triangle, look down water column for vertical scale
                if i < j:
                    l_theta = upper_interp(ptmp[i, profile])
                    cov[i + profile_1, j] = np.exp(-1 * (ptmp[j, profile] - ptmp[i, profile]) ** 2
                                                   / l_theta ** 2)

                # belongs in the lower triangle, look up water column for vertical scale
                elif i > j:
                    l_theta = lower_interp(ptmp[i, profile])
                    cov[i + profile_1, j] = np.exp(-1 * (ptmp[j, profile] - ptmp[i, profile]) ** 2
                                                   / l_theta ** 2)

                # it is in the leading diagonal, so make it equal to 1
                else:
                    cov[i + profile_1, j] = 1

                # if we don't have a value, make it equal to 1
                if np.isnan(cov[i + profile_1, j]):
                    cov[i + profile_1, j] = 1

    # set up matrix to hold horizontal covariance
    h_cov = np.ones((ptmp_columns, ptmp_columns)) * np.nan

    for profile in range(0, ptmp_columns):
        h_cov[profile, :] = covarxy_pv(coord_float[profile], coord_float,
                                       config['MAPSCALE_LONGITUDE_SMALL'],
                                       config['MAPSCALE_LATITUDE_SMALL'],
                                       config['MAPSCALE_PHI_SMALL'],
                                       config['MAPSCALE_USE_PV'])

    h_cov = h_cov[:, 0:ptmp_columns]

    # build final covariance matrix, using horizontal and vertical covariance

    n_cov = np.tile(cov, [1, ptmp_columns])

    # Have to find the covariance for each profile against all other profiles
    for profile in range(0, ptmp_columns):

        lower = profile * ptmp_rows
        upper = (profile + 1) * ptmp_rows

        # go through each profile
        for profile_1 in range(0, ptmp_columns):
            lower_1 = profile_1 * ptmp_rows
            upper_1 = (profile_1 + 1) * ptmp_rows
            n_cov[lower:upper, lower_1:upper_1] = h_cov[profile, profile_1] * \
                                                  n_cov[lower:upper, lower_1:upper_1]

    return n_cov
