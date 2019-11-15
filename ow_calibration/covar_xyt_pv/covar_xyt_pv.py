"""
-----Covariance of 2D space and time (and potential vorticity)-----

Written by: A. Wong
When: 29/04/2005
Converted to python by: Edward Small
When: 13/11/2019

Calculates covariance of each point against every other point using the
Squared Exponential (SE) covariance function:

SE(x,y) = exp(-(x-y)^2/2l) where (x-y) is the difference between points (could be distance,
time, etc), and l is the characteristic length scale (how close points have to be
to influence each other significantly).

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import math
import numpy as np


# pylint: disable=too-many-arguments
def covar_xyt_pv(points1, points2, lat, long, age, phi, p_v):
    """
    Calculates how "close" two sets of points are to each other, taking into account
    space, time and (if wanted) potential vorticity. The closer two values are to
    each other, the closer their value will be to 1. Points that differ greatly will be
    nearer to 0.
    :param points1: m*4 matrix containing the latitude, longitude, date,
    and depth of each data point
    :param points2: n*4 matrix containing the latitude, longitude, date,
    and depth of each data point
    :param lat: float, the characteristic latitude
    :param long: float, the characteristic longitude
    :param age: float, the characteristic time scale
    :param phi: float, the characteristic cross-isobaric scale (for depth dependence)
    :param p_v: int, flag for using vorticity (1=include)
    :return: m*n matrix containing the covariance of each point
    """

    # create the m*n covariance matrix filled with 0's
    try:
        points_covar = np.full((points1.__len__(), points2.__len__()), 0, float)

    except AttributeError:
        raise AttributeError("A set of points has no length associated with it")

    for i in range(0, points1.__len__()):
        for j in range(0, points2.__len__()):

            # calculate the absolute difference between points over chracteristic length scale
            lat_covar = ((points1[i][0] - points2[j][0]) / lat) ** 2
            long_covar = ((points1[i][1] - points2[j][1]) / long) ** 2
            age_covar = 0
            p_v_covar = 0

            if age != 0:
                age_covar = ((points1[i][2] - points2[j][2]) / age) ** 2

            # pylint: disable=fixme
            # TODO: ARGODEV-163
            # use the potential vorticity function made in ARGODEV-150
            if p_v == 1:
                print("pv not yet included. Phi: ", phi)

            points_covar[i][j] = math.exp(-(lat_covar + long_covar + age_covar + p_v_covar))

    return points_covar
