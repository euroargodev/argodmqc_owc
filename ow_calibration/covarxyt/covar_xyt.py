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

import numpy as np
import math


def covar_xyt(x, y, lat, long, age, phi, p_v):
    """
    Calculates how "close" two sets of points are to each other, taking into account
    space, time and (if wanted) potential vorticity. The closer two values are to
    each other, the closer their value will be to 1. Points that differ greatly will be
    nearer to 0.
    :param x: m*4 matrix containing the latitude, longitude, date, and depth of each data point
    :param y: n*4 matrix containing the latitude, longitude, date, and depth of each data point
    :param lat: float, the characteristic latitude
    :param long: float, the characteristic longitude
    :param age: float, the characteristic time scale
    :param phi: float, the chracteristic cross-sobaric scale (for depth dependence)
    :param p_v: int, flag for using vorticity (1=include)
    :return: m*n matrix containing the covariance of each point
    """

    # create the m*n covariance matrix filled with 0's
    grid_covar = np.full((x.__len__(), y.__len__()), 0, float)

    for i in range(0, x.__len__()):
        for j in range(0, y.__len__()):

            # calculate the absolute difference between points over length scale
            lat_covar = ((x[i][0] - y[j][0]) / lat) ** 2
            print("latitude: ", lat_covar)
            long_covar = ((x[i][1] - y[j][1]) / long) ** 2
            print("longitude: ", long_covar)
            age_covar = 0
            p_v_covar = 0

            if age != 0:
                age_covar = ((x[i][2] - y[j][2]) / age) ** 2
                print("age: ", age_covar)

            # TODO: ARGODEV-163
            # use the potential vorticity function made in ARGODEV-150
            if p_v == 1:
                print("pv not yet included")

            grid_covar[i][j] = math.exp(-(lat_covar + long_covar + age_covar + p_v_covar))

    return grid_covar


x = np.array([[-0.057996, 0.053195, 1.9740875, 5.229838],
              [-0.0564902, 0.0631170, 1.9870367, 4.6300392],
              [-0.05208, 0.0619770, 1.9941118, 4.6536932]]) * 10 ** 3
y = x
lat = 4
long = 8
age = 20
phi = 0.5
p_v = 0

ans = covar_xyt(x, y, lat, long, age, phi, p_v)

print(ans)
