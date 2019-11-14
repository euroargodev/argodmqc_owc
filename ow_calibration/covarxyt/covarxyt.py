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

