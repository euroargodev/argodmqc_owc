"""
-----find best historical casts-----

Written by: Breck Owens and Annie Wong
When: xx/12/2006
Converted to python by: Edward Small
When: 30/10/2019

Find ln_max_casts unique historical points that are most strongly correlated with the float profile

Rewritten in December 2006 to only use latitude, longitude, data, and water depth as arguments
and to return index of the station list

N.B. Change to code on the xx/06/2013: add age_large when computing correlation_large
- C Cabanes

N.B Change during conversion to python on the 31/10/2019: Potential Vorticity and Correlation
are calculated multiple times, so I moved them into their own function. These functions can be
vectorised using the numpy library (numpy.vectorize(function)) to use the functions on arrays.
- Edward Small

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
import math


def barotropic_potential_vorticity(lat, z):
    """
    Calculates barotropic potential vorticity (pv)
    :param lat: latitude
    :param z: depth
    """
    earth_angular_velocity = 2 * 7.292 * 10 ** -5
    lat_radians = lat * math.pi / 180

    pv = (earth_angular_velocity * math.sin(lat_radians)) / z

    if pv == 0:
        pv = 1 * 10 ** -5

    return pv


def spatial_correlation(
        longitude_1, longitude_2, ellipse_longitude, latitude_1, latitude_2,
        ellipse_latitude, dates_1, dates_2, ellipse_age, phi, pv_1=0, pv_2=0
):
    """
    Calculates the spatial correlation between two points.
    Can be done with or without potential vorticity
    :param longitude_1: longitude of point 1
    :param longitude_2: longitude if point 2
    :param ellipse_longitude: longitudinal size of ellipse
    :param latitude_1: latitude of point 1
    :param latitude_2: latitude of point 2
    :param ellipse_latitude: latitudinal size of ellipse
    :param dates_1: dates of data for point 1
    :param dates_2: dates of data for point 2
    :param ellipse_age: age of data wanted in ellipse
    :param phi: potential gradient
    :param pv_1: potential vorticity of point 1
    :param pv_2: potential vorticity of point 2
    :return: spatial correlation between points
    """
    if pv_1 == 0 or pv_2 == 0:
        correlation = (longitude_1 - longitude_2) ** 2 / ellipse_longitude + \
                      (latitude_1 - latitude_2) ** 2 / ellipse_latitude + \
                      (dates_1 - dates_2) ** 2 / ellipse_age

    if pv_1 != 0 or pv_2 != 0:
        correlation = correlation + ((pv_2 - pv_1) / math.sqrt(pv_2**2+pv_1**2)/phi)**2

    return correlation
