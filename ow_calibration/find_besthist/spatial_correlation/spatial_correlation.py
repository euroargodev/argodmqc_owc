"""
-----Spatial correlation-----

Written by: Edward Small
When: 7/11/2019

Calculates how closely correlated two points are in space and time

Used to belong in "find_besthist", but was refactored and removed
to its own file for neatness.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import math


# pylint: disable=too-many-arguments
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

    pv_correlation = 0
    correlation = ((longitude_1 - longitude_2) ** 2 / (ellipse_longitude ** 2)) + \
                  ((latitude_1 - latitude_2) ** 2 / (ellipse_latitude ** 2)) + \
                  ((dates_1 - dates_2) ** 2 / (ellipse_age ** 2))

    if pv_1 != 0 and pv_2 != 0:
        pv_correlation = ((pv_2 - pv_1) / math.sqrt(pv_2 ** 2 + pv_1 ** 2) / phi) ** 2

    return correlation + pv_correlation
