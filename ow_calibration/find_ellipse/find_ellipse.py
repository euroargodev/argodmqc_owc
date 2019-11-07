"""
-----Find Ellipse-----

Written by: Edward Small
When: 7/11/2019

Calculates whether or not a value exists inside of an ellipse
of specifed size. If the answer is <1, it exists in the ellipse

Used to belong in "find_besthist", but was refactored and removed
to its own file for neatness.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import math

#pylint: disable=too-many-arguments
def find_ellipse(data_long, ellipse_long, ellipse_size_long,
                 data_lat, ellipse_lat, ellipse_size_lat,
                 phi, data_pv=0, ellipse_pv=0):
    """
    Finds whether a data point exists inside an ellipse
    :param data_long: longitude of the data point
    :param ellipse_long: longitude of the centre of the ellipse
    :param ellipse_size_long: size of the ellipse in the longitudinal direction
    :param data_lat: latitude of the data point
    :param ellipse_lat: latitude of the centre of the ellipse
    :param ellipse_size_lat: size of the ellipse in the latitudinal direction
    :param phi: cross-isobath scale for ellipse
    :param data_pv: potential vorticity of the data point
    :param ellipse_pv: potential vorticity of the centre of the ellipse
    :return: float. If <1, it exists inside the ellipse
    """

    total_pv = 0
    if data_pv != 0 and ellipse_pv != 0:
        total_pv = (ellipse_pv - data_pv) / math.sqrt(ellipse_pv ** 2 + data_pv ** 2) / phi

    ellipse = math.sqrt((data_long - ellipse_long) ** 2 / (ellipse_size_long * 3) ** 2 + \
                        (data_lat - ellipse_lat) ** 2 / (ellipse_size_lat * 3) ** 2 + \
                        total_pv ** 2)

    return ellipse
