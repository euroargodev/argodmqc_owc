"""
-----Potential Vorticity-----

Written by: Edward Small
When: 7/11/2019

Calculates the potential vorticity for a given latitude and z

Used to belong in "find_besthist", but was refactored and removed
to its own file for neatness.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import math


def potential_vorticity(lat, z_value):
    """
    Calculates barotropic potential vorticity (pv)
    :param lat: latitude
    :param z_value: depth
    """
    earth_angular_velocity = 2 * 7.292 * 10 ** -5
    lat_radians = lat * math.pi / 180

    p_v = (earth_angular_velocity * math.sin(lat_radians)) / z_value

    if p_v == 0:
        p_v = 1 * 10 ** -5

    return p_v
