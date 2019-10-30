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

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
import math


def calc_barotropic_potential_vorticity(lat, z):
    """
    Calculates potential vorticity (PV)
    :param lat: latitude
    :param z: depth
    """

    return (2 * 7.292 * 10 ** -5 * math.sin(lat * math.pi / 180)) / z
