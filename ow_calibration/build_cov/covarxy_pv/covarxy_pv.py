"""
-----Build Covariance Matrix-----

Written by: Edward Small
When: 10/03/2020

Finds the correlation between spatial and temporal data, and uses this
to construct the covariance

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np

# pylint: disable=too-many-arguments
def covarxy_pv(input_coords, coords, long, lat, phi, use_pv):
    """
    Returns a matrix for the horizontal covariance
    :param input_coords: the input coordinates of the the float profile
    :param coords: coordinates for all the float profiles
    :param long: longitude scale
    :param lat: latitude scale
    :param phi: potential gradient
    :param use_pv: whether or not to use potential vorticity
    :return: horizontal covariance matrix
    """

    # Derive the planetary vorticity at each point

    # Get the depth for each data point
    z_input_coords = input_coords[2]
    z_coords = coords[:, 2]

    # define a vectorized function to calculation potential vorticity
    potential_vorticity = np.vectorize(lambda latitude, depth:
                                       (2 * 7.292 * 10 ** -5 *
                                        np.sin(latitude * np.pi / 180)) / depth)

    # calculate potential vorticity
    pv_input_coords = potential_vorticity(input_coords[0], z_input_coords)
    pv_coords = potential_vorticity(coords[:, 0], z_coords)

    # calculate correlation
    cor_term = ((input_coords[0] - coords[:, 0]) / lat) ** 2 + \
               ((input_coords[1] - coords[:, 1]) / long) ** 2

    # include potential vorticity in correlation, if the user has asked for it

    if use_pv and pv_input_coords.any() and pv_coords.any() != 0:
        cor_term = cor_term + ((pv_input_coords - pv_coords) /
                               np.sqrt(pv_input_coords ** 2 + pv_coords ** 2) /
                               phi) ** 2

    cov_term = np.exp(-cor_term.transpose())

    return cov_term
