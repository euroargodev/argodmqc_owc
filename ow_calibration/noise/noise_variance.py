"""
-----noise variance-----

Written by:
When:
Converted to python by: Edward Small
When: 19/11/2019

Finds the variance for the noise of each salinity measurement by comparing it
to the noise of all the other measurements. Can be thought of as the average
difference between measured salinity and expected salinity.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np


def noise_variance(sal, lat, long):
    """
    calculates the variance in noise of salinity at different pressures
    :param sal: m*n matrix containing m layers and n casts
    :param lat: vector of latitudes of each cast
    :param long: vector of longitudes of each cast
    :return: the variance in the noise
    """
    # set up our numpy matrix for memory efficiency
    sal_diff = np.zeros(sal.shape, dtype=float)

    # find the difference in salinities between each point and its closest other point
    for i in range(0, lat.__len__()):
        # find the nearest spatial point to lat[i], long[i]
        distances = (long - long[i]) ** 2 + (lat - lat[i]) ** 2

        # find the smallest distance between this point and another point
        # we do this by first finding the instances where there is some distance
        # between points (> 0), and then finding the minimum of these instances
        min_distance = np.min(distances[np.nonzero(distances)])

        # find index of the minimum distance
        min_index = np.argwhere(distances == min_distance)

        # store the differences in salinities between these two points
        sal_diff[i] = sal[i] - sal[min_index]


sal = np.array([34.4988, 34.3267, 34.0346])
lat = np.array([-57.9960, -56.4902, -52.0800])
long = np.array([53.1950, 63.1170, 61.9770])

noise_variance(sal, lat, long)
