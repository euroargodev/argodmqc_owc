"""
-----noise_variance-----

Written by:
When:
Converted to python by: Edward Small
When: 19/11/2019

Finds the variance for the noise variance of each salinity measurement by comparing it
to the noise_variance of all the other measurements. Can be thought of as the average
difference between measured salinity and expected salinity:

noise_variance = (sum(x - u)^2) / 2*N where
x is the current observation
u is the closest observation (spatially)
N is the number of elements

This is because we assume that the noise variance is uncorrelated over distance, that is has
uniform variance, and that the signal has a longer correlation distance than the
data separation (WJO 2003, Delayed-Mode Calibration of Autonomous CTD Profiling
Float Salinity Data by θS Climatology).

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np


def noise_variance(sal, lat, long):
    """
    calculates the variance in noise_variance of salinity at different pressures
    :param sal: m*n matrix containing m layers and n casts
    :param lat: vector of latitudes of each cast
    :param long: vector of longitudes of each cast
    :return: the variance in the noise_variance
    """
    # set up our numpy matrix for memory efficiency
    sal_noise = np.zeros(sal.shape, dtype=float)

    # find the difference in salinities between each point and its closest other point (x-u)
    for i in range(0, sal.__len__()):
        # find the nearest spatial point to lat[i], long[i]
        distances = (long - long[i]) ** 2 + (lat - lat[i]) ** 2

        # find the smallest distance between this point and another point
        # we do this by first finding the instances where there is some distance
        # between points (> 0), and then finding the minimum of these instances
        min_distance = np.min(distances[np.nonzero(distances)])

        # find index of the minimum distance
        min_index = np.argwhere(distances == min_distance)

        # store the differences in salinities between these two points
        sal_noise[i] = sal[i] - sal[min_index]

    # make sure we have unique points
    index = np.argwhere(sal_noise != 0)

    # find the variance in the noise by summing the difference squared
    # and dividing it
    sal_noise_var = sum(sal_noise[index] ** 2) / (2 * index.__len__())

    return sal_noise_var
