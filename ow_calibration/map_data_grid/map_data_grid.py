"""
-----Map Data to Grid-----

Written by: Breck Owens
When: xx/11/2007
Converted to python by: Edward Small
When: 22/11/2019

An optimal mapping routine, taking data measured in arbitrary geographic locations and
mapping these data onto a more regular grid. As should happen in every mapping problem
the data are both mapped onto the prescribed grid, and the data locations (so that
the mapped field can be checked with the original data to ensure that the statistics
are valid and consistent).

Before objectively mapping the data, a mean, using the correlation scales to define the
weights for the sum (see Bretherton, etal, 1975) is removed.  The error estimate includes
the contributions from both the mapping and the mean estimation.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
from ow_calibration.map_data_grid.covar_xyt_pv.covar_xyt_pv import covar_xyt_pv


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
def map_data_grid(sal, grid_pos, data_pos, lat, long, age,
                  signal_variance, noise_variance, phi, map_pv_use):
    """
    maps historical float data onto a single float
    :param sal: array of salinities of the historical float data
    :param grid_pos: array containing single float data [lat, long, age, depth]
    :param data_pos: n*4 array containing historical float data [lat, long, age, depth]
    :param lat: scalar latitude
    :param long: scalar longitude
    :param age: scalar age
    :param signal_variance: scalar signal variance
    :param noise_variance: scalar noise variance
    :param phi: scalar cross isobaric scale
    :param map_pv_use: flag for including vorticity (1=include)
    :return: tuple containing mapped fields, error estimates of mapped fields,
             mapped fields on original locations, and their error estimates
    """

    # create the data-data covariance matrix
    data_pos_covar = covar_xyt_pv(data_pos, data_pos, lat, long, age, phi, map_pv_use)
    data_data_covar = np.linalg.inv(signal_variance * data_pos_covar +
                                    noise_variance * np.identity(data_pos.__len__()))

    # estimate the mean field and weights
    sum_data_data_covar = sum(sum(data_data_covar))
    mean_field = sum(np.dot(data_data_covar, sal.transpose())) / sum_data_data_covar
    weight = np.dot(data_data_covar, (sal - mean_field))

    # calculate the objectively mapped fields on position data and grid data
    pos_data_covar = signal_variance * data_pos_covar
    data_weight_covar = np.dot(pos_data_covar, weight) + mean_field

    # include the error in the mean (from Brethrerton, 1975)
    dot_covar_diag = np.diag(np.dot(
        np.dot(pos_data_covar, data_data_covar),
        np.transpose(pos_data_covar)))
    covar_sum = np.sum(np.dot(pos_data_covar, data_data_covar), axis=1)
    data_weight_covar_error = np.sqrt(signal_variance - dot_covar_diag +
                                      ((1 - covar_sum) ** 2) / sum_data_data_covar)

    # now map to the data to the regular grid
    grid_data_covar = (signal_variance * covar_xyt_pv(data_pos, grid_pos, lat, long,
                                                      age, phi, map_pv_use)).transpose()
    grid_weight_covar = np.dot(grid_data_covar, weight) + mean_field
    dot_covar_diag = np.diag(np.dot(
        np.dot(grid_data_covar, data_data_covar), np.transpose(grid_data_covar)))
    covar_sum = np.sum(np.dot(grid_data_covar, data_data_covar), axis=1)
    grid_weight_covar_error = np.sqrt(signal_variance - dot_covar_diag +
                                      ((1 - covar_sum) ** 2) / sum_data_data_covar)

    return grid_weight_covar[0], grid_weight_covar_error[0], \
           data_weight_covar, data_weight_covar_error
