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


# pylint: disable=too-many-arguments
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
