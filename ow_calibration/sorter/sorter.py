"""
-----Sorter Function-----

Written by: Edward Small
When: 05/06/2020

Function to fill out arrays for the piecewise-fit

Used to find the interval in a linear fit

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np


def sorter(msites, sites):
    """
    Flag points as inside our outside boundaries msites[0] to msites[1]
    :param msites: boundaries
    :param sites: points
    :return: Array containing values defining points as inside or outside boundaries
    """

    # put arrays together
    all_sites = np.concatenate((msites.flatten(), sites.flatten()))

    index = np.argsort(all_sites, kind='stable')

    pointer = (np.argwhere(index > msites.__len__() - 1).flatten() - \
               np.arange(1, sites.__len__() + 1))

    return pointer