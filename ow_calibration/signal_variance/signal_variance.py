"""
-----Signal Variance-----

Written by: X X
When: xx/xx/xxxx
Converted to python by: Edward Small
When: 20/11/2019

Calculates an estimate of the signal variance at a given level using:

(sum({di - <d>}^2))/N

where di is a data point in d (which is a collection of salinities at a given level,
 <d> is the mean of all the data points in d, and N is the number of data points
 in d.

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import math
import numpy as np


def signal_variance(sal):
    """
    Calculates signal variance
    :param sal: vector of salinities at a given level
    :return: float estimate of the variance of the signal of the given data
    """

    # remove 0's and NaN's from our data set
    sal_no_nan = []
    for i in range(0, sal.__len__()):
        if sal[i] != 0 and not math.isnan(sal[i]):
            sal_no_nan.append(sal[i])

    # check that we have got some valid values. If not, raise an exception
    num_sal = sal_no_nan.__len__()
    if num_sal == 0:
        raise Exception("Received no valid salinity values when calculating signal variance")

    # approximate the signal variance
    sal_mean = np.mean(sal_no_nan)
    signal = (sal_no_nan - sal_mean)
    variance = np.sum(signal ** 2) / num_sal

    return variance
