"""
-----Change Dates-----

Written by: A. Wong
When: 29/05/2001
Converted to python by: Edward Small
When: 08/11/2019

Converts calendar date from one input into year, month, day, hour, minute.
Passes these variables to cal2dec to get a decimalised date and returns it

N.B. Change to code on the xx/06/2013: Speed up this function - C Cabanes

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
from ow_calibration.change_dates.cal2dec import cal2dec


def change_dates(cal_dates):
    """
    Seperates date into year, month, day, hour, minute
    :param cal_dates: array of date as one continuous int
    eg (01/02/1974, hour 3, minute 5, second 44 = 19740102030544)
    :return: decimalised date
    """

    # might need to return zeroes for bad dates
    dec_dates = np.full(len(cal_dates), 0)

    # cannot interate through integers, so convert calendar dates into array of strings
    cal_dates_string = list(map(str, cal_dates))


change_dates([123, 456, 789])
