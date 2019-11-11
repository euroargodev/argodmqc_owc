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
from ow_calibration.change_dates.cal2dec.cal2dec import cal2dec


def change_dates(cal_dates):
    """
    Seperates date into year, month, day, hour, minute
    :param cal_dates: array of date as one continuous int
    eg (01/02/1974, hour 3, minute 5, second 44 = 19740102030544)
    :return: decimalised date
    """

    # might need to return zeroes for bad dates
    dec_dates = np.full(cal_dates.__len__(), 0, dtype=float)

    # cannot iterate through integers, so convert calendar dates into array of strings
    cal_dates_string = list(map(str, cal_dates))

    # go through and decimalise dates
    for i in range(0, cal_dates_string.__len__()):
        try:
            # select the separate date entities (year, month, day, etc)
            hour = 0
            minute = 0
            year = int(cal_dates_string[i][0:4])
            month = int(cal_dates_string[i][4:6])
            day = int(cal_dates_string[i][6:8])

            if cal_dates_string[i].__len__() > 9:
                hour = int(cal_dates_string[i][8:10])

            if cal_dates_string[i].__len__() > 11:
                minute = int(cal_dates_string[i][10:12])

            if 13 > month > 0:
                if 32 > day > 0:
                    day = year + (cal2dec(month - 1, day, hour, minute) / 365)
                    dec_dates[i] = day

        except ValueError:
            print("Date is incorrect length or format")
            continue

    return dec_dates
