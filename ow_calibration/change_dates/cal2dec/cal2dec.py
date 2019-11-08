"""
-----Cal 2 Dec-----

Written by: P. Robbins
When: xx/xx/1994
Converted to python by: Edward Small
When: 27/09/2019

Converts a calendar date (month, day, hour, minute) to a decimal date (float)

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np


def cal2dec(pa_month, pa_day, pa_hour=0, pa_minute=0):
    """
    :param pa_month: Month in the year (where 0 is Janurary and 11 is Decemeber)
    :param pa_day: Day in the month
    :param pa_hour: Hour in the day
    :param pa_minute: Minute in the hour
    :return: decimalised version of the date
    """

    ln_cumulative_months = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

    try:
        dec_date = float(
            ln_cumulative_months[pa_month] + pa_day - 1 + pa_hour / 24 + pa_minute / 60 / 24
        )
        if dec_date > 366:
            raise ValueError("Day is out of scope of the year")
        return dec_date

    except IndexError:
        raise ValueError('Month is out of scope')