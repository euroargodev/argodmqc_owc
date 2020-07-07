""" Core utility functions
"""

import math
import numpy as np


def sorter(msites, sites):
    """ Function to fill out arrays for the piecewise-fit

        Used to find the interval in a linear fit.

        Flag points as inside our outside boundaries msites[0] to msites[1]

        Parameters
        ----------
        msites: boundaries
        sites: points

        Returns
        -------
        Array containing values defining points as inside or outside boundaries
    """

    # put arrays together
    all_sites = np.concatenate((msites.flatten(), sites.flatten()))

    index = np.argsort(all_sites, kind='stable')

    pointer = (np.argwhere(index > msites.__len__() - 1).flatten() - \
               np.arange(1, sites.__len__() + 1))

    return pointer


def wrap_longitude(grid_long):
    """ Allows longitude (0-360) to wrap beyond the 360 mark, for mapping purposes.

        Makes sure that, if the longitude is near the boundary (0 or 360) that we
        wrap the values beyond
        360 so it appears nicely on a map

        This is a refactor between get_region_data and get_region_hist_locations to
        avoid duplicate code

        Parameters
        ----------
        grid_long: array of longitude values

        Returns
        -------
        array of longitude values that can extend past 360
    """

    neg_long = np.argwhere(grid_long < 0)
    grid_long[neg_long] = grid_long[neg_long] + 360

    # if we have data close to upper boundary (360), then wrap some of the data round
    # so it appears on the map
    top_long = np.argwhere(grid_long >= 320)
    if top_long.__len__() != 0:
        bottom_long = np.argwhere(grid_long <= 40)
        grid_long[bottom_long] = 360 + grid_long[bottom_long]

    return grid_long


def potential_vorticity(lat, z_value):
    """ Calculates barotropic potential vorticity (pv)

        Calculates the potential vorticity for a given latitude and z

        Used to belong in "find_besthist", but was refactored and removed
        to its own file for neatness.

        Parameters
        ----------
        lat: latitude

        Returns
        -------
        z_value: depth
    """
    earth_angular_velocity = 2 * 7.292 * 10 ** -5
    lat_radians = lat * math.pi / 180

    p_v = (earth_angular_velocity * math.sin(lat_radians)) / z_value

    if p_v == 0:
        p_v = 1 * 10 ** -5

    return p_v


def cal2dec(pa_month, pa_day, pa_hour=0, pa_minute=0):
    """ Converts a calendar date (month, day, hour, minute) to a decimal date (float)

        Parameters
        ----------
        pa_month: Month in the year (where 0 is Janurary and 11 is Decemeber)
        pa_day: Day in the month
        pa_hour: Hour in the day
        pa_minute: Minute in the hour

        Returns
        -------
        decimalised version of the date
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


def change_dates(cal_dates):
    """ Seperates date into year, month, day, hour, minute

        Converts calendar date from one input into year, month, day, hour, minute.
        Passes these variables to cal2dec to get a decimalised date and returns it

        N.B. Change to code on the xx/06/2013: Speed up this function - C Cabanes

        Parameters
        ----------
        cal_dates: array of date as one continuous int
        eg (01/02/1974, hour 3, minute 5, second 44 = 19740102030544)

        Returns
        -------
        decimalised date
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
