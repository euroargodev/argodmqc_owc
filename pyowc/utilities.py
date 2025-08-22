"""Core utility functions"""

import json
import math

import numpy as np
from pydantic import BaseModel, model_validator


class ConfigFileNotFoundError(Exception):
    """Raised when the config file is not found."""

class Config(BaseModel):
    """A representation of a config, to ensure correct types are passed to the rest of the app."""

    HISTORICAL_DIRECTORY: str
    HISTORICAL_CTD_PREFIX: str
    HISTORICAL_BOTTLE_PREFIX: str
    HISTORICAL_ARGO_PREFIX: str

    FLOAT_SOURCE_DIRECTORY: str
    FLOAT_SOURCE_POSTFIX: str
    FLOAT_MAPPED_DIRECTORY: str
    FLOAT_MAPPED_PREFIX: str
    FLOAT_MAPPED_POSTFIX: str

    FLOAT_CALIB_DIRECTORY: str
    FLOAT_CALIB_PREFIX: str
    FLOAT_CALSERIES_PREFIX: str
    FLOAT_CALIB_POSTFIX: str

    FLOAT_PLOTS_DIRECTORY: str
    FLOAT_PLOTS_FORMAT: str

    CONFIG_DIRECTORY: str
    CONFIG_COASTLINES: str
    CONFIG_WMO_BOXES: str
    CONFIG_SAF: str

    CONFIG_MAX_CASTS: int = 300
    MAP_USE_PV: int = 0
    MAP_USE_SAF: int = 0

    MAPSCALE_LONGITUDE_LARGE: int = 8
    MAPSCALE_LONGITUDE_SMALL: int = 4
    MAPSCALE_LATITUDE_LARGE: int = 4
    MAPSCALE_LATITUDE_SMALL: int = 2

    MAPSCALE_PHI_LARGE: float = 0.1
    MAPSCALE_PHI_SMALL: float = 0.02

    MAPSCALE_AGE_LARGE: int = 20
    MAPSCALE_AGE_SMALL: int = 10

    MAP_P_EXCLUDE: int = 100
    MAP_P_DELTA: int = 250

    THETA_BOUNDS: list[list[int]]


    @model_validator(mode="before")
    def check_for_empty_strings(cls, values):
        """Check for any empty strings."""
        for key, value in values.items():
            if isinstance(value, str) and not value:
                raise ValueError(f"The key {key} has an empty value in your config file!")
        return values


def load_configuration_from_json_file(filepath: str) -> dict:
    """Read the config JSON file and return the contents as a dict.

    Function will raise a ConfigFileNotFoundError for an invalid path.

    Args:
        filepath(str): The path to the JSON file on the users machine.

    Returns:
        A dictionary of values to run the dmqc.

    """
    try:
        with open(filepath) as file:
            model = Config(**json.loads(file.read()))
            return model.model_dump()
    except FileNotFoundError:
        raise ConfigFileNotFoundError(f"Config file not found: - {filepath}") from None

# pylint: disable=too-many-arguments
def spatial_correlation(
    longitude_1,
    longitude_2,
    ellipse_longitude,
    latitude_1,
    latitude_2,
    ellipse_latitude,
    dates_1,
    dates_2,
    ellipse_age,
    phi,
    pv_1=0,
    pv_2=0,
):
    """Calculates the spatial correlation between two points.

    Calculates how closely correlated two points are in space and time

    Used to belong in "find_besthist", but was refactored and removed
    to its own file for neatness.

    Can be done with or without potential vorticity

    Parameters
    ----------
    longitude_1: longitude of point 1
    longitude_2: longitude if point 2
    ellipse_longitude: longitudinal size of ellipse
    latitude_1: latitude of point 1
    latitude_2: latitude of point 2
    ellipse_latitude: latitudinal size of ellipse
    dates_1: dates of data for point 1
    dates_2: dates of data for point 2
    ellipse_age: age of data wanted in ellipse
    phi: potential gradient
    pv_1: potential vorticity of point 1
    pv_2: potential vorticity of point 2

    Returns:
    -------
    Spatial correlation between points
    """
    pv_correlation = 0
    correlation = (
        ((longitude_1 - longitude_2) ** 2 / (ellipse_longitude**2))
        + ((latitude_1 - latitude_2) ** 2 / (ellipse_latitude**2))
        + ((dates_1 - dates_2) ** 2 / (ellipse_age**2))
    )

    if pv_1 != 0 and pv_2 != 0:
        pv_correlation = ((pv_2 - pv_1) / math.sqrt(pv_2**2 + pv_1**2) / phi) ** 2

    return correlation + pv_correlation


def sorter(msites, sites):
    """Function to fill out arrays for the piecewise-fit

    Used to find the interval in a linear fit.

    Flag points as inside our outside boundaries msites[0] to msites[1]

    Parameters
    ----------
    msites: boundaries
    sites: points

    Returns:
    -------
    Array containing values defining points as inside or outside boundaries
    """
    # put arrays together
    all_sites = np.concatenate((msites.flatten(), sites.flatten()))

    index = np.argsort(all_sites, kind="stable")

    pointer = np.argwhere(index > msites.__len__() - 1).flatten() - np.arange(1, sites.__len__() + 1)

    return pointer


def wrap_longitude(grid_long):
    """Allows longitude (0-360) to wrap beyond the 360 mark, for mapping purposes.

    Makes sure that, if the longitude is near the boundary (0 or 360) that we
    wrap the values beyond
    360 so it appears nicely on a map

    This is a refactor between get_region_data and get_region_hist_locations to
    avoid duplicate code

    Parameters
    ----------
    grid_long: array of longitude values

    Returns:
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
    """Calculates barotropic potential vorticity (pv)

    Calculates the potential vorticity for a given latitude and z

    Used to belong in "find_besthist", but was refactored and removed
    to its own file for neatness.

    Parameters
    ----------
    lat: latitude

    Returns:
    -------
    z_value: depth
    """
    earth_angular_velocity = 2 * 7.292 * 10**-5
    lat_radians = lat * math.pi / 180

    p_v = (earth_angular_velocity * math.sin(lat_radians)) / z_value

    if p_v == 0:
        p_v = 1 * 10**-5

    return p_v


def cal2dec(pa_month, pa_day, pa_hour=0, pa_minute=0):
    """Converts a calendar date (month, day, hour, minute) to a decimal date (float)

    Parameters
    ----------
    pa_month: Month in the year (where 0 is Janurary and 11 is Decemeber)
    pa_day: Day in the month
    pa_hour: Hour in the day
    pa_minute: Minute in the hour

    Returns:
    -------
    decimalised version of the date
    """
    ln_cumulative_months = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

    try:
        dec_date = float(ln_cumulative_months[pa_month] + pa_day - 1 + pa_hour / 24 + pa_minute / 60 / 24)
        if dec_date > 366:
            raise ValueError("Day is out of scope of the year") from None
        return dec_date

    except IndexError:
        raise ValueError("Month is out of scope") from None


def change_dates(cal_dates):
    """Seperates date into year, month, day, hour, minute

    Converts calendar date from one input into year, month, day, hour, minute.
    Passes these variables to cal2dec to get a decimalised date and returns it

    N.B. Change to code on the xx/06/2013: Speed up this function - C Cabanes

    Parameters
    ----------
    cal_dates: array of date as one continuous int
    eg (01/02/1974, hour 3, minute 5, second 44 = 19740102030544)

    Returns:
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

            if 13 > month > 0 and 32 > day > 0:
                day = year + (cal2dec(month - 1, day, hour, minute) / 365)
                dec_dates[i] = day

        except ValueError:
            print("Date is incorrect length or format")
            continue

    return dec_dates
