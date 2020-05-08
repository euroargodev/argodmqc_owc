"""
-----Find 25 Boxes-----

Written by: A. Wong
When: 16/08/2004
Converted to python by: Edward Small
When: 24/09/2019

Finds the 5 x 5 = 25 WMO boxes with the float profile in the centre
The WMO box numbers, between 90N and 90S are stored in /data/constants/wmo_boxes.mat
The structure of the matrix is as such:

Column 1 - box number
Column 2 - Do we have CTD data (1 = yes, 0 = no)
Column 3 - Do we have bottle data (1 = yes, 0 = no)
Column 4 - do we have Argo data (1 = yes, 0 = no)

N.B. Change to code on the xx/11/2014: extend la_x so interp2 does not think that longitudes in
the range [5W 5E] are out-of-bound with matlab version >=R2012b - C Cabanes

N.B. Change during conversion to python on the 01/10/2019. Struggled to find an interpolation
function that exactly mirrored Matlab's, so I wrote my own - Edward Small

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np


def nearest_neighbour(x_axis, y_axis, table, x_input, y_input):
    """
    :param x_axis: x-axis values
    :param y_axis: y-axis values
    :param table: grid data
    :param x_input: data to interpolate
    :param y_input: data to interpolate
    :return: nearest neighbour
    """
    x_output = np.abs(x_axis - x_input).argmin()
    y_output = np.abs(y_axis - y_input).argmin()
    return table[y_output, x_output]


def find_25boxes(pn_float_long, pn_float_lat, pa_wmo_boxes):
    """
    :param pn_float_long: float longitude, float
    :param pn_float_lat: float latitude, float
    :param pa_wmo_boxes: wmo boxes (explained above), data frame
    :return: (explained above), 25x4 matrix

    First, we need to create a look-up table in the form of

        | -5    5   15  ...  355  365
    ----|----------------------------
     85 | 631   1   19  ...  631   1
     75 | 632   2   20  ...  632   2
     65 | 633   3   21  ...  633   3
     ...|...   ...  ... ...  ...  ...
     -85| 648   18  36  ...  648   18

     We do this using 3 matrices:

     - A 1-D matrix for the x axis (la_lookup_x)
     - A 1-D matrix for the y axis (la_lookup_y)
     - A 2-D matrix for the grid data (la_lookup_no)
    """

    la_lookup_x = np.arange(-5, 366, 10, int)

    la_lookup_y = np.arange(85, -86, -10).transpose()

    la_lookup_no = np.full((1, 648), np.arange(1, 649), dtype=int).reshape(36, 18)
    la_lookup_no = np.insert(
        la_lookup_no, 0, la_lookup_no[la_lookup_no.shape[0] - 1]
    ).reshape(37, 18)
    la_lookup_no = np.insert(la_lookup_no, 666, la_lookup_no[1]).reshape(38, 18).transpose()

    # Set up longitudinal and latitudinal values
    ln_x = []
    ln_y = []
    ln_x.append(pn_float_long + .01)
    ln_x.append(pn_float_long + 10.01)
    ln_x.append(pn_float_long - 9.99)
    ln_x.append(pn_float_long + 20.01)
    ln_x.append(pn_float_long - 19.99)

    ln_y.append(pn_float_lat + .01)
    ln_y.append(pn_float_lat + 10.01)
    ln_y.append(pn_float_lat - 9.99)
    ln_y.append(pn_float_lat + 20.01)
    ln_y.append(pn_float_lat - 19.99)

    # wrap longitudinal values
    if ln_x[2] < 0:
        ln_x[2] += 360

    if ln_x[4] < 0:
        ln_x[4] += 360

    if ln_x[0] >= 360:
        ln_x[0] -= 360

    if ln_x[1] >= 360:
        ln_x[1] -= 360

    if ln_x[3] >= 360:
        ln_x[3] -= 360

    if not np.isnan(pn_float_long) and not np.isnan(pn_float_lat):
        ln_i = []
        for i in range(0, 5):
            for j in range(0, 5):
                ln_i.append(
                    nearest_neighbour(la_lookup_x, la_lookup_y, la_lookup_no, ln_x[j], ln_y[i])
                )

    else:
        ln_i = np.full(25, np.nan)

    pa_wmo_numbers = np.full((25, 4), np.nan)
    for i in range(0, 25):
        if not np.isnan(ln_i[i]):
            pa_wmo_numbers[i] = pa_wmo_boxes.get('la_wmo_boxes')[ln_i[i] - 1]

    return pa_wmo_numbers
