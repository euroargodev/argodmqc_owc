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
the range [5W 5E] are out-of-bound with matlab version >=R2012b

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
import scipy.io as scipy


def find_25boxes(pn_float_long, pn_float_lat, pa_wmo_boxes):

    """
    :param pn_float_long: float longitude, float
    :param pn_float_lat: float latitude, float
    :param pa_wmo_boxes: wmo boxes (explained above), matrix
    :return: (explained above), 25x4 matrix
    """

    """
    Create matrix with 18 rows of
    [-5, 5, 15, ..., 355, 365]
    """
    la_x = np.arange(-5, 366, 10, int) #Goes to 366 to ensure inclusion of 365
    la_lookup_x = np.full((18, 38), la_x, dtype=int)

    """
    Create matrix with 18 columns of
    [
    85,
    75,
    65,
    ...,
    -75,
    -85
    ]
    """
    la_y = np.arange(85,-86, -10)
    la_lookup_y = np.full((38, 18), la_y, dtype=int).transpose()
    vector_y = np.concatenate(la_lookup_y.transpose(), axis=0)

    """
    Create an 18x36 matrix of
    [631  1  19  37  ...  631  1
     632  2  20  38  ...  632  2
     ...                      ...
     648 18  36  54  ...  648  18]
     """
    la_lookup_no = np.full((1, 648), np.arange(1,649), dtype=int).reshape(36, 18)
    la_lookup_no = np.insert(la_lookup_no, 0, la_lookup_no[la_lookup_no.shape[0]-1]).reshape(37, 18)
    la_lookup_no = np.insert(la_lookup_no, 666, la_lookup_no[1]).reshape(38, 18)
    la_lookup_no = la_lookup_no.transpose()


wmo_boxes = scipy.loadmat('../../data/constants/wmo_boxes.mat')
print(wmo_boxes)
longitude = 57.1794
latitude = -59.1868
find_25boxes(longitude,latitude,wmo_boxes)