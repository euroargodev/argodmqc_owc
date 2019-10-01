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
from scipy.interpolate import griddata, interp2d, LinearNDInterpolator, RectBivariateSpline, interpn, NearestNDInterpolator
from scipy.ndimage import map_coordinates
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

def lookupNearest(x, y, data, x0, y0):
    xi = np.abs(x-x0).argmin()
    yi = np.abs(y-y0).argmin()
    return data[yi, xi]

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
    la_x = np.arange(-5, 366, 10, int)  # Goes to 366 to ensure inclusion of 365
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
    la_y = np.arange(85, -86, -10)
    la_lookup_y = np.full((38, 18), la_y, dtype=int).transpose()
    vector_y = np.concatenate(la_lookup_y.transpose(), axis=0)
    """
    Create an 18x36 matrix of
    [631  1  19  37  ...  631  1
     632  2  20  38  ...  632  2
     ...                      ...
     648 18  36  54  ...  648  18]
     """
    la_lookup_no = np.full((1, 648), np.arange(1, 649), dtype=int).reshape(36, 18)
    la_lookup_no = np.insert(la_lookup_no, 0, la_lookup_no[la_lookup_no.shape[0] - 1]).reshape(37, 18)
    la_lookup_no = np.insert(la_lookup_no, 666, la_lookup_no[1]).reshape(38, 18).transpose()

    ln_x1 = pn_float_long + .01
    ln_x2 = pn_float_long + 10.01
    ln_x3 = pn_float_long - 9.99
    ln_x4 = pn_float_long + 20.01
    ln_x5 = pn_float_long - 19.99

    ln_y1 = pn_float_lat + .01
    ln_y2 = pn_float_lat + 10.01
    ln_y3 = pn_float_lat - 9.99
    ln_y4 = pn_float_lat + 20.01
    ln_y5 = pn_float_lat - 19.99

    # wrap longitudinal values
    if ln_x3 < 0:
        ln_x3 += 360

    if ln_x5 < 0:
        ln_x5 += 360

    if ln_x1 >= 360:
        ln_x1 -= 360

    if ln_x2 >= 360:
        ln_x2 -= 360

    if ln_x4 >= 360:
        ln_x4 -= 360

    print(la_lookup_x.shape)
    print(la_lookup_y)

    test_x = la_lookup_x[0]
    test_y = la_lookup_y.transpose()[0][::-1]

    test = RectBivariateSpline(test_x, test_y, la_lookup_no.transpose(), kx=1, ky=1, s=0)
    test2 = interp2d(test_x, test_y, la_lookup_no, kind='linear')
    print(np.linalg.matrix_rank((test_x, test_y)))
    #test3 = map_coordinates(la_lookup_no, [ln_x1, ln_y1], [test_x, test_y], mode='nearest', order=1)
    #test4 = interpn(la_lookup_x, la_lookup_no.transpose(), (ln_x1, ln_y1), method='nearest')
    #test5 = NearestNDInterpolator((test_x, test_y), la_lookup_no.transpose())
    test6 = lookupNearest(test_x, test_y, la_lookup_no, ln_x1, ln_y1)
    print(test6)
    print(test2(ln_x1, ln_y1))
    print(test.ev(ln_x1, ln_y1))
    #print(test3)


wmo_boxes = scipy.loadmat('../../data/constants/wmo_boxes.mat')
longitude = 57.1794
latitude = -59.1868
find_25boxes(60, -40, wmo_boxes)


