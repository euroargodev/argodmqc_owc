"""
-----break point fit-----

Written by: Ned Campion
When: xx/09/2005
Converted to python by: Edward Small
When: 01/06/2020

Routine to get least squares fit for p[iecewise linear fit with break points
at prescribed points

y = A(1) + A(2)*x+eps                     for x(1) <= x <= b(1)
y = A(1) + A(2)*b(1) + A(3)*(x-b(1))+eps  for b(1) <= x <= b(2)
...                                            ...
y = A(1) + A(2)*(b(2)-b(1)) + ...
                       A(m+2)*(x-b(m)+eps for b(m) <= x <= x(n)
where x = vector of observed indepedent variables [n]
      y = vector of observed dependent variables  [n]
      W_i = inverse of weights for fit [n]
            if W_i is empty, then equal weighting for each point
      b = vector of break points as values of x   [m]

      A = fitting coefficients                    [m+1]
      yg = fitted estimate of y

For information on how to use this file, check the README at either:

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import numpy as np
import scipy.linalg
from ow_calibration.sorter.sorter import sorter


# pylint: disable=too-many-locals
def brk_pt_fit(x_obvs, y_obvs, w_i, breaks=None):
    """
    Get least-squares estimates for piecewise linear fit with breakpoints at prescribed points
    :param x_obvs: indendent variables
    :param y_obvs: dependent variables
    :param w_i: inverse of weights for fit
    :param breaks: vector of break points
    :return: Matrix relating observations to the fit parameters of the linear fit
    """

    if breaks is None:
        breaks = []

    b_length = breaks.__len__()
    x_length = x_obvs.__len__()
    y_length = y_obvs.__len__()

    # shake we have the same number of dependent and independent variables
    if x_length != y_length:
        residual = y_obvs
        fit_param = np.zeros((b_length + 2, 1))
        print("ERROR in brk_pt_fit: input vectors for brk_pt_fit did not match")
        return fit_param, residual

    # check inputs are flat
    x_obvs = x_obvs.flatten()
    y_obvs = y_obvs.flatten()

    # make the first point the first break point
    btem = np.concatenate(([x_obvs[0]], breaks))

    # form matrix
    # include intercept as well as the trends between each point

    trends = np.zeros((x_length, b_length + 2))
    trends[:, 0] = np.ones(x_length)
    ixb = sorter(btem, x_obvs)

    for j in range(b_length + 1):
        ixb_j = np.argwhere(ixb == j)  # point to x values greater than break point
        trends[ixb_j, j + 1] = x_obvs[ixb_j] - btem[j]
        ixb_g_j = np.argwhere(ixb > j)  # point to values less than the

        if ixb_g_j.__len__() > 0:
            trends[ixb_g_j, j + 1] = btem[j + 1] - btem[j]

    # Get least squares estimate. Use weights, if we have them
    if w_i.__len__() > 0:
        ls_est = np.dot(np.dot(trends.T, w_i), trends)

    else:
        ls_est = np.dot(trends.T, trends)

    if np.linalg.det(ls_est) == 0:
        fit_param = np.zeros((b_length + 2, 1))
        residual = y_obvs
        print("ERROR in brk_pt_fit: DET(A) == 0")
        return fit_param, residual

    # calculate fit parameters
    if w_i.__len__() > 0:
        fit_param = np.dot(np.dot(scipy.linalg.solve(ls_est, trends.T), w_i), y_obvs)

    else:
        fit_param = scipy.linalg.solve(ls_est, np.dot(trends.T, y_obvs))


    # calculate fit estimate
    residual = y_obvs - np.dot(trends, fit_param)

    return fit_param, residual
