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

from ow_calibration.sorter.sorter import sorter
import numpy as np


def brk_pt_fit(x, y, w_i, b=[]):
    """
    Get least-squares estimates for piecewise linear fit with breakpoints at prescribed points
    :param x: indendent variables
    :param y: dependent variables
    :param w_i: inverse of weights for fit
    :param b: vector of break points
    :return: Matrix relating observations to the fit parameters of the linear fit
    """

    b_length = b.__len__()
    x_length = x.__len__()
    y_length = y.__len__()

    # shake we have the same number of dependent and independent variables
    if x_length != y_length:
        residual = 999
        a = np.zeros((b_length + 2, 1))
        print("ERROR: input vectors for brk_pt_fit did not match")
        return a, residual

    # check inputs are flat
    x = x.flatten()
    y = y.flatten()

    # make the first point the first break point
    btem = np.concatenate(([x[0]], b))

    # form matrix
    # include intercept as well as the trends between each point

    e = np.zeros((x_length, b_length + 2))
    e[:, 0] = np.ones(x_length)
    ixb = sorter(btem, x)

    for j in range(b_length + 1):
        ib = np.argwhere(ixb == j)  # point to x values greater than break point
        e[ib, j+1] = x[ib] - btem[j]
        ii = np.argwhere(ixb > j)  # point to values less than the

        if ii.__len__() > 0:
            e[ii, j+1] = btem[j+1] - btem[j]

    # Get least squares estimate. Use weights, if we have them
    if w_i.__len__() > 0:
        b = np.dot(np.dot(e.T, w_i), e)

    else:
        b = np.dot(e.T, e)
