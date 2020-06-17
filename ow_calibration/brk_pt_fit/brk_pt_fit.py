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

def brk_pt_fit(x, y, w_i, b=[]):
    """
    Get least-squares estimates for piecewise linear fit with breakpoints at prescribed points
    :param x: indendent variables
    :param y: dependent variables
    :param w_i: inverse of weights for fit
    :param b: vector of break points
    :return: Matrix relating observations to the fit parameters of the linear fit
    """