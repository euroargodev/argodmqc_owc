"""
% INPUT:
% x,y = observations
% n_err = rms a priori error estimate for each data point
% lvcov = matrix expressing vertical and horizontal covariances (SEE WJO)
%
% Additional parameters that can be passed to the routine.
%        The data is passed as a pair of arguments, the name of the
%        variable and its value
%
% initial_breaks    1st guess for break points
%                   [Default equally spaced between 2nd and next to last
%                   points]
% max_no_breaks     maximum number of break points to be evaluated.
%                   Valid range includes -1 and 0, where
%                   -1 = mean offset only
%                    0 = mean offset and linear trend
%                   [Default = 4];
% number_breaks     specify one or more number of break points to be
%                   fitted
%                    0 = linear fit with no break points
%                   -1 = mean offset
% breaks            evaluate fit for a specified number of break points,
%
%                   if number_breaks > length(breaks), then we are
%                   specifying some breaks, but will compute the
%                   other breaks.
% nloops            specify the number of loops to compute error
%                   of piecewise fit
%                   [Default = 200]
%
% OUTPUT:
% xfit              profile number (unique of x)
% condslope         fit estimate for each profile
% condslope_err     estimated rms error for each profile
% time_deriv        estimated change per profile
% time_deriv_err    estimated rms error in change per profile
% sta_mean          mean difference between estimate and actual values
%                   averaged over each profile
% sta_rms           rms difference between estimates and actual values
%                   averaged over each profile
% NDF               The effective number of independent observations when
%                   the off-diagonal coariance (lvcov) is taken into
%                   account
%
%------------------------------------------------------------------------
%

-----Fit Condition-----

Written by: Ned Campion
When: xx/11/2007
Converted to python by: Edward Small
When: 14/06/2020

To decide which fit is optimal, we will use the small sample variation of the
Akaike Information Criterion to choose the fit.   Having chosen the
number of parameters, we then use the F-test to see if the reduction in
variance relative to the original variance is statistically significant.

We have also used the correlation matrix for the horizontal and vertical
scales to estimate the number of effective degrees of freedom for the
fits and for estimating the uncertainties of the fits
This function implements a non-linear fit of a piecewise linear fit.  The
methodology is described in:
Jones, R.H. and I. Dey, 1995, Determining one or more change points.
Chemistry and Physics of Lipids, 76, 1-6.

Cecile Cabanes, 2017: force the fit to an offset only if NDF <13. and display
a warning : to track change see change config 129 """

import numpy as np


def fit_cond(x, y, n_err, lvcov, param, br):
    """
    Get optimal fit
    :param x: observations
    :param y: observations
    :param n_err: error estimate for each observation
    :param lvcov: covariance matrix
    :param param: parameter for our breaks
    :param br: break points
    :return:
            xfit              profile number (unique of x)
            condslope         fit estimate for each profile
            condslope_err     estimated rms error for each profile
            time_deriv        estimated change per profile
            time_deriv_err    estimated rms error in change per profile
            sta_mean          mean difference between estimate and actual values
                              averaged over each profile
            sta_rms           rms difference between estimates and actual values
                              averaged over each profile
                              the off-diagonal coariance (lvcov) is taken into
                              account
    """

    # define global variables needed for ine fitting
    global A, breaks, nbr1, ubrk_g, xf, yf, w_i, xblim

    # Set up some default values
    max_brk_dflt = 4
    max_brk_in = []
    max_break = []
    nbr1 = -1
    brk_init = []  # guesses for break point
    setbreaks = 0

    nloops = 200  # number of loops to fit profile error

    # parameters for optimisation
    max_fun_evals = 1000
    tol_fun = 1 * 10 ** (-6)

    # form x and y variables for the fit
    xfit = np.unique(x)
    nfit = xfit.__len__()

    # exclude bad points before the fit. Need to reformat matrices so that they match matlab
    x = x.T.flatten()
    y = y.T.flatten()
    n_err = n_err.T.flatten()

    good = np.argwhere(np.isfinite(y) & np.isfinite(x))

    x = x[good].flatten()
    y = y[good].flatten()
    n_err = n_err[good].flatten()

    temp_lvcov = np.empty((good.__len__(), good.__len__()))

    itx = 0
    ity = 0
    for i in good:
        for j in good:
            temp_lvcov[ity, itx] = lvcov[i, j]
            itx += 1
        itx = 0
        ity += 1

    lvcov = temp_lvcov

    npts = x.__len__()

    # check that we actually have some good data

    if npts == 0:
        print("Failed to find any good data")

        condslope = np.nan
        condslope_err = np.nan
        time_deriv = np.nan
        time_deriv_err = np.nan
        sta_mean = np.nan
        sta_rms = np.nan
        NDF = []
        fit_coef = []
        fit_breaks = []

        return (xfit, condslope, condslope_err, time_deriv, time_deriv_err, sta_mean,
                sta_rms, NDF, fit_coef, fit_breaks)

    #condition the series so that the fit is well behaved
    # sort by the independent variable

    x = np.sort(x)
    sorted_index = np.argsort(x, kind='stable')
    y = y[sorted_index]
    n_err = n_err[sorted_index]

    # scale x from -1 to 1

    x_0 = (x[npts - 1] + x[0])/2

    if x[0] != x[npts - 1]:
        x_scale = (x[npts - 1] - x[0])/2

    else:
        x_scale = 1

    # remove the mean of y and scale by the standard deviation

    y_0 = np.mean(y)
    y_scale = np.std(y)

    if y_scale == 0:
        y_scale = 1

    # calculate x and y used for fitting routine
    xf = (x - x_0)/x_scale
    yf = (y - y_0)/y_scale
    n_err = n_err/y_scale
    xfit = (xfit - x_0)/x_scale

    # get profile times that will be used as independent variables to get
    # error statstics. xp could be different to xfit if there is a profile
    # with no good data

    x_unique, index_unique = np.unique(xf, return_index=True)
    n_prof = x_unique.__len__()

    # convert errors from rms to variance
    err_var = (n_err) ** 2

    # weights for weighted least squares fit

    w_i = np.diag(np.mean(err_var) / err_var)

    # use correlation matrix to compute degrees of freedom

    ndf = np.sum(np.ones((npts, 1)) / (np.dot(lvcov, np.ones((npts,1)))))
    






