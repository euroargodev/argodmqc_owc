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
import scipy.optimize as scipy
from ow_calibration.brk_pt_fit.brk_pt_fit import brk_pt_fit
from ow_calibration.sorter.sorter import sorter


def nlbpfun(ubrk_i):
    """
    find residual
    :param ubrk_i: input
    :return: residual
    """

    global A, breaks, nbr1, ubrk_g, xf, yf, w_i, xblim

    if nbr1 > 1:
        ubrk = ubrk_g[0:nbr1 - 1]
        for i in range(nbr1, ubrk_g.__len__()):
            ubrk.append(ubrk_i)

    else:
        ubrk = ubrk_i

    m_b = ubrk.__len__()
    fnumer = np.zeros(ubrk.shape)
    fnumer[0] = np.exp(ubrk[0])

    for i in range(1, m_b):
        fnumer[i] = fnumer[i - 1] + np.exp(ubrk(i))

    fdenom = 1 + fnumer[m_b - 1]

    ftem = (xblim[1] - xblim[0]) / fdenom

    breaks = xblim[0] + ftem * fnumer

    if np.argwhere(np.diff(breaks) == 0).__len__() > 0:
        difference = np.argwhere(np.diff(breaks) == 0)
        breaks[difference + 1] = breaks[difference + 1] + 0.00001

    A, residual = brk_pt_fit(xf, yf, w_i, breaks)

    return residual


def fit_cond(x, y, n_err, lvcov, *args):
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
    tol = 1e-06
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

    # condition the series so that the fit is well behaved
    # sort by the independent variable

    x = np.sort(x)
    sorted_index = np.argsort(x, kind='stable')
    y = y[sorted_index]
    n_err = n_err[sorted_index]

    # scale x from -1 to 1

    x_0 = (x[npts - 1] + x[0]) / 2

    if x[0] != x[npts - 1]:
        x_scale = (x[npts - 1] - x[0]) / 2

    else:
        x_scale = 1

    # remove the mean of y and scale by the standard deviation

    y_0 = np.mean(y)
    y_scale = np.std(y)

    if y_scale == 0:
        y_scale = 1

    # calculate x and y used for fitting routine
    xf = (x - x_0) / x_scale
    yf = (y - y_0) / y_scale
    n_err = n_err / y_scale
    xfit = (xfit - x_0) / x_scale

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

    ndf = np.sum(np.ones((npts, 1)) / (np.dot(lvcov, np.ones((npts, 1)))))

    # calculate the residual sum of squares for the initial series

    rss_0 = np.sum((yf ** 2) / err_var)

    if x_unique.__len__() > 3:
        # find 2nd and 2nd to last profile and use them as limits for the break points
        xblim = [x_unique[1], x_unique[n_prof - 2]]

    else:
        # too few profiles
        xblim = [1, 1]

    no_args = args.__len__()

    if np.remainder(no_args, 2) != 0:
        raise ValueError("FIT_COND ERROR - inputs are incorrect")

    if no_args > 0:
        for n in range(int(no_args / 2)):
            parm = args[2 * (n - 1)]
            value = args[(n * 2) + 1]

            if type(parm) != str:
                raise ValueError("FIT_COND ERROR - inputs are incorrect")

            param = str.lower(parm)

            if param == 'initial_breaks':
                # initial guess for breakpoints
                brk_init = value

                # rescale
                brk_init = (brk_init - x_0) / x_scale
                brk_init = (brk_init - xblim[0]) / np.diff(xblim)

            elif param == 'max_no_breaks':
                if value.__len__() > 0:
                    max_brk_in = value[0]
                    nbr1 = -1

            elif param == 'number_breaks':
                pbrk = value[0]
                nbr1 = pbrk
                max_brk_in = pbrk

            elif param == nloops:
                nloops = value[0]

            elif param == 'breaks':
                if value.__len__() > 0:
                    breaks = value[0]
                    breaks = (breaks - x_0) / x_scale
                    nbr = breaks.__len__()
                    setbreaks = 1

            else:
                raise ValueError("Paramater " + param + " not found in parameter list")

    # intialise variable for search over number of break points
    b_pts = np.ones((max_brk_in, max_brk_in + 1)) * np.nan
    b_A = np.ones((max_brk_in + 2, max_brk_in + 1)) * np.nan
    rss = np.ones((1, max_brk_in + 2)) * np.nan
    aic = np.ones((1, max_brk_in + 2)) * np.nan

    # check to see if we have set break points
    if setbreaks:
        if max_brk_in == 0:
            max_brk_in = nbr
            nbr1 = nbr

        # we have fixed break points
        elif max_brk_in > nbr:
            nbr1 = nbr + 1
            A, residual = brk_pt_fit(xf, yf, w_i, breaks)
            b_pts[0:nbr, nbr + 1] = breaks.T
            b_A[0:nbr + 2, nbr + 2] = A[0:nbr + 2]
            rss[0, nbr + 2] = np.sum(residual ** 2 / err_var)
            no_param = 2 * (nbr + 1)
            aic[0, nbr + 2] = ndf * np.log(rss[0, nbr + 2] / npts) + \
                              ndf * (ndf + A) / (ndf - A - 2)

        # we have the same number of specified breaks
        else:
            nbr1 = nbr

        max_brk = max_brk_in
        pbrk = np.arange(nbr1, max_brk + 1)

    else:
        # no break points set
        if isinstance(max_brk_in, list):
            max_brk_in = max_brk_dflt

        max_brk = max_brk_in
        pbrk = np.arange(nbr1, max_brk + 1)

    if ndf < 2 * (max_brk + 2) + 1:
        if ndf > 2 * (nbr1 + 2) + 1:
            pbrk = np.arange(nbr1, np.floor((ndf - 1) / 2 - 2) + 1)
            print("WARNING: only have " + str(ndf) + " degrees of freedom")
            print("Maximum breakpoints to be tried: " + str(np.max(pbrk)))

        else:
            if setbreaks == 1:
                pbrk = nbr
                max_brk = nbr
                nbr1 = nbr
                print("WARNING: Only have " + str(ndf) + " degrees of freedom")
                print("Estimate fit with fixed breakpoints")

            else:
                pbrk = -1
                print("WARNING: Only have " + str(ndf) + " degrees of freedom")
                print("Estimate offset only")

    nbr = 1
    if nbr == -1:
        # offset only
        # since this is an error weighted average, yx won't necessarily be 0
        ones_column = np.ones((npts, 1))
        b_A[0, 0] = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(ones_column.T, w_i),
                                                              ones_column)),
                                         ones_column.T),
                                  w_i),
                           yf)

        residual = yf - (ones_column * b_A[0, 0]).flatten()
        rss[0, 0] = np.sum(residual ** 2 / err_var)
        aic[0, 0] = ndf * np.log(rss[0, 0] / npts) + ndf * (ndf + 1) / (ndf - 3)

    elif nbr == 0:
        # linear fit, no break points
        A, residual = brk_pt_fit(xf, yf, w_i)
        b_A[0:2, 1] = A[0:2]
        rss[0, 1] = np.sum(residual ** 2 / err_var)
        aic[0, 1] = ndf * np.log(rss[0, 1] / npts) + ndf * (ndf + 2) / (ndf - 4)

    else:
        nbr2 = brk_init.__len__()

        # Check if there are enough initial guesses
        if nbr2 >= nbr:
            if brk_init.shape[0] > brk_init.shape[1]:
                brk_init = brk_init.T

            b_guess = brk_init[0:nbr]

        # first guess for breaks as evenly distributed between 2nd and 2nd to last point
        else:
            b_guess = -1 + 2 * np.arange(1, nbr + 1) / (nbr + 1)

        b_g = np.concatenate(([-1], b_guess))
        ubrk_g = []

        for n in range(nbr):
            ubrk_g.append(np.log((b_g[n + 1] - b_g[n]) / (1 - b_g[nbr])))

        if setbreaks:
            # break points are already set
            if nbr1 == max_brk:
                A, residual = brk_pt_fit(xf, yf, w_i, breaks)

            # fit over limited number of breaks
            else:
                optim = scipy.least_squares(nlbpfun, ubrk_g[nbr1:nbr],
                                            method='lm', ftol=tol)
                ubrk = optim['x'][0]
                residual = optim['fun']

                ubrk = np.concatenate((ubrk_g[0:nbr1 - 1], ubrk))
        # get non-linear least squares for break points
        else:
            optim = scipy.least_squares(nlbpfun, ubrk_g,
                                        method='lm', ftol=tol)
            ubrk = optim['x'][0]
            residual = optim['fun']

        b_pts[0:nbr, nbr] = breaks.T
        b_A[0:nbr + 2, nbr + 1] = A[0:nbr + 2]
        rss[0, nbr + 1] = np.sum(residual ** 2 / err_var)
        p = 2 * (nbr + 1)
        aic[0, nbr + 1] = ndf * np.log(rss[0, nbr + 1] / npts) + ndf * (ndf + p) / (ndf - p - 2)

    if setbreaks and nbr1 == max_brk:
        best = pbrk + 2

    # decide which fit to use (offset, linear, piecewise)
    else:
        if nbr1 > 1:
            pbrk = np.arange((nbr1 - 1), max_brk)

        good = pbrk + 1
        aic_b = np.min(aic[0, good])
        best = np.argmin(aic[0, good])

        if isinstance(good, np.ndarray):
            best = good[best] + 1
        else:
            best = good + 1

    # best fit includes break points
    if best > 2:
        no_breaks = (best - 1) * 2

    else:
        no_breaks = best

    if setbreaks & nbr1 == max_brk:
        comment = "Fit evaluated "

    else:
        comment = "Best model found with "

    if best > 2:
        comment = comment + str(best - 2) + " break points"

    elif best == 2:
        comment = comment + "linear fit"

    else:
        comment = comment + "offset value only"

    print(comment)

    if best > 2:
        breaks = b_pts[np.arange(0, best - 2), best - 2].T

    else:
        breaks = []

    A = b_A[0:best, best]
    btem = np.concatenate(([xf[0]], breaks))
    E = np.zeros((npts, best))
    E[:, 0] = np.ones(npts).T
    ixb = sorter(btem, xf)

    if best > 1:
        for j in range(best - 1):
            ib = np.argwhere(ixb == j)
            E[ib, j + 1] = xf[ib] - btem[j]
            ii = np.argwhere(ixb > j)

            if ii.__len__() > 0:
                E[ii, j + 1] = btem[j + 1] - btem[j]

    # get uncertainnties in fit parameters
    B_i = (np.dot(np.dot(E.T, w_i), E)) ** -1
    P = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(B_i, E.T),
                                           w_i), np.diag(err_var)),
                             w_i),
                      E),
               B_i)
    P_1 = np.diag(P)
    P = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(B_i, E.T),
                                           w_i),
                                    np.diag(err_var)),
                                    lvcov),
                             w_i),
                      E),
               B_i)
    P_2 = np.diag(P)

    # reduce matrix to have only one value per profile
    btem = np.concatenate[xfit[0], breaks]
    print(btem)
