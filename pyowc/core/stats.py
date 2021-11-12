""" Key functions to run statistics

        Parameters
        ----------

        Returns
        -------

"""
import copy
import numpy as np
from scipy import interpolate
from scipy.optimize import least_squares
from scipy import linalg
from scipy.spatial import cKDTree  # pylint: disable=no-name-in-module

from ..utilities import sorter, potential_vorticity


#pylint: disable=too-many-lines
#pylint: disable=global-variable-undefined
#pylint: disable=redefined-outer-name
#pylint: disable=invalid-name
#pylint: disable=too-many-locals
#pylint: disable=too-many-branches
#pylint: disable=too-many-statements
def fit_cond(x, y, n_err, lvcov, *args):
    """ Get optimal fit

        To decide which fit is optimal, we will use the small sample variation of the
        Akaike Information Criterion.   Having chosen the
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
        a warning : to track change see change config 129

        Parameters
        ----------
        x: observations
        y: observations
        n_err: error estimate for each observation
        lvcov: covariance matrix
        param: parameter for our breaks
        br: break points

        Returns
        -------
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
    def nlbpfun(ubrk_i):
        """ Find residual

            Parameters
            ----------
            ubrk_i: input

            Returns
            -------
            residual
        """

        # TODO Address the use of globals in this file.
        global A, breaks, nbr1, ubrk_g, xf, yf, w_i, xblim  # pylint: disable=global-variable-not-assigned

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
            fnumer[i] = fnumer[i - 1] + np.exp(ubrk[i])

        fdenom = 1 + fnumer[m_b - 1]

        ftem = (xblim[1] - xblim[0]) / fdenom

        breaks = xblim[0] + ftem * fnumer

        if np.argwhere(np.diff(breaks) == 0).__len__() > 0:
            difference = np.argwhere(np.diff(breaks) == 0)
            breaks[difference + 1] = breaks[difference + 1] + 0.00001

        A, residual = brk_pt_fit(xf, yf, w_i, breaks)

        return residual

    # define global variables needed for line fitting
    global A, breaks, nbr1, ubrk_g, xf, yf, w_i, xblim


    # Set up some default values
    tol = 1e-06
    max_brk_dflt = 4
    max_brk_in = []
    nbr1 = -1
    brk_init = [] # guesses for break point
    breaks = np.array([])
    setbreaks = 0

    nloops = 200  # number of loops to fit profile error

    # parameters for optimisation
    max_fun_evals = 1000

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

    if x_unique.__len__() > 3:
        # find 2nd and 2nd to last profile and use them as limits for the break points
        xblim = [x_unique[1], x_unique[n_prof - 2]]

    else:
        # too few profiles
        xblim = [1, 1]

    no_args = args.__len__()

    if np.remainder(no_args, 2) != 0:
        raise ValueError("FIT_COND ERROR - inputs are incorrect") from None

    if no_args > 0:
        for n in range(int(no_args / 2)):
            parm = args[2 * n]
            value = args[(n * 2) + 1]

            if not isinstance(parm, str):
                raise ValueError("FIT_COND ERROR - inputs are incorrect") from None

            param = str.lower(parm)

            if param == 'initial_breaks':
                # initial guess for breakpoints
                brk_init = value

                # rescale
                brk_init = (brk_init - x_0) / x_scale
                brk_init = (brk_init - xblim[0]) / np.diff(xblim)

            elif param == 'max_no_breaks':
                max_brk_in = value
                nbr1 = -1

            elif param == 'number_breaks':
                pbrk = value
                nbr1 = pbrk
                max_brk_in = pbrk

            elif param == 'nloops':
                nloops = value

            elif param == 'breaks':
                if value.__len__() > 0:
                    breaks = value
                    breaks = (breaks - x_0) / x_scale
                    nbr = breaks.__len__()
                    setbreaks = 1

            else:
                raise ValueError("Paramater " + param + " not found in parameter list") from None

    # intialise variable for search over number of break points
    max_brk_in = int(max_brk_in)

    if max_brk_in == -1:
        b_pts = np.ones((1, 1)) * np.nan
        b_A = np.ones((1, 1)) * np.nan

    else:
        b_pts = np.ones((max_brk_in, max_brk_in + 1)) * np.nan
        b_A = np.ones((max_brk_in + 2, max_brk_in + 2)) * np.nan

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
                              ndf * (ndf + no_param) / (ndf - no_param - 2)

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
                pbrk = np.array([nbr])
                max_brk = nbr
                nbr1 = nbr
                print("WARNING: Only have " + str(ndf) + " degrees of freedom")
                print("Estimate fit with fixed breakpoints")

            else:
                pbrk = np.array([-1])
                print("WARNING: Only have " + str(ndf) + " degrees of freedom")
                print("Estimate offset only")

    for nbr in pbrk:
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

            nbr = int(nbr)
            for n in range(nbr):
                ubrk_g.append(np.log((b_g[n + 1] - b_g[n]) / (1 - b_g[nbr])))

            if setbreaks != 0:
                # break points are already set
                if nbr1 == max_brk:
                    A, residual = brk_pt_fit(xf, yf, w_i, breaks)

                # fit over limited number of breaks
                else:
                    optim = least_squares(nlbpfun, ubrk_g[nbr1:nbr],
                                          method='lm', ftol=tol, max_nfev=max_fun_evals)
                    ubrk = optim['x'][0]
                    residual = optim['fun']

                    ubrk = np.concatenate((ubrk_g[0:nbr1 - 1], ubrk))
            # get non-linear least squares for break points
            else:
                ubrk_g = np.array(ubrk_g)
                optim = least_squares(nlbpfun, ubrk_g,
                                      method='lm', ftol=tol, max_nfev=max_fun_evals)
                ubrk = optim['x'][0]
                residual = optim['fun']

            if breaks.__len__() > 0:
                b_pts[0:nbr, nbr] = breaks.T

            else:
                b_pts[0:nbr, nbr] = np.nan

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

        good = np.array(pbrk + 1, dtype=int)
        best = np.argmin(aic[0, good])

        if isinstance(good, np.ndarray):
            best = good[best] + 1
        else:
            best = good + 1

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

    best = int(best)
    A = b_A[0:best, best - 1]
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
    B_i = np.linalg.inv((np.dot(np.dot(E.T, w_i), E)))
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
    btem = np.concatenate(([xfit[0]], breaks))
    E = np.zeros((xfit.__len__(), best))
    E[:, 0] = np.ones((xfit.__len__(), 1)).T
    ixb = sorter(btem, xfit)

    if best >= 2:
        for j in range(best - 1):
            # point to x values greater than break j
            ib = np.argwhere(ixb == j)
            E[ib, j + 1] = xfit[ib] - btem[j]
            # point to x values less than the one just below x
            ii = np.argwhere(ixb > j)

            if ii.__len__() > 0:
                E[ii, j + 1] = btem[j + 1] - btem[j]

    # fit the values
    yfit = np.dot(E, A)

    # factor to increase monte carlo error estimate, taking into account
    # that we have correlated noise
    P_3 = np.dot(E, P_2) / np.dot(E, P_1)

    real_A = copy.deepcopy(A)
    real_E = copy.deepcopy(E)
    real_breaks = copy.deepcopy(breaks)
    real_xf = copy.deepcopy(xf)
    real_yf = copy.deepcopy(yf)
    ubrk_g = []

    if best == 1:
        err = np.ones((nfit, 1)) * P_2[0]

    else:
        err = 0

        for i in range(nloops):

            yf = real_yf + n_err * np.random.randn(yf.size)

            if best == 2:
                # E for linear case is already calculated
                A, residual = brk_pt_fit(xf, yf, w_i)

            elif setbreaks:
                # E stays fixed if breaks are specified
                A, residual = brk_pt_fit(xf, yf, w_i, breaks)

            else:
                # give an initial guess as the fitted break points to speed up calculation
                nbr = real_breaks.__len__()
                b_g = np.concatenate(([-1], real_breaks))

                for n in range(nbr):
                    ubrk_g.append(np.log((b_g[n + 1] - b_g[n]) / (1 - b_g[nbr])))

                optim = least_squares(nlbpfun, ubrk_g,
                                      method='lm', ftol=tol, max_nfev=max_fun_evals)

                ubrk = optim['x'][0]
                residual = optim['fun']

                btem = np.concatenate([xfit[0]], breaks)
                E = np.zeros((xfit.__len__(), best))
                E[:, 0] = np.ones((xfit.__len__(), 1)).T
                ixb = sorter(btem, xfit)

                for j in range(best - 1):
                    # pointer to x values greater than break point j
                    ib = np.argwhere(ixb == j)
                    E[ib, j + 1] = xfit[ib] - btem[j]
                    # pointer to break points less than the one just below x
                    ii = np.argwhere(ixb > j)

                    if ii.__len__() > 0:
                        E[ii, j + 1] = btem[j + 1] - btem[j]

            err = err + (yfit - np.dot(E, A)) ** 2

        err = err / nloops

        # rescale error to reflect the decrease to the off diagonal covariances

        err = err * P_3

    A = copy.deepcopy(real_A)
    E = copy.deepcopy(real_E)
    breaks = copy.deepcopy(real_breaks)
    xf = copy.deepcopy(real_xf)
    yf = copy.deepcopy(real_yf)

    # get residual statistics for each profile

    w = np.diag(w_i)
    sta_mean = np.ones((1, nfit)) * np.nan
    sta_rms = np.ones((1, nfit)) * np.nan
    ip_1 = np.concatenate(([-1], index_unique.flatten()))

    for n in range(nfit):
        index = np.argwhere(x_unique == xfit[n])

        if index.__len__() > 0:
            w_values = w[np.arange(ip_1[index][0, 0] + 1, ip_1[index + 1][0, 0] + 1)]
            yf_values = yf[np.arange(ip_1[index][0, 0] + 1, ip_1[index + 1][0, 0] + 1)]

            sta_mean[0, n] = np.sum(w_values * yf_values) / np.sum(w_values)

            sta_rms[0, n] = np.sum(w_values * (sta_mean[0, n] - yf_values) ** 2) / np.sum(w_values)
            sta_rms[0, n] = np.sqrt(sta_rms[0, n])

    # convert back to original units

    x_unique = x_unique * x_scale + x_0
    xfit = xfit * x_scale + x_0

    # pass the coeffecients and break points back
    break_pts = []
    if breaks.__len__() > 0:
        break_pts = breaks * x_scale + x_0

    A[0] = A[0] * y_scale + y_0
    P[0] = P[0] * y_scale

    if A.__len__() > 1:
        A[1:best] = A[1:best] * y_scale / x_scale
        P[1:best] = P[1:best] * y_scale / x_scale

    yfit = yfit * y_scale + y_0
    err = np.sqrt(err) * y_scale
    yerr = err
    n_err = n_err * y_scale
    sta_mean = sta_mean * y_scale + y_0
    sta_rms = sta_rms * y_scale

    # get time derivatives and time derivatives errors
    ixb = sorter(np.concatenate(([np.min(xfit)], break_pts)), xfit)

    if best == 1:
        time_deriv = np.zeros((nfit, 1))
        time_deriv_err = np.ones((nfit, 1)) * np.nan

    elif best == 2:
        time_deriv = A[1] * np.ones((nfit, 1))
        time_deriv_err = P_2[1] * np.ones((nfit, 1))

    else:
        time_deriv = np.ones((nfit, 1)) * np.nan
        time_deriv_err = np.ones((nfit, 1)) * np.nan
        for j in range(best - 1):
            ib = np.argwhere(ixb == j)
            time_deriv[ib] = A[j + 1]
            time_deriv_err[ib] = P_2[j + 1]

    condslope = yfit.T
    condslope_err = yerr.T

    # return fit parameters
    fit_coef = A
    fit_breaks = break_pts

    return (xfit, condslope, condslope_err, time_deriv, time_deriv_err,
            sta_mean, sta_rms, ndf, fit_coef, fit_breaks)


def signal_variance(sal):
    """ Calculates signal variance

        Calculates an estimate of the signal variance at a given level using:

        (sum({di - <d>}^2))/N

        where di is a data point in d (which is a collection of salinities at a given level,
         <d> is the mean of all the data points in d, and N is the number of data points
         in d.

        Parameters
        ----------
        sal: vector of salinities at a given level

        Returns
        -------
        float estimate of the variance of the signal of the given data
    """
    sal = np.array(sal)
    non_zero_values = sal[sal != 0.0]

    if np.all(np.isnan(non_zero_values)):
        raise RuntimeError("Received no valid salinity values when calculating signal variance") from None

    return np.nanvar(non_zero_values)


def noise_variance(sal, lat, long):
    """ Calculates the variance in noise_variance of salinity at different pressures

        Finds the variance for the noise variance of each salinity measurement by comparing it
        to the noise_variance of all the other measurements. Can be thought of as the average
        difference between measured salinity and expected salinity:

        noise_variance = (sum({x - y}^2)) / 2*N where
        x is the current observation
        y is the closest observation (spatially)
        N is the number of elements

        This is because we assume that the noise variance is uncorrelated over distance, that it has
        uniform variance, and that the signal has a longer correlation distance than the
        data separation (WJO 2003, Delayed-Mode Calibration of Autonomous CTD Profiling
        Float Salinity Data by θS Climatology).

        Parameters
        ----------
        sal: m*n matrix containing m layers and n casts
        lat: vector of latitudes of each cast
        long: vector of longitudes of each cast

        Returns
        -------
        the variance in the noise_variance
    """
    locations = np.column_stack((lat, long))

    kdtree = cKDTree(locations)  # pylint: disable=not-callable

    # query the second nearest neighbour to exclude self
    distances, min_distances_indices = kdtree.query(locations, k=[2])

    if np.all(distances == 0.0):
        print("WARNING: no unique points")
        return 0.0

    # remove the last (singleton) axis which results from querying only the second nearest neighbour
    min_distances_indices = np.squeeze(min_distances_indices, axis=-1)

    sal_noise = sal - sal[min_distances_indices]
    sal_noise_var = np.mean(np.square(sal_noise[sal_noise != 0.0])) / 2

    return sal_noise_var


#pylint: disable=too-many-arguments
def covar_xyt_pv(points1, points2, lat, long, age, phi, map_pv_use):
    """ Calculates how "close" two sets of points are to each other,

        Calculates how "close" two sets of points are to each other, taking into account
        space, time and (if wanted) potential vorticity. The closer two values are to
        each other, the closer their value will be to 1. Points that differ greatly will be
        nearer to 0.

        Calculates covariance of each point against every other point using the
        Squared Exponential (SE) covariance function:

        SE(x,y) = exp(-(x-y)^2/2l) where (x-y) is the difference between points (could be distance,
        time, etc), and l is the characteristic length scale (how close points have to be
        to influence each other significantly).

        Parameters
        ----------
        points1: m*4 model matrix containing the latitude, longitude, date, and depth of each data point
        points2: n*4 data matrix containing the latitude, longitude, date, and depth of each data point
        lat: float, the characteristic latitude
        long: float, the characteristic longitude
        age: float, the characteristic time scale
        phi: float, the characteristic cross-isobaric scale (for depth dependence)
        map_pv_use: int, flag for using vorticity (1=include)

        Returns
        -------
        m*n matrix containing the covariance of each point
    """
    points1 = np.atleast_2d(points1)
    points2 = np.atleast_2d(points2)

    long_covar = ((points1[:, np.newaxis, 0] - points2[np.newaxis, :, 0]) / long)**2
    lat_covar = ((points1[:, np.newaxis, 1] - points2[np.newaxis, :, 1]) / lat)**2
    age_covar = 0.0
    p_v_covar = 0.0

    # pylint: disable=fixme
    if age != 0:
        age_covar = ((points1[:, np.newaxis, 2] - points2[np.newaxis, :, 2]) / age)**2

    if map_pv_use == 1:
        # define a vectorized function to calculation potential vorticity
        potential_vorticity_vec = np.vectorize(potential_vorticity)

        pv_lat1 = potential_vorticity_vec(points1[:, 0], points1[:, 3])
        pv_lat2 = potential_vorticity_vec(points2[:, 0], points2[:, 3])

        p_v_covar = ((pv_lat1[:, np.newaxis] - pv_lat2[np.newaxis, :]) /
                     np.sqrt(pv_lat1[:, np.newaxis]**2 + pv_lat2[np.newaxis, :]**2) / phi)**2

    return np.exp(-(lat_covar + long_covar + age_covar + p_v_covar))


#pylint: disable=too-many-locals
def build_cov(ptmp, coord_float, config):
    """ Build the covariance matrix

        Function builds a square covariance matrix that has n*n tiles, and each
        tile is of size m*m

        The vertical covariance is the building tile. It contains a 1 down its diagonal,
        which decreases exponentially in the off-diagonals (representing the vertical covariance
        between water masses).

        We assume that each profile is independent from the other profiles

        N.B. Change to code on the xx/11/2007: Undefined - Breck Owens

        N.B. Change to code on the xx/06/2013: Take into account the horizontal covariance between
        different mapped profiles. This lateral covariance takes into account the fact that a mapped
        profile on an Argo position is built from a set of historical profiles that are not very
        different from the set used to build a mapped profile for the next or previous Argo position.
        This lateral covariance between two mapped profiles is constructed using a Guassian function
        and the large spatial scales - Cecile Cabanes

        N.B. Change to code on the xx/xx/2017: Use the small spatial scales instead of the large spatial
        scales to build the lateral covariance. This has been found to be the best comprimise between the
        informative errors and large enough NDF for AIC criterion, at least for the scales defined for the
        North Atlantic basin - Cecile Cabanes

        Parameters
        ----------
        ptmp: matrix of potential temperatures
        coord_float_config: the x, y, z position of the float

        Returns
        -------
        covariance matrix
    """

    # Set up theta boundaries for water masses

    ptboundaries = np.array([30, 24, 18, 12, 8, 4, 2.5, 1, -2])
    ptscale_down = np.array([6, 6, 6, 4, 4, 1.5, 1.5, 1, 1])
    ptscale_up = np.array([6, 6, 6, 6, 4, 4, 1.5, 1.5, 1])

    # Set up the building tile = vertical covariance matrix

    # The upper triangle of the matrix (top right values) are the covariance of each ptlevel with
    # every ptlevel below it, looking down the water column from the diagonal
    # The lower triangle of the matrix (bottom left values) are the covariance of each ptlevel with
    # every ptlevel above it, looking up the water column from the diagonal

    ptmp_rows = ptmp.shape[0]
    ptmp_columns = ptmp.shape[1]

    # set up the covariance matrix
    cov = np.zeros((ptmp_rows * ptmp_columns, ptmp_rows))

    # set up interpolation grids
    upper_interp = interpolate.interp1d(ptboundaries, ptscale_down)
    lower_interp = interpolate.interp1d(ptboundaries, ptscale_up)

    # go through each profile
    for profile in range(0, ptmp_columns):

        profile_1 = profile * ptmp_rows

        # go through each level
        for i in range(0, ptmp_rows):
            for j in range(0, ptmp_rows):

                # belongs in the upper triangle, look down water column for vertical scale
                if i < j:
                    l_theta = upper_interp(ptmp[i, profile])
                    cov[i + profile_1, j] = np.exp(-1 * (ptmp[j, profile] - ptmp[i, profile]) ** 2
                                                   / l_theta ** 2)

                # belongs in the lower triangle, look up water column for vertical scale
                elif i > j:
                    l_theta = lower_interp(ptmp[i, profile])
                    cov[i + profile_1, j] = np.exp(-1 * (ptmp[j, profile] - ptmp[i, profile]) ** 2
                                                   / l_theta ** 2)

                # it is in the leading diagonal, so make it equal to 1
                else:
                    cov[i + profile_1, j] = 1

                # if we don't have a value, make it equal to 1
                if np.isnan(cov[i + profile_1, j]):
                    cov[i + profile_1, j] = 1

    # set up matrix to hold horizontal covariance
    h_cov = np.ones((ptmp_columns, ptmp_columns)) * np.nan

    for profile in range(0, ptmp_columns):
        h_cov[profile, :] = covarxy_pv(coord_float[profile], coord_float,
                                       config['MAPSCALE_LONGITUDE_SMALL'],
                                       config['MAPSCALE_LATITUDE_SMALL'],
                                       config['MAPSCALE_PHI_SMALL'],
                                       config['MAP_USE_PV'])

    h_cov = h_cov[:, 0:ptmp_columns]

    # build final covariance matrix, using horizontal and vertical covariance

    n_cov = np.tile(cov, [1, ptmp_columns])

    # Have to find the covariance for each profile against all other profiles
    for profile in range(0, ptmp_columns):

        lower = profile * ptmp_rows
        upper = (profile + 1) * ptmp_rows

        # go through each profile
        for profile_1 in range(0, ptmp_columns):
            lower_1 = profile_1 * ptmp_rows
            upper_1 = (profile_1 + 1) * ptmp_rows
            n_cov[lower:upper, lower_1:upper_1] = h_cov[profile, profile_1] * \
                                                  n_cov[lower:upper, lower_1:upper_1]

    return n_cov


#pylint: disable=too-many-arguments
def covarxy_pv(input_coords, coords, long, lat, phi, use_pv):
    """ Returns a matrix for the horizontal covariance

        Finds the correlation between spatial and temporal data, and uses this
        to construct the covariance

        Parameters
        ----------
        input_coords: the input coordinates of the the float profile
        coords: coordinates for all the float profiles
        long: longitude scale
        lat: latitude scale
        phi: potential gradient
        use_pv: whether or not to use potential vorticity

        Returns
        -------
        horizontal covariance matrix
    """

    # Derive the planetary vorticity at each point

    # Get the depth for each data point
    z_input_coords = input_coords[2]
    z_coords = coords[:, 2]

    # define a vectorized function to calculation potential vorticity
    potential_vorticity = np.vectorize(lambda latitude, depth:
                                       (2 * 7.292 * 10 ** -5 *
                                        np.sin(latitude * np.pi / 180)) / depth)

    # calculate potential vorticity
    pv_input_coords = potential_vorticity(input_coords[0], z_input_coords)
    pv_coords = potential_vorticity(coords[:, 0], z_coords)

    # calculate correlation
    cor_term = ((input_coords[0] - coords[:, 0]) / lat) ** 2 + \
               ((input_coords[1] - coords[:, 1]) / long) ** 2

    # include potential vorticity in correlation, if the user has asked for it

    if use_pv and pv_input_coords.any() and pv_coords.any() != 0:
        cor_term = cor_term + ((pv_input_coords - pv_coords) /
                               np.sqrt(pv_input_coords ** 2 + pv_coords ** 2) /
                               phi) ** 2

    cov_term = np.exp(-cor_term.transpose())

    return cov_term


#pylint: disable=too-many-locals
def brk_pt_fit(x_obvs, y_obvs, w_i, breaks=None):
    """ Get least-squares estimates for a piecewise linear fit with breakpoints at prescribed points

        Routine to get least squares fit for piecewise linear fit with break points at prescribed points

        y = A(1) + A(2)*x+eps                     for x(1) <= x <= b(1)
        y = A(1) + A(2)*b(1) + A(3)*(x-b(1))+eps  for b(1) <= x <= b(2)
        ...                                            ...
        y = A(1) + A(2)*(b(2)-b(1)) + ... + A(m+2)*(x-b(m)+eps for b(m) <= x <= x(n)

        where:
        x = vector of observed indepedent variables [n]
        y = vector of observed dependent variables  [n]
        W_i = inverse of weights for fit [n] if W_i is empty, then equal weighting for each point
        b = vector of break points as values of x   [m]
        A = fitting coefficients                    [m+1]
        yg = fitted estimate of y

        Parameters
        ----------
        x_obvs: indendent variables
        y_obvs: dependent variables
        w_i: inverse of weights for fit
        breaks: vector of break points

        Returns
        -------
        fit_param, residual: Matrix relating observations to the fit parameters of the linear fit
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
        fit_param = np.dot(np.dot(linalg.solve(ls_est, trends.T), w_i), y_obvs)

    else:
        fit_param = linalg.solve(ls_est, np.dot(trends.T, y_obvs))


    # calculate fit estimate
    residual = y_obvs - np.dot(trends, fit_param)

    return fit_param, residual
