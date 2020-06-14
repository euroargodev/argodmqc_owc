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

global A, breaks, nbr1, ubrk_g, xf, yf, w_i, xblim


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
