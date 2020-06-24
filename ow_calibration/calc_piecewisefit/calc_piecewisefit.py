"""
-----Calculate Piece wise fit-----

Written by: Annie Wong
When: xx/10/2008
Converted to python by: Edward Small
When: 16/06/2020

Calculate the fit of each break and calibrate salinities

https://github.com/ArgoDMQC/matlab_owc
https://gitlab.noc.soton.ac.uk/edsmall/bodc-dmqc-python
"""

import copy
import gsw
import numpy as np
import scipy.io as scipy
import scipy.interpolate as interpolate
from ow_calibration.build_cov.build_cov import build_cov
from ow_calibration.find_10thetas.find_10thetas import find_10thetas
from ow_calibration.fit_cond.fit_cond import fit_cond
from ow_calibration.tbase_decoder.tbase_decoder import get_topo_grid


# pylint: disable=invalid-name
# pylint: disable=too-many-locals
# pylint: disbale=unused-variable
# pylint: disable=too-many-branches
# pylint: disable=too-many-satements
def calc_piecewisefit(float_dir, float_name, system_config):
    """
    calibrate salinities
    :param float_dir: float directory name
    :param float_name: name of float
    :param system_config: configuration parameter set up
    :return: Nothing, save output
    """

    # load in the source data

    float_source_data = scipy.loadmat(system_config['FLOAT_SOURCE_DIRECTORY'] +
                                      float_dir + float_name +
                                      system_config['FLOAT_SOURCE_POSTFIX'])

    lat = float_source_data['LAT']
    long = float_source_data['LONG']
    sal = float_source_data['SAL']
    ptmp = float_source_data['PTMP']
    pres = float_source_data['PRES']
    profile_no = float_source_data['PROFILE_NO']
    x_in = np.tile(profile_no, (10, 1))

    # load in the mapped data
    float_mapped_data = scipy.loadmat(system_config['FLOAT_MAPPED_DIRECTORY'] +
                                      float_dir + system_config['FLOAT_MAPPED_PREFIX'] +
                                      float_name + system_config['FLOAT_MAPPED_POSTFIX'])

    mapped_sal = float_mapped_data['la_mapped_sal']
    mapsalerror = float_mapped_data['la_mapsalerrors']
    mapped_ptmp = float_mapped_data['la_ptmp']
    selected_hist = float_mapped_data['selected_hist']

    # retrieve XYZ of float position used to build covariance
    if selected_hist.__len__() > 0:
        if long.shape[0] > 1:
            long = long.flatten()

        if lat.shape[0] > 1:
            lat = lat.flatten()

        if np.any(long > 180):
            long_1 = copy.deepcopy(long) - 360

        else:
            long_1 = copy.deepcopy(long)

        elev, x_grid, y_grid = get_topo_grid(np.nanmin(long_1) - 1, np.nanmax(long_1) + 1,
                                             np.nanmin(lat) - 1, np.nanmax(lat) + 1)

        grid_interp = interpolate.interp2d(x_grid[0, :], y_grid[:, 0],
                                           elev, kind='linear')

        z_grid = []
        for i in range(long_1[0].__len__()):
            z_grid.append(grid_interp(long_1[0][i], lat[0][i]))

        z_grid = -np.array(z_grid)
        coord_float = np.column_stack((long.T, lat.T, z_grid))

    # load the calibration settings

    float_calseries = scipy.loadmat(system_config['FLOAT_CALIB_DIRECTORY'] + float_dir +
                                    system_config['FLOAT_CALSERIES_PREFIX'] + float_name +
                                    system_config['FLOAT_MAPPED_POSTFIX'])

    calseries = float_calseries['calseries']
    max_breaks = float_calseries['max_breaks']
    breaks = float_calseries['breaks']
    use_theta_gt = float_calseries['use_theta_gt']
    use_theta_lt = float_calseries['use_theta_lt']
    use_pres_gt = float_calseries['use_pres_gt']
    use_pres_lt = float_calseries['use_pres_lt']
    use_percent_gt = float_calseries['use_percent_gt']

    m, n = pres.shape

    cal_sal = np.ones((m, n)) * np.nan
    cal_sal_err = np.ones((m, n)) * np.nan
    cal_cond = np.ones((m, n)) * np.nan
    cal_cond_err = np.ones((m, n)) * np.nan
    pcond_factor = np.ones((1, n)) * np.nan
    pcond_factor_err = np.ones((1, n)) * np.nan
    time_deriv = np.ones((1, n)) * np.nan
    time_deriv_err = np.ones((1, n)) * np.nan
    sta_mean = np.ones((1, n)) * np.nan
    sta_rms = np.ones((1, n)) * np.nan
    sta_sal = np.ones((m, n)) * np.nan
    sta_sal_err = np.ones((m, n)) * np.nan
    fceof = []
    fbreaks = []

    sstatus = 1
    unique_cal = np.unique(calseries)
    # bad profiles are flagged as zero
    bad = np.argwhere(unique_cal == 0)

    if bad.__len__() > 0:
        unique_cal[bad] = np.delete(unique_cal, bad)

    n_seq = unique_cal.__len__()
    if n_seq == 1 and max_breaks.__len__() > 1:
        print("Error in specificying number of possible break points")
        print(str(max_breaks), " specified, should be ",
              str([max_breaks.__len__(), n_seq]))

    # we have multiple cal series, make sure that break information is provideed for all segments
    elif n_seq > 1:
        # only one max break specified, specify this break for all segements
        if max_breaks.__len__() == 1:
            max_breaks = np.ones((n_seq, 1)) * max_breaks

        # error in specification of max breaks
        elif max_breaks.__len__() != n_seq:
            print("Error in specifying the number of possible break points")
            print(str(max_breaks), " specified, should be 1 or ",
                  str([max_breaks.__len__(), n_seq]))
            sstatus = 0

    if breaks.__len__() > 0:
        ns, nb = breaks.shape

        # error in specifying breaks
        if ns != n_seq:
            print("Error in specifying break points")
            print("For multiple cal series, need to specify breaks for each series")
            print("Have ", str(n_seq), " or ", str(ns), " sets of breaks")

            sstatus = 0

        for n in range(n_seq):
            nb = np.argwhere(np.isfinite(breaks[n, :])).__len__()

            if nb > max_breaks[n]:
                print("Error, for cal series ", str(unique_cal[n]), "max number of breaks ",
                      str(max_breaks[n]), " less than ", str(nb), "prescribed breaks")
                sstatus = 0

            elif nb < max_breaks[n]:
                print("Specified ", str(nb), " breaks. Will search up to ",
                      str(max_breaks[n]), " breaks")

            else:
                print(str(nb), "fixed breaks prescribed")

    # set_calseries returned a bad status variable, write out file with NaNs
    if sstatus == 0:
        float_calib_filename = (system_config['FLOAT_CALIB_DIRECTORY'] + float_dir +
                                system_config['FLOAT_CALIB_PREFIX'] + float_name +
                                system_config['FLOAT_CALIB_POSTFIX'])

        scipy.savemat(float_calib_filename,
                      {'cal_SAL': cal_sal,
                       'cal_SAL_err': cal_sal_err,
                       'pcond_factor': pcond_factor,
                       'pcond_factor_err': pcond_factor_err,
                       'cal_COND': cal_cond,
                       'cal_COND_err': cal_cond_err,
                       'time_deriv': time_deriv,
                       'time_deriv_err': time_deriv_err,
                       'sta_mean': sta_mean,
                       'sta_rms': sta_rms,
                       'sta_SAL': sta_sal,
                       'sta_SAL_err': sta_sal_err,
                       'PROFILE_NO': profile_no,
                       'fcoef': fceof,
                       'fbreaks': fbreaks})

        return

    # loop through sequences of calseries

    for i in range(n_seq):
        calindex = np.argwhere(calseries == unique_cal[i])[:, 1]
        k = calindex.__len__()

        # chose 10 float theta levels to use for the piecewise linear fit
        unique_coord_float = coord_float[calindex, :]
        unique_sal = sal[:, calindex]
        unique_ptmp = ptmp[:, calindex]
        unique_pres = pres[:, calindex]
        unique_mapped_ptmp = mapped_ptmp[:, calindex]
        unique_mapped_sal = mapped_sal[:, calindex]
        unique_mapsalerrors = mapsalerror[:, calindex]

        ten_sal = np.ones((10, k)) * np.nan
        ten_ptmp = np.ones((10, k)) * np.nan
        ten_pres = np.ones((10, k)) * np.nan
        ten_mapped_sal = np.ones((10, k)) * np.nan
        ten_mapsalerrors = np.ones((10, k)) * np.nan

        # make deep copies for calibration layter
        unique_sal_1 = copy.deepcopy(unique_sal)
        unique_ptmp_1 = copy.deepcopy(unique_ptmp)
        unique_pres_1 = copy.deepcopy(unique_pres)
        unique_mapped_ptmp_1 = copy.deepcopy(unique_ptmp)
        unique_mapped_sal_1 = copy.deepcopy(unique_mapped_sal)
        unique_mapsalerrors_1 = copy.deepcopy(unique_mapsalerrors)

        theta, p, index, var_s_th, th = find_10thetas(copy.deepcopy(unique_sal),
                                                      copy.deepcopy(unique_ptmp),
                                                      copy.deepcopy(unique_pres),
                                                      copy.deepcopy(unique_ptmp),
                                                      use_theta_lt[0, 0], use_theta_gt[0, 0],
                                                      use_pres_lt[0, 0], use_pres_gt[0, 0],
                                                      use_percent_gt[0, 0])

        index = np.array(index, dtype=int)
        pp = np.argwhere(np.isnan(index) == 0)
        # only proceed if we have valied theta levels
        if pp.__len__() > 0:
            for ipr in range(k):
                jj = np.argwhere(index[:, ipr] >= 0)
                if jj.__len__() > 0:
                    ten_sal[0:jj.__len__(), ipr] = unique_sal[index[jj, ipr], ipr].flatten()
                    ten_ptmp[0:jj.__len__(), ipr] = unique_ptmp[index[jj, ipr], ipr].flatten()
                    ten_pres[0:jj.__len__(), ipr] = unique_pres[index[jj, ipr], ipr].flatten()
                    ten_mapped_sal[0:jj.__len__(), ipr] = unique_mapped_sal[index[jj, ipr],
                                                                            ipr].flatten()
                    ten_mapsalerrors[0:jj.__len__(), ipr] = unique_mapsalerrors[index[jj, ipr],
                                                                                ipr].flatten()
            # calculate potential conductivites and errors for mapped values and float values
            # calculate pcond error by perturbing salinity
            # (avoids problems caused by non-linearity)

            # constant for conductivity at sal=35, temp=15 and pres=0
            sw_c3515 = 42.914

            icond = gsw.conversions.C_from_SP(ten_sal,
                                              ten_ptmp,
                                              0)
            mapped_cond = gsw.conversions.C_from_SP(ten_mapped_sal,
                                                    ten_ptmp,
                                                    0)

            mapped_cond1 = gsw.conversions.C_from_SP(ten_mapped_sal + ten_mapsalerrors / 100,
                                                     ten_ptmp, 0)

            mapconderrors = 100 * np.abs(mapped_cond - mapped_cond1)

            # independent variable for pieve wise fit (profile number)
            x = x_in[:, calindex]
            y = mapped_cond / icond
            err = mapconderrors / icond

            # calculate off-diagonal terms for error estimate

            covariance = build_cov(ten_ptmp, unique_coord_float, system_config)

            # if no break points are set
            if breaks.__len__() == 0:
                (xfit, condslope, condslope_err, time_deriv, time_deriv_err,
                 sta_mean, sta_rms, ndf, fit_coef, fit_breaks) = fit_cond(x, y, err,
                                                                          covariance,
                                                                          'max_no_breaks',
                                                                          max_breaks[i][0])
                pcond_factor[0][calindex] = condslope
                pcond_factor_err[0][calindex] = condslope_err
                time_deriv[calindex] = time_deriv
                time_deriv_err[calindex] = time_deriv_err
                sta_mean[0][calindex], = sta_mean
                sta_rms[0][calindex] = sta_rms

            else:
                breaks_in = breaks[i, :]
                breaks_in = breaks_in[np.argwhere(np.isfinite(breaks_in))]

                if max_breaks[i]:
                    (xfit, condslope, condslope_err, time_deriv, time_deriv_err,
                     sta_mean, sta_rms, ndf, fit_coef, fit_breaks) = fit_cond(x, y, err,
                                                                              covariance,
                                                                              'breaks',
                                                                              breaks_in,
                                                                              'max_no_breaks',
                                                                              max_breaks[i][0])
                    pcond_factor[0][calindex] = condslope
                    pcond_factor_err[0][calindex] = condslope_err
                    time_deriv[calindex] = time_deriv
                    time_deriv_err[calindex] = time_deriv_err
                    sta_mean[0][calindex], = sta_mean
                    sta_rms[0][calindex] = sta_rms

                else:
                    (xfit, condslope, condslope_err, time_deriv, time_deriv_err,
                     sta_mean, sta_rms, ndf, fit_coef, fit_breaks) = fit_cond(x, y, err,
                                                                              covariance,
                                                                              'breaks',
                                                                              breaks_in)
                    pcond_factor[0][calindex] = condslope
                    pcond_factor_err[0][calindex] = condslope_err
                    time_deriv[calindex] = time_deriv
                    time_deriv_err[calindex] = time_deriv_err
                    sta_mean[0][calindex], = sta_mean
                    sta_rms[0][calindex] = sta_rms

            # apply calibrations to float data

            if pcond_factor[0][calindex].__len__() > 0:
                unique_cond = gsw.conversions.C_from_SP(unique_sal_1, unique_ptmp_1, 0)
                cal_cond[:, calindex] = np.dot(np.ones((m, 1)),
                                               pcond_factor[:, calindex]) * unique_cond
                cal_sal[:, calindex] = gsw.conversions.SP_from_C(cal_cond[:, calindex],
                                                                 unique_ptmp_1,
                                                                 0)
                cal_cond_err[:, calindex] = np.dot(np.ones((m, 1)),
                                                   pcond_factor_err[:, calindex]) * unique_cond
                cal_sal1 = gsw.conversions.SP_from_C((cal_cond[:, calindex] +
                                                      cal_cond_err[:, calindex]),
                                                     unique_ptmp, 0)

                cal_sal_err[:, calindex] = np.abs(cal_sal[:, calindex] - cal_sal1[:, calindex])

                # estimate the error in salinity for station by fit

                sta_cond = np.dot(np.ones((m, 1)), sta_mean[:, calindex]) * unique_cond
                sta_sal[:, calindex] = gsw.conversions.SP_from_C(sta_cond, unique_ptmp, 0)
                sta_cond_err = np.dot(np.ones((m, 1)), sta_rms[:, calindex]) * unique_cond
                sta_sal1 = gsw.conversions.SP_from_C(sta_cond + sta_cond_err, unique_ptmp, 0)
                sta_sal_err[:, calindex] = np.abs(sta_sal[:, calindex] - sta_sal1)

                for n in range(fit_coef.__len__()):
                    fceof.append(fit_coef[0])

                if fit_breaks.__len__() > 0:
                    fbreaks.append(fit_breaks)

    # save calibration data

    float_calib_name = (system_config['FLOAT_CALIB_DIRECTORY'] + float_dir +
                        system_config['FLOAT_CALIB_PREFIX'] + float_name +
                        system_config['FLOAT_CALIB_POSTFIX'])

    scipy.savemat(float_calib_name, {'cal_SAL': cal_sal,
                                     'cal_SAL_err': cal_sal_err,
                                     'pcond_factor': pcond_factor,
                                     'pcond_factor_err': pcond_factor_err,
                                     'cal_COND': cal_cond,
                                     'cal_COND_err': cal_cond_err,
                                     'time_deriv': time_deriv,
                                     'time_deriv_err': time_deriv_err,
                                     'sta_mean': sta_mean,
                                     'sta_rms': sta_rms,
                                     'sta_SAL': sta_sal,
                                     'sta_SAL_err': sta_sal_err,
                                     'PROFILE_NO': profile_no,
                                     'fcoef': fceof,
                                     'fbreaks': fbreaks})
