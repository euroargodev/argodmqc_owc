""" Functions to create specific plots
"""
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import geopandas as gdp # pylint: disable=import-error
from pyowc import core
from scipy.interpolate import interpolate


# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
def trajectory_plot(bath, reef, floats, climatology, float_name, config):
    """ Plot locations of all the data used in the analysis

        function for plotting locations of all the data used in the analysis, including:
        
        historical climatology
        profile locations and order
        
        Can also plot reef data and bathymetry by passing in a 1 (True) into the function

        Parameters
        ----------
        bath: should we plot bathymetry (1 == yes)
        reef: should we plot reefs (1 == yes)
        floats: float location data frame
        climatology: climatology location dataframe
        float_name: name of float

        Returns
        -------
        Nothing
    """

    # load in the coastline data
    coastline = config['CONFIG_DIRECTORY'] + "coastline/ne_10m_coastline.shp"
    map_coast = gdp.read_file(coastline)
    traj_map = map_coast.plot(color='black', label='coastline')

    # if wanted, load in bathymetric data and plot it
    if bath == 1:
        bathymetry0 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_L_0.shp"
        bathymetry200 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_K_200.shp"
        bathymetry1000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_J_1000.shp"
        bathymetry2000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_I_2000.shp"
        bathymetry3000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_H_3000.shp"
        bathymetry4000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_G_4000.shp"
        bathymetry5000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_F_5000.shp"
        bathymetry6000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_E_6000.shp"
        bathymetry7000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_D_7000.shp"
        bathymetry8000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_C_8000.shp"
        bathymetry9000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_B_9000.shp"
        bathymetry10000 = config['CONFIG_DIRECTORY'] + "bathymetry/ne_10m_bathymetry_A_10000.shp"
        map_bath0 = gdp.read_file(bathymetry0)
        map_bath200 = gdp.read_file(bathymetry200)
        map_bath1000 = gdp.read_file(bathymetry1000)
        map_bath2000 = gdp.read_file(bathymetry2000)
        map_bath3000 = gdp.read_file(bathymetry3000)
        map_bath4000 = gdp.read_file(bathymetry4000)
        map_bath5000 = gdp.read_file(bathymetry5000)
        map_bath6000 = gdp.read_file(bathymetry6000)
        map_bath7000 = gdp.read_file(bathymetry7000)
        map_bath8000 = gdp.read_file(bathymetry8000)
        map_bath9000 = gdp.read_file(bathymetry9000)
        map_bath10000 = gdp.read_file(bathymetry10000)
        traj_map = map_bath0.plot(ax=traj_map, color='#BEBEBE', label='>200m', linewidth=2)
        traj_map = map_bath200.plot(ax=traj_map, color='#B8B8B8', linewidth=2)
        traj_map = map_bath1000.plot(ax=traj_map, color='#B0B0B0', label='1000m', linewidth=2)
        traj_map = map_bath2000.plot(ax=traj_map, color='#A9A9A9')
        traj_map = map_bath3000.plot(ax=traj_map, color='#A8A8A8')
        traj_map = map_bath4000.plot(ax=traj_map, color='#A0A0A0')
        traj_map = map_bath5000.plot(ax=traj_map, color='#989898')
        traj_map = map_bath6000.plot(ax=traj_map, color='#909090', label='6000m')
        traj_map = map_bath7000.plot(ax=traj_map, color='#888888')
        traj_map = map_bath8000.plot(ax=traj_map, color='#808080')
        traj_map = map_bath9000.plot(ax=traj_map, color='#787878')
        traj_map = map_bath10000.plot(ax=traj_map, color='#707070')

    # if we want reef data, load it in and plot it
    if reef == 1:
        reef = config['CONFIG_DIRECTORY'] + "reefs/ne_10m_reefs.shp"
        map_reef = gdp.read_file(reef)
        traj_map = map_reef.plot(ax=traj_map, color='green', label='reef')

    # set up the latitude and longitude data
    geo_floats = gdp.GeoDataFrame(floats,
                                  geometry=gdp.points_from_xy(floats.Longitude,
                                                              floats.Latitude))
    geo_climatology = gdp.GeoDataFrame(climatology,
                                       geometry=gdp.points_from_xy(climatology.Longitude,
                                                                   climatology.Latitude))

    traj_map = geo_floats.plot(ax=traj_map, color='red', marker="+", label='profile')
    geo_climatology.plot(ax=traj_map, color='#00008B', marker="s",
                         markersize=12, label='climatology')
    plt.plot(floats['Longitude'], floats['Latitude'], color='red', linestyle='-')
    plt.title(("Locations of float " + float_name + " with historical data"))
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.axis([-180, 180, -90, 90])
    plt.xlim(np.min(climatology['Longitude']) - 10, np.max(climatology['Longitude']) + 10)
    plt.ylim(np.min(climatology['Latitude']) - 10, np.max(climatology['Latitude']) + 10)

    for i, txt in enumerate(floats['number']):
        plt.annotate(txt, (floats['Longitude'][i], floats['Latitude'][i]))

    plt.legend(loc='center left', bbox_to_anchor=(0.85, 0.85))
    plt.draw()


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
def theta_sal_plot(sal, theta, map_sal, map_theta, map_errors, index, title='uncalibrated'):
    """ Create the salinity theta curve

        Parameters
        ----------
        index: index of theta levels with least variance
        title: Addition to the title
        sal: float salinity
        theta: float potential temperature
        map_sal: mapped salinity
        map_theta: mapped potential temperature
        map_errors: mapped salinity errors

        Returns
        -------
        Nothing
    """

    # set up plot
    plt.subplots()
    color_n = sal.__len__()
    colors = pl.cm.jet(np.linspace(0, 1, color_n))

    # can only fit 30 profiles on legend
    n_legend = np.arange(0, 30, np.ceil(color_n / 30))

    for i in range(sal.__len__()):
        # plot salinities
        if i in n_legend:
            plt.plot(sal[i], theta[i], color=colors[i], label=i)
        else:
            plt.plot(sal[i], theta[i], color=colors[i])

        good = np.argwhere(~np.isnan(index[:, i]))
        good_index = np.array(index[good, i], dtype=int)

        for j in good_index:
            plt.errorbar(map_sal[j, i],
                         map_theta[j, i],
                         xerr=map_errors[j, i],
                         marker='o', color=colors[i], fillstyle='none')

    # neaten up plot
    plt.legend(loc='center right', bbox_to_anchor=(1.12, 0.5))
    plt.title(title + " float data with mapped salinity and objective errors")
    plt.xlabel("Salinity (PSS-78)")
    plt.ylabel(r"$\theta$ $^\circ$C")
    plt.show()


# pylint: disable=too-many-arguments
def t_s_profile_plot(sal, ptmp, pres, sal_var, theta_levels, tlevels, plevels, float_name):
    """ Plots profile plots

        Parameters
        ----------
        sal: salinity
        ptmp: potential temperature
        pres: pressure
        sal_var: salininty variance on theta levels
        theta_levels: theta levels
        tlevels: temp at theta levels
        plevels: pressure at theta levels
        float_name: name of the float

        Returns
        -------
        Nothing
    """

    plt.figure(1)

    # plot t-s profile
    plt.subplot(222)

    plt.plot(sal, ptmp, color='b', linestyle='-', linewidth=0.5)
    plt.xlabel("PSS-78")
    plt.title(" OW chosen levels - " + float_name)

    for i in tlevels:
        plt.axhline(y=i, color=(0, 1, 0), linestyle='-', linewidth=0.5)

    # plot s_variance on t
    plt.subplot(221)

    plt.plot(sal_var, theta_levels, color='b', linewidth=0.5)
    plt.xlabel("Salinity variance")
    plt.ylabel(r"Potential temp ($^{\circ}$C")
    plt.title("Salinity Variance on Theta")

    for i in tlevels:
        plt.axhline(y=i, color=(0, 1, 0), linestyle='-', linewidth=0.5)

    # plot p-t profile
    plt.subplot(223)
    plt.plot(ptmp, -1 * pres, color='b', linewidth=0.5)
    plt.xlabel(r"$^{\circ}$C")
    plt.ylabel("Pressure (dbar)")
    plt.title("OW chosen levels - " + float_name)

    for i in plevels:
        plt.axhline(y=-i, color=(0, 1, 0), linestyle='-', linewidth=0.5)

    # plot p-s profile
    plt.subplot(224)
    plt.plot(sal, -1 * pres, color='b', linewidth=0.5)
    plt.xlabel("PSS-78")
    plt.title("OW chosen levels - " + float_name)

    for i in plevels:
        plt.axhline(y=-i, color=(0, 1, 0), linestyle='-', linewidth=0.5)

    plt.tight_layout(pad=1)
    plt.show()


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
# pylint: disable=too-many-statements
# pylint: disable=too-many-branches
def sal_var_plot(levels, sal, pres, ptmp, map_sal, map_sal_errors,
                 map_ptmp, cal_sal, cal_sal_errors, boundaries, profile_no, float_name):
    """ Create the salinity variance plot for each level

        Parameters
        ----------
        float_name: name of the float
        profile_no: profile numbers
        map_ptmp: mapped potential temperature
        levels: number of levels to plot
        sal: float salinity
        pres: float pressure
        ptmp: float potential temperature
        map_sal: mapped salinity
        map_sal_errors: mapped salinity errors
        cal_sal: calibrated salinity
        cal_sal_errors: calibrated salinity errors
        boundaries: pressure and temperature boundaries

        Returns
        -------
        Nothing
    """

    # set up matrices
    profile_no = profile_no.flatten()
    no_profiles = pres.shape[1]
    s_int = np.nan * np.ones((levels, no_profiles))
    s_map = np.nan * np.ones((levels, no_profiles))
    s_map_err = np.nan * np.ones((levels, no_profiles))
    s_cal = np.nan * np.ones((levels, no_profiles))
    s_cal_err = np.nan * np.ones((levels, no_profiles))
    thetalevel_index = np.nan * np.ones((levels, no_profiles))

    # Find levels on which we should plot
    use_theta_lt = boundaries[0]
    use_theta_gt = boundaries[1]
    use_pres_lt = boundaries[2]
    use_pres_gt = boundaries[3]
    use_percent_gt = boundaries[4]

    thetas = core.finders.find_10thetas(copy.deepcopy(sal),
                           copy.deepcopy(ptmp),
                           copy.deepcopy(pres),
                           copy.deepcopy(map_ptmp),
                           use_theta_lt, use_theta_gt,
                           use_pres_lt, use_pres_gt,
                           use_percent_gt)

    bad = np.argwhere(np.isnan(map_ptmp))
    for i in bad:
        pres[i[0], i[1]] = np.nan
        sal[i[0], i[1]] = np.nan
        ptmp[i[0], i[1]] = np.nan
        map_sal[i[0], i[1]] = np.nan
        map_sal_errors[i[0], i[1]] = np.nan
        cal_sal[i[0], i[1]] = np.nan
        cal_sal_errors[i[0], i[1]] = np.nan

    if use_theta_lt != 0 and use_theta_gt == 0:
        good = np.argwhere(ptmp > use_theta_lt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_theta_lt == 0 and use_theta_gt != 0:
        good = np.argwhere(ptmp < use_theta_gt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_theta_lt != 0 and use_theta_gt != 0:
        if use_theta_gt > use_theta_lt:
            good = np.argwhere(use_theta_gt > ptmp > use_theta_lt)

        else:
            good = np.argwhere(ptmp < use_theta_gt or ptmp > use_theta_lt)

        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_pres_lt != 0 and use_pres_gt == 0:
        good = np.argwhere(pres > use_pres_lt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_pres_lt == 0 and use_pres_gt != 0:
        good = np.argwhere(pres < use_pres_gt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_pres_lt != 0 and use_pres_gt != 0:
        if use_pres_gt > use_pres_lt:
            good = np.argwhere(use_pres_gt > pres > use_pres_lt)

        else:
            good = np.argwhere(pres < use_pres_gt or pres > use_pres_lt)

        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    for i in range(no_profiles):
        for j in range(levels):

            if np.nanmax(ptmp[:, i]) > thetas[0][j] > np.nanmin(ptmp[:, i]):
                diff_theta = np.abs(ptmp[:, i] - thetas[0][j])

                if np.argwhere(~np.isnan(diff_theta)).__len__() == 0:
                    thetalevel_index[j, i] = np.nan

                else:
                    thetalevel_index[j, i] = np.nanmin(np.argwhere(diff_theta ==
                                                                   np.nanmin(diff_theta)))

    # Build s matrix for plotting
    for i in range(levels):
        for j in range(no_profiles):
            theta_index = thetalevel_index[i, j]

            if ~np.isnan(theta_index):
                theta_index = int(theta_index)
                inter = np.arange(np.max([theta_index - 1, 0]),
                                  np.min([theta_index + 1, pres.shape[0] - 1]) + 1,
                                  dtype=int)

                ptmp_diff = ptmp[theta_index, j] - ptmp[inter, j]

                if ptmp[theta_index, j] > thetas[0][i]:
                    pos_diff = np.argwhere(ptmp_diff > 0)

                    if pos_diff.__len__() > 0:
                        min_diff = np.argwhere(ptmp_diff == np.nanmin(ptmp_diff[pos_diff]))
                        k_index = inter[min_diff]

                    else:
                        k_index = theta_index

                if ptmp[theta_index, j] < thetas[0][i]:
                    neg_diff = np.argwhere(ptmp_diff < 0)

                    if neg_diff.__len__() > 0:
                        min_diff = np.argwhere(-ptmp_diff == np.nanmin(-ptmp_diff[neg_diff]))
                        k_index = inter[min_diff]

                    else:
                        k_index = theta_index

                else:
                    k_index = theta_index

                if ((k_index != theta_index and ~np.isnan(sal[theta_index, j])) and
                        ~np.isnan(sal[k_index, j]) and ~np.isnan(ptmp[theta_index, j]) and
                        ~np.isnan(ptmp[k_index, j])):

                    interp_ptmp_sal = interpolate.interp1d([ptmp[theta_index, j],
                                                            ptmp[k_index, j]],
                                                           [sal[theta_index, j],
                                                            sal[k_index, j]])

                    s_int[i, j] = interp_ptmp_sal(thetas[0][i][0])

                else:
                    s_int[i, j] = sal[theta_index, j]

                if ((k_index != theta_index and ~np.isnan(map_sal[theta_index, j])) and
                        ~np.isnan(map_sal[k_index, j]) and ~np.isnan(ptmp[theta_index, j]) and
                        ~np.isnan(ptmp[k_index, j])):

                    interp_map_sal = interpolate.interp1d([ptmp[theta_index, j],
                                                           ptmp[k_index, j]],
                                                          [map_sal[theta_index, j],
                                                           map_sal[k_index, j]])
                    s_map[i, j] = interp_map_sal(thetas[0][i][0])

                    interp_map_sal_err = interpolate.interp1d([ptmp[theta_index, j],
                                                               ptmp[k_index, j]],
                                                              [map_sal_errors[theta_index, j],
                                                               map_sal_errors[k_index, j]])

                    s_map_err[i, j] = interp_map_sal_err(thetas[0][i][0])

                else:
                    s_map[i, j] = map_sal[theta_index, j]
                    s_map_err[i, j] = map_sal_errors[theta_index, j]

                if ((k_index != theta_index and ~np.isnan(cal_sal[theta_index, j])) and
                        ~np.isnan(cal_sal[k_index, j]) and ~np.isnan(ptmp[theta_index, j]) and
                        ~np.isnan(ptmp[k_index, j])):

                    interp_cal_sal = interpolate.interp1d([ptmp[theta_index, j],
                                                           ptmp[k_index, j]],
                                                          [cal_sal[theta_index, j],
                                                           cal_sal[k_index, j]])
                    interp_cal_sal_err = interpolate.interp1d([ptmp[theta_index, j],
                                                               ptmp[k_index, j]],
                                                              [cal_sal_errors[theta_index, j],
                                                               cal_sal_errors[k_index, j]])

                    s_cal[i, j] = interp_cal_sal(thetas[0][i][0])
                    s_cal_err[i, j] = interp_cal_sal_err(thetas[0][i][0])

                else:
                    s_cal[i, j] = cal_sal[theta_index, j]
                    s_cal_err[i, j] = cal_sal_errors[theta_index, j]

    # plot data (one plot for each theta level, as selected by user)
    for i in range(levels):
        plt.figure(1)
        plt.plot(profile_no, s_int[i, :], marker='o', color='b',
                 label='uncalibrated float')
        plt.plot(profile_no, s_map[i, :], color='r',
                 linewidth=4, zorder=0)
        plt.plot(profile_no, s_cal[i, :], color=(0, 1, 0),
                 label='calibrated float w/ 1xerr')
        plt.errorbar(profile_no, s_map[i, :], yerr=s_map_err[i, :], color='r', capsize=2)
        plt.fill_between(profile_no, s_cal[i, :] - s_cal_err[i, :],
                         s_cal[i, :] + s_cal_err[i, :], color=(0, 1, 0))
        plt.plot(profile_no, s_map[i, :], color='r',
                 label='mapped salinity')

        plt.xlabel("Profile number")
        plt.ylabel("PSS-78")

        pl.title(float_name +
                 r" salinities with error on $\theta$=" +
                 str(np.round(thetas[0][i][0], 5)) + r"$\circ$C")

        plt.legend()
        plt.show()


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
def cal_sal_curve_plot(sal, cal_sal, cal_sal_err, sta_sal, sta_sal_err, sta_mean,
                       pcond_factor, pcond_factor_err, profile_no, float_name):
    """ Create the calibrated salinity curve plot

        Parameters
        ----------
        sta_mean: mean of the differences
        profile_no: profile numbers
        sta_sal_err: error in mean difference
        cal_sal_err: error in calibrated salinity
        sal: float salinity
        cal_sal: calibrated salinity
        sta_sal: mean difference between salinity and calculated salinity per profile
        pcond_factor: slope
        pcond_factor_err: slope error
        float_name: name of the float

        Returns
        -------
        Nothing
    """

    # only plot this if we have calibrated salinity
    if np.argwhere(~np.isnan(cal_sal)).__len__() > 0:
        # set up matrices
        profiles = sal.shape[1]
        sal_offset = cal_sal - sal
        avg_sal_offset = np.nan * np.ones((1, profiles)).flatten()
        avg_sal_offset_err = np.nan * np.ones((1, profiles)).flatten()
        sta_offset = sta_sal - sal
        avg_sta_offset = np.nan * np.ones((1, profiles)).flatten()
        avg_sta_offset_err = np.nan * np.ones((1, profiles)).flatten()

        for i in range(profiles):
            avg_sal_offset[i] = np.nanmean(sal_offset[:, i])
            avg_sal_offset_err[i] = np.nanmean(cal_sal_err[:, i])
            avg_sta_offset[i] = np.nanmean(sta_offset[:, i])
            avg_sta_offset_err[i] = np.nanmean(sta_sal_err[:, i])

        # begin plotting data (two plots)

        # potential conductivity plot
        plt.figure(1)
        plt.subplot(211)

        plt.plot(profile_no[0], pcond_factor[0], color=(0, 1, 0), linestyle='-', linewidth=1)
        good = np.argwhere(np.isfinite(sta_mean))
        plt.plot(profile_no[0, good[:, 1]], sta_mean[0, good[:, 1]],
                 color=(1, 0, 0), linestyle='-', linewidth=1,
                 label='1-1 profile fit')

        plt.errorbar(profile_no[0], pcond_factor[0], yerr=2 * pcond_factor_err[0],
                     color=(0, 0, 1), linestyle='', capsize=2, label='2 x cal error')
        plt.errorbar(profile_no[0], pcond_factor[0], yerr=pcond_factor_err[0],
                     color=(0, 1, 0), linestyle='', capsize=2, label='1 x cal error')

        plt.legend()
        plt.ylabel('r')
        plt.title(float_name +
                  " potential conductivity (mmho/cm) multiplicative correction r with errors")

        # vertically averaged salinity plot
        plt.subplot(212)
        plt.figure(1)

        plt.plot(profile_no[0], avg_sal_offset, color=(0, 1, 0), linestyle='-', linewidth=1)
        good = np.argwhere(np.isfinite(avg_sta_offset))
        plt.plot(profile_no[0, good[:, 0]], avg_sta_offset[good[:, 0]],
                 color=(1, 0, 0), linestyle='-', linewidth=1, label='1-1 profile fit')
        plt.errorbar(profile_no[0], avg_sal_offset, yerr=2 * avg_sal_offset_err,
                     color=(0, 0, 1), linestyle='', capsize=2, label='2 x cal error')
        plt.errorbar(profile_no[0], avg_sal_offset, yerr=avg_sal_offset_err,
                     color=(0, 1, 0), linestyle='', capsize=2, label='1 x cal error')

        plt.legend()
        plt.ylabel(r'$\Delta$ S')
        plt.xlabel("Profile number")
        plt.title(float_name +
                  r" vertically averaged salinity (PSS-78) additive " +
                  r"correction $\Delta$ S with errors")

        plt.show()