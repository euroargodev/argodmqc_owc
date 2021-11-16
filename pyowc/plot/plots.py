""" Functions to create specific plots
"""
import os
import copy
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from scipy.interpolate import interpolate
import shapefile

from ..core.finders import find_10thetas


def draw_shapes_as_patches(ax, shapes, **kwargs):
    """Draw a series of shapes as PathPatches on a given Matplotlib axis."""
    for shape in shapes:
        points = np.array(shape.points)

        codes = Path.LINETO * np.ones(len(points), dtype=Path.code_type)
        codes[shape.parts] = Path.MOVETO

        path = Path(points, codes)
        patch = PathPatch(path, **kwargs)
        ax.add_patch(patch)


# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
# pylint: disable=too-many-arguments
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

    ax = plt.gca()
    ax.set_aspect(1)

    # load in the coastline data
    coastline = os.path.sep.join([config['CONFIG_DIRECTORY'], "coastline", "ne_10m_coastline.shp"])
    with shapefile.Reader(coastline) as shp:
        shapes = [shape for shape in shp.shapes() if shape.shapeType == shapefile.POLYLINE]
        draw_shapes_as_patches(ax, shapes, linewidth=0.5, edgecolor="black", facecolor="None")

    # if wanted, load in bathymetric data and plot it
    if bath:
        root_path = os.path.join(config['CONFIG_DIRECTORY'], "bathymetry")
        bathymetry_config = {
            "ne_10m_bathymetry_L_0": {"facecolor": "#F0F0F0"},
            "ne_10m_bathymetry_K_200.shp": {"facecolor": "#D2D2D2"},
            "ne_10m_bathymetry_J_1000.shp": {"facecolor": "#B4B4B4"},
            "ne_10m_bathymetry_I_2000.shp": {"facecolor": "#969696"},
            "ne_10m_bathymetry_H_3000.shp": {"facecolor": "#737373"},
            "ne_10m_bathymetry_G_4000.shp": {"facecolor": "#646464"},
            "ne_10m_bathymetry_F_5000.shp": {"facecolor": "#505050"},
            "ne_10m_bathymetry_E_6000.shp": {"facecolor": "#464646"},
            "ne_10m_bathymetry_D_7000.shp": {"facecolor": "#323232"},
            "ne_10m_bathymetry_C_8000.shp": {"facecolor": "#1E1E1E"},
            "ne_10m_bathymetry_B_9000.shp": {"facecolor": "#0A0A0A"},
        }

        for filename, plot_config in bathymetry_config.items():
            with shapefile.Reader(os.path.join(root_path, filename)) as shp:
                shapes = [shape for shape in shp.shapes() if shape.shapeType == shapefile.POLYGON]
                draw_shapes_as_patches(ax, shapes, linewidth=0.0, **plot_config)

    # if we want reef data, load it in and plot it
    if reef:
        reef = os.path.sep.join([config['CONFIG_DIRECTORY'], "reefs", "ne_10m_reefs.shp"])
        with shapefile.Reader(reef) as shp:
            shapes = [shape for shape in shp.shapes() if shape.shapeType == shapefile.POLYGON]
            draw_shapes_as_patches(ax, shapes, linewidth=0.0, facecolor="green", label="Reef")

    # set up the latitude and longitude data
    ax.plot(
        "Longitude",
        "Latitude",
        data=floats,
        color="red",
        marker="o",
        markersize=0.5,
        label="Float Profiles",
        linestyle="-",
    )
    ax.plot(
        "Longitude",
        "Latitude",
        data=climatology,
        color='mediumblue',
        marker="o",
        markersize=0.5,
        label="Climatology",
        linestyle="None",
    )

    plt.title(("Locations of float " + float_name + " with historical data"))
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.axis([-180, 180, -90, 90])
    plt.xlim(np.min(climatology['Longitude']) - 20, np.max(climatology['Longitude']) + 20)
    plt.ylim(np.min(climatology['Latitude']) - 15, np.max(climatology['Latitude']) + 15)

    color = plt.get_cmap('jet')

    # annotate float data (every 3, plus first and last float))
    row_count = len(floats.index)

    for row in floats.itertuples():
        idx = row.Index
        if idx == 0 or idx % 3 == 0 or idx == row_count - 1:
            plt.annotate(
                row.number,
                (row.Longitude, row.Latitude),
                color=color((idx + 1) / row_count),
                size=5,
            )

    plt.legend(loc=4, prop={'size': 6})

    save_format = config['FLOAT_PLOTS_FORMAT']
    plot_loc = os.path.sep.join([config['FLOAT_PLOTS_DIRECTORY'], float_name])
    plt.savefig(plot_loc + "_trajectory." + save_format, format=save_format)


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
def theta_sal_plot(sal, theta, map_sal, map_theta, map_errors,
                   index, profiles, config, float_name, title='uncalibrated'):
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
        profiles: profile numbers array
        config: user configuration
        float_name: name of the float

        Returns
        -------
        Nothing
    """

    # set up plot
    plt.figure(figsize=(10, 6))
    #plt.subplots()
    color_n = sal.__len__()
    colors = pl.cm.jet(np.linspace(0, 1, color_n))

    # can only fit 30 profiles on legend
    n_legend = np.arange(0, profiles[profiles.__len__() - 1], np.ceil(color_n / 30))

    for i in range(sal.__len__()):
        # plot salinities
        if i in n_legend:
            plt.plot(sal[i], theta[i], color=colors[i], label=profiles[i])
        else:
            plt.plot(sal[i], theta[i], color=colors[i])

        good = np.argwhere(~np.isnan(index[:, i]))
        good_index = np.array(index[good, i], dtype=int)

        for j in good_index:
            plt.errorbar(map_sal[j, i],
                         map_theta[j, i],
                         xerr=map_errors[j, i],
                         marker='o', color=colors[i], fillstyle='none', zorder=10)

    # neaten up plot
    plt.legend(loc='center right', bbox_to_anchor=(1.12, 0.5))
    plt.title(title + " float data with mapped salinity and objective errors")
    plt.xlabel("Salinity (PSS-78)")
    plt.ylabel(r"$\theta$ $^\circ$C")

    save_format = config['FLOAT_PLOTS_FORMAT']
    plot_loc = os.path.sep.join([config['FLOAT_PLOTS_DIRECTORY'], float_name])
    plt.savefig(plot_loc + "_" + title + "_theta_sal." + save_format, format=save_format)
    plt.show()


# pylint: disable=too-many-arguments
def t_s_profile_plot(sal, ptmp, pres, sal_var, theta_levels, tlevels, plevels, float_name, config):
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
        config: user configuration

        Returns
        -------
        Nothing
    """

    #plt.figure(1)
    plt.figure(1, figsize=(7, 9))
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

    save_format = config['FLOAT_PLOTS_FORMAT']
    plot_loc = os.path.sep.join([config['FLOAT_PLOTS_DIRECTORY'], float_name])
    plt.savefig(plot_loc + "_salinity_profile." + save_format, format=save_format, bbox_inches='tight')

    plt.show()


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
# pylint: disable=too-many-statements
# pylint: disable=too-many-branches
def sal_var_plot(levels, sal, pres, ptmp, map_sal, map_sal_errors,
                 map_ptmp, cal_sal, cal_sal_errors, boundaries, profile_no,
                 float_name, config):
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
        config: user configuration

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

    thetas = find_10thetas(copy.deepcopy(sal),
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

    if use_theta_lt.__len__() > 0 and use_theta_gt.__len__() == 0:
        good = np.argwhere(ptmp > use_theta_lt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_theta_lt.__len__() == 0 and use_theta_gt.__len__() > 0:
        good = np.argwhere(ptmp < use_theta_gt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_theta_lt.__len__() > 0 and use_theta_gt.__len__() > 0:
        theta_range_lt = (ptmp < use_theta_gt)
        theta_range_gt = (ptmp > use_theta_lt)

        if use_theta_lt < use_theta_gt:
            # exclude middle band
            good = np.argwhere(np.logical_and(theta_range_lt, theta_range_gt))

        else:
            good = np.argwhere(np.logical_or(theta_range_gt, theta_range_lt))

        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_pres_lt.__len__() > 0 and use_pres_gt.__len__() == 0:
        good = np.argwhere(pres > use_pres_lt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_pres_lt.__len__() == 0 and use_pres_gt.__len__() > 0:
        good = np.argwhere(pres < use_pres_gt)
        pres[good] = np.nan
        sal[good] = np.nan
        ptmp[good] = np.nan
        map_sal[good] = np.nan
        map_sal_errors[good] = np.nan
        cal_sal[good] = np.nan
        cal_sal_errors[good] = np.nan

    if use_pres_lt.__len__() > 0 and use_pres_gt.__len__() > 0:
        pres_range_lt = (pres < use_pres_gt)
        pres_range_gt = (pres > use_pres_lt)

        if use_pres_lt < use_pres_gt:
            # exclude middle band
            good = np.argwhere(np.logical_and(pres_range_lt, pres_range_gt))

        else:
            good = np.argwhere(np.logical_or(pres_range_gt, pres_range_lt))

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

                if ptmp[theta_index, j] == thetas[0][i]:
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
        #plt.figure(1)
        plt.figure(figsize=(12, 6))
        good_index = np.isfinite(s_cal[i])
        plt.errorbar(profile_no, s_map[i], yerr=s_map_err[i], color='r', capsize=2)
        plt.fill_between(profile_no[good_index], s_cal[i, good_index] + s_cal_err[i, good_index],
                         s_cal[i, good_index] - s_cal_err[i, good_index], color=(0, 1, 0))
        plt.plot(profile_no, s_int[i, :], marker='o', color='b',
                 label='uncalibrated float')
        plt.plot(profile_no, s_map[i, :], color='r',
                 linewidth=4, zorder=0)
        plt.plot(profile_no, s_cal[i, :], color=(0, 1, 0),
                 label='calibrated float w/ 1xerr', zorder=0)
        plt.plot(profile_no, s_map[i, :], color='r',
                 label='mapped salinity')

        plt.xlabel("Profile number")
        plt.ylabel("PSS-78")

        pl.title(float_name +
                 r" salinities with error on $\theta$=" +
                 str(np.round(thetas[0][i][0], 5)) + r"$\circ$C")

        plt.legend()

        save_format = config['FLOAT_PLOTS_FORMAT']
        plot_loc = os.path.sep.join([config['FLOAT_PLOTS_DIRECTORY'], float_name])
        plt.savefig(plot_loc + "_salinity_variance_" + str(i + 1) + "." + save_format,
                    format=save_format, bbox_inches='tight')

        plt.ylim((np.nanmin(s_int) - 0.05, np.nanmax(s_int) + 0.05))
        plt.show()


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=no-member
def cal_sal_curve_plot(sal, cal_sal, cal_sal_err, sta_sal, sta_sal_err, sta_mean,
                       pcond_factor, pcond_factor_err, profile_no, float_name, config):
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
        config: user configuration

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
        plt.figure(figsize=(12, 7))
        plt.subplot(211)

        plt.plot(profile_no[0], pcond_factor[0], color=(0, 1, 0), linestyle='-', linewidth=2)
        good = np.argwhere(np.isfinite(sta_mean))
        plt.plot(profile_no[0, good[:, 1]], sta_mean[0, good[:, 1]],
                 color=(1, 0, 0), linestyle='-', linewidth=1,
                 label='1-1 profile fit')

        plt.errorbar(profile_no[0], pcond_factor[0], yerr=2 * pcond_factor_err[0],
                     color=(0, 0, 1), linestyle='', capsize=2, label='2 x cal error')
        plt.errorbar(profile_no[0], pcond_factor[0], yerr=pcond_factor_err[0],
                     color=(0, 1, 0), linestyle='', capsize=2, label='1 x cal error')

        # add line at x=1 to help operators

        plt.plot([-1000, 1000], [1, 1], color=(0, 0, 0), linewidth=1.8)
        plt.xlim((np.nanmin(profile_no), np.nanmax(profile_no)))

        plt.legend()
        plt.ylabel('r')
        plt.title(float_name +
                  " potential conductivity (mmho/cm) multiplicative correction r with errors")

        # vertically averaged salinity plot
        plt.subplot(212)
        plt.figure(1)

        plt.plot(profile_no[0], avg_sal_offset, color=(0, 1, 0), linestyle='-', linewidth=2)
        good = np.argwhere(np.isfinite(avg_sta_offset))
        plt.plot(profile_no[0, good[:, 0]], avg_sta_offset[good[:, 0]],
                 color=(1, 0, 0), linestyle='-', linewidth=1, label='1-1 profile fit')
        plt.errorbar(profile_no[0], avg_sal_offset, yerr=2 * avg_sal_offset_err,
                     color=(0, 0, 1), linestyle='', capsize=2, label='2 x cal error')
        plt.errorbar(profile_no[0], avg_sal_offset, yerr=avg_sal_offset_err,
                     color=(0, 1, 0), linestyle='', capsize=2, label='1 x cal error')

        # add line at x=0 to help operators

        plt.plot([-1000, 1000], [0, 0], color=(0, 0, 0), linewidth=1.8)
        plt.xlim((np.nanmin(profile_no), np.nanmax(profile_no)))

        plt.legend()
        plt.ylabel(r'$\Delta$ S')
        plt.xlabel("Profile number")
        plt.title(float_name +
                  r" vertically averaged salinity (PSS-78) additive " +
                  r"correction $\Delta$ S with errors")

        save_format = config['FLOAT_PLOTS_FORMAT']
        plot_loc = os.path.sep.join([config['FLOAT_PLOTS_DIRECTORY'], float_name])
        plt.savefig(plot_loc + "_salinity_curve." + save_format, format=save_format, bbox_inches='tight')

        plt.show()


def sal_anom_plot(sal, ptmp, profile_no, config, float_name, title='uncalibrated'):
    """ Create the salinity anomoly plot

        Parameters
        ----------
        sal: salinity
        ptmp: potential temperature,
        pres: pressure
        profile_no: profile numbers
        config: user configuration
        float_name: name of the float
        title: Addition to the title

        Returns
        -------
        Nothing
    """

    # Set up values and allocate memory for matrices

    theta_base = np.arange(0.1, 30.1, 0.1)
    #ptmp_0 = gsw.conversions.pt0_from_t(sal, ptmp, pres)
    sal_int = np.nan * np.ones((len(theta_base), sal.shape[1]))
    sal_anom = np.nan * np.ones((len(theta_base), sal.shape[1]))
    prof_range = np.nan*np.ones((2, 1), dtype=int).flatten()

    # find anomaly

    good_ptmp = np.where(~np.isnan(ptmp.T[:, :]))
    prof_range[0] = int(np.nanmin(good_ptmp[0]))
    prof_range[1] = int(np.nanmax(good_ptmp[0]))

    for k in range(int(prof_range[0]), int(prof_range[1] + 1)):
        temp = ptmp[:, k]
        sal1 = sal[:, k]
        sal_temp = np.vstack((temp, sal1)).T
        sal_temp_sorted = sal_temp[sal_temp[:, 0].argsort()]

        # make sure the values are unique

        unique = np.unique(sal_temp_sorted[:, 0], return_index=True)[1]
        good_temp = sal_temp_sorted[unique, 0]
        good_sal = sal_temp_sorted[unique, 1]

        good_data = np.argwhere(np.logical_and(~np.isnan(good_temp), ~np.isnan(good_sal)))

        if good_data.__len__() >= 3:
            sal_int_interp = interpolate.interp1d(good_temp[good_data].flatten(), good_sal[good_data].flatten(),
                                                  bounds_error=False)
            sal_int[:, k] = sal_int_interp(theta_base[:])

    # find median and calculate anomaly
    sal_1 = sal_int.T
    sal_size = sal_1.shape[1]
    sal_med = np.nan * np.ones((1, sal_size)).flatten()
    sal_std = np.nan * np.ones((1, sal_size)).flatten()

    for i in range(sal_size):
        good = np.argwhere(~np.isnan(sal_1[:, i]))
        if good.__len__() > 0:
            sal_med[i] = np.nanmedian(sal_1[good, i])
            sal_std[i] = np.nanstd(sal_1[good, i])

    for k in range(int(prof_range[0]), int(prof_range[1] + 1)):
        sal_anom[:, k] = sal_int[:, k] - sal_med

    # create plot

    levels = [-0.1, -0.06, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04, 0.06, 0.1]

    for bounds in config['THETA_BOUNDS']:
        plt.figure(figsize=(10, 6))
        c_bounds = plt.contourf(profile_no[0], theta_base, sal_anom, levels=levels, cmap='seismic')
        plt.contourf(profile_no[0], theta_base, sal_anom, levels=[0.1, 1000], colors='red')
        plt.contourf(profile_no[0], theta_base, sal_anom, levels=[-1000, -0.1], colors='blue')
        plt.contour(profile_no[0], theta_base, sal_anom, levels=[0], colors='black')
        plt.contour(profile_no[0], theta_base, sal_anom, levels=levels,
                    colors='black', linestyles='solid', linewidths=0.25)
        plt.colorbar(c_bounds, label="Difference", ticks=levels)
        plt.xlabel("profile number")
        plt.ylabel("theta")
        plt.ylim((bounds[0], bounds[1]))
        plt.title(title + " salinity anomaly on theta bounds " +
                  str(bounds[0]) + " - " + str(bounds[1]) + " for " + float_name)

        save_format = config['FLOAT_PLOTS_FORMAT']
        plot_loc = os.path.sep.join([config['FLOAT_PLOTS_DIRECTORY'], float_name])
        plt.savefig(plot_loc + "_" + title + "_salinity_anomaly_" + str(bounds) + "." + save_format,
                    format=save_format, bbox_inches='tight')
        plt.show()
