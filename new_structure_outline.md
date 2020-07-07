Overall new structure of the library

# New positionning of the existing functions

pyowc/core
    utils.py: change_dates, cal2dec, potential_vorticity, wrap_longitudes, sorter
    stats.py: brk_pt_fit, build_cov, covarxy_pv, calc_piecewisefit, covar_xyt_pv, noise_variance, signal_variance, spatial_correlation, fit_cond
    finders.py: find_10thetas, find_25boxes, find_besthit, find_ellipse

pyowc/data
    fetchers.py: get_region_data, get_region_hist_locations, get_data, get_topo_grid
    wrangling.py: interp_climatology, map_data_grid 

pyowc/calibration.py: set_calseries, update_salinity_mapping

pyowc/configuration.py: load_configuration

pyowc/errors.py  # List of specific errors

pyowc/plot
    dashboard.py: plot_diagnostics
    plots.py: cal_sal_curve_plot, sal_var_plot, t_s_profile_plot, theta_sal_plot, trajectory_plot
    utils.py: create_dataframe

pyowc/tests  # Contains all the tests !


