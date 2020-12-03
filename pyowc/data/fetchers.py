""" Functions to access, load data

        Parameters
        ----------

        Returns
        -------

"""

import warnings
import os
import struct
import numpy as np
from scipy.io import loadmat
import scipy.interpolate as interpolate

from ..utilities import wrap_longitude, change_dates
from ..configuration import print_cfg

#pylint: disable=too-many-locals
def get_topo_grid(min_long, max_long, min_lat, max_lat, config):
    """  Find depth grid over given area using tbase.int file

        The old matlab version of this uses an old .int file from NOAA which contains
        5 arcminute data of global terrain. Whilst other more complete data sets now
        exists, the data files for them are very large, so this file can be used for now,
        and perhaps we can update to better data when we either move the code online, or
        find a way to susbet it/compress it reasonably.

        The .int file is very weird, and stores 16bit integers in binary. The below
        code opens the file and converts the old binary back to the numbers we expect.
        It then finds global terrain over a specified area before closing the file

        Parameters
        ----------
        min_long: minimum longitudinal value for grid
        max_long: maximum longitudinal value for grid
        min_lat: minimum latidunal value for grid
        max_lat: maximum latidunal value for grid

        Returns
        -------
        Matrices containing a uniform grid of latitudes and longitudes, along with the depth at these points
    """
    # manipulate input values to match file for decoding
    blat = int(np.max((np.floor(min_lat * 12), -90 * 12 + 1)))
    tlat = int(np.ceil(max_lat * 12))
    llong = int(np.floor(min_long * 12))
    rlong = int(np.ceil(max_long * 12))

    # use these values to form the grid
    lgs = np.arange(llong, rlong + 1, 1) / 12
    lts = np.flip(np.arange(blat, tlat + 1, 1) / 12, axis=0)

    if rlong > 360 * 12 - 1:
        rlong = rlong - 360 * 12
        llong = llong - 360 * 12

    if llong < 0:
        rlong = rlong + 360 * 12
        llong = llong + 360 * 12

    decoder = [llong, rlong, 90 * 12 - blat, 90 * 12 - tlat]

    # get the amount of elevation values we need
    nlat = int(round(decoder[2] - decoder[3])) + 1
    nlong = int(round(decoder[1] - decoder[0])) + 1

    # initialise matrix to hold z values
    topo = np.zeros((nlat, nlong))

    # Open the binary file
    elev_file = open(os.path.sep.join([config['CONFIG_DIRECTORY'], "tbase.int"]), "rb")

    # decode the file, and get values
    for i in range(nlat):
        elev_file.seek((i + decoder[3]) * 360 * 12 * 2 + decoder[0] * 2)

        for j in range(nlong):
            topo[i, j] = struct.unpack('h', elev_file.read(2))[0]

    # make the grid
    longs, lats = np.meshgrid(lgs, lts)

    # close the file
    elev_file.close()

    return topo, longs, lats


def get_data(wmo_box, data_type, config, pa_float_name):
    """ Gets all the data from regions specified by the WMO boxes

        Uses the WMO box numbers to fetch the relevant data for the region.

        This is a refactor between get_region_data and get_region_hist_locations to avoid
        duplicate code

        Relies on configuration parameters:
            - config['HISTORICAL_DIRECTORY']
            - config['HISTORICAL_CTD_PREFIX']
            - config['HISTORICAL_BOTTLE_PREFIX']
            - config['HISTORICAL_ARGO_PREFIX']

        Parameters
        ----------
        wmo_box: 2D array containing the name of the WMO boxes that cover the area of interest, and flags for whether we want to use argo, bottle, and/or CTD data
        data_type: Which type of data we are checking: 1 (CTD), 2 (BOTTLE) or 3 (ARGO)
        config: Dictionary containing configuration settings. Used to find locations of folders and file containing data
        pa_float_name: String containing the name of the float being profiled

        Returns
        -------
        Array containing all (ctd/bottle/argo) data in a wmo box
    """
    # warnings.warn(print_cfg(config))
    # warnings.warn("pa_float_name: %s (%s)" % (str(pa_float_name), type(pa_float_name)))

    data = []
    box_name = str(int(wmo_box[0]))
    # check if we should use this data. If so, get the data
    if wmo_box[data_type] == 1 and data_type == 1:
        data = loadmat(os.path.sep.join([config['HISTORICAL_DIRECTORY'],
                                         config['HISTORICAL_CTD_PREFIX'] + box_name + '.mat']))

    if wmo_box[data_type] == 1 and data_type == 2:
        data = loadmat(os.path.sep.join([config['HISTORICAL_DIRECTORY'],
                                         config['HISTORICAL_BOTTLE_PREFIX'] + box_name + '.mat']))

        # if data dimensions don't match, transpose to avoid indexing issues
        if (data['lat'].size == data['pres'].size and
                data['lat'].shape[1] == data['pres'].shape[1]):
            data['pres'] = data['pres'].T
            data['ptmp'] = data['ptmp'].T
            data['sal'] = data['sal'].T
            data['temp'] = data['temp'].T

    if wmo_box[data_type] == 1 and data_type == 3:
        data = loadmat(os.path.sep.join([config['HISTORICAL_DIRECTORY'],
                                         config['HISTORICAL_ARGO_PREFIX'] + box_name + '.mat']))

        # remove the argo float being analysed from the data
        not_use = []
        for i in range(0, data['lat'][0].__len__()):
            if str(data['source'][0][i]).find(pa_float_name) != -1:
                not_use.append(i)

        data['lat'] = [np.delete(data['lat'], not_use)]
        data['long'] = [np.delete(data['long'], not_use)]
        data['dates'] = [np.delete(data['dates'], not_use)]
        data['sal'] = [np.delete(data['sal'], not_use, axis=1)]
        data['ptmp'] = [np.delete(data['ptmp'], not_use, axis=1)]
        data['pres'] = [np.delete(data['pres'], not_use, axis=1)]

    return data


#pylint: disable=bare-except
def get_region_hist_locations(pa_wmo_numbers, pa_float_name, config):
    """ Uses the WMO boxes and to return all of the historical data in the given area, excluding the float that is currently being analysed.

        Function returns the needed spatial and temporal data from each of the historical
        data points inside the WMO boxes it is given, excluding the current float being processed.
        The WMO boxes are passed in the form of [WMO name, use ctd data, use bottle data, use argo data],
        and the configuration is also passed in, which is used to locate the different data locally. The
        configuration can be changed in load_configuration.py

        N.B. Change to code on xx/xx/2015 - Only load longitude, latitude, and dates

        Parameters
        ----------
        pa_wmo_numbers: 2D array containing the name of the WMO boxes that cover the area of interest, and flags for
            whether we want to use argo, bottle, and/or CTD data
        pa_float_name: string of the name of the float currently being processed
        config: Dictionary containing configuration settings. Used to find locations of folders and file containing data

        Returns
        -------
        The latitude, longitude and age of each data point we want
    """

    # set up matrices to hold data
    grid_lat = []
    grid_long = []
    grid_dates = []

    # go through each of the WMO boxes
    for wmo_box in pa_wmo_numbers:

        # go through each of the columns denoting whether we should use CTD (1), bottle (2), and/or argo (3)
        for data_type in range(1, 4):

            # get the data
            try:
                data = get_data(wmo_box, data_type, config, pa_float_name)

                # if we have data, combine it with the other data then reset it
                if data:
                    grid_lat = np.concatenate([grid_lat, data['lat'][0]])
                    grid_long = np.concatenate([grid_long, data['long'][0]])
                    grid_dates = np.concatenate([grid_dates, data['dates'][0]])
                    data = []

            except:
                pass

    if grid_lat.__len__() == 0:

        raise ValueError("get_region_hist_locations found no data for your specification. "
                         "Are your wmo_boxes files set up correctly?\n%s" % str(config)) from None



    grid_long = wrap_longitude(grid_long)
    # decimalise dates
    grid_dates = change_dates(grid_dates)

    return grid_lat, grid_long, grid_dates


#pylint: disable=too-many-arguments
#pylint: disable=too-many-locals
#pylint: disable=too-many-branches
#pylint: disable=too-many-statements
#pylint: disable=too-many-nested-blocks
#pylint: disable=bare-except
def get_region_data(pa_wmo_numbers, pa_float_name, config, index, pa_float_pres):
    """ Get the historical pressure, salinity, and temperature of selected casts

        Function returns the needed spatial, temporal, and ocean characteristic (salinity, temp, pressure)
        data from each of the historical data points inside the WMO boxes it is given, excluding the
        current float being processed. It then filters out the data that was not chosen as a good fit.

        The WMO boxes are passed in the form of [WMO name, use ctd data,
        use bottle data, use argo data], and the configuration is also passed in, which is used to locate
        the different data locally. The configuration can be changed in load_configuration.py

        Isabelle Gaboury, 26 Sep. 2017: Added check on the dimensions of bottle data.


        Parameters
        ----------
        pa_wmo_numbers: 2D array containing the name of the WMO boxes that cover the area of interest, and flags for
            whether we want to use argo, bottle, and/or CTD data
        pa_float_name: string of the name of the float currently being processed
        config: Dictionary containing configuration settings. Used to find locations of folders and file containing data
        index: array of indices of selected historical casts
        pa_float_pres: array of pressures for the float being processed

        Returns
        -------
        The salinity, potential temperature, pressure, latitude, longitude, and age of each historical cast selected to use
    """

    # maximum depth to retrieve data from (deepest float measurement + MAP_P_DELTA
    max_pres = np.nanmax(pa_float_pres) + config["MAP_P_DELTA"]

    # set up empty arrays to hold the data to return
    grid_sal = []
    grid_ptmp = []
    grid_pres = []
    grid_lat = []
    grid_long = []
    grid_dates = []
    data = []

    # set up current maximum depth and number of columns
    max_depth = 0
    how_many_cols = 0

    # set up variable to save beginning index for each set of data
    starting_index = 0

    # go through each of the WMO boxes
    for wmo_box in pa_wmo_numbers:

        # go through each of the columns denoting whether we should use CTD, bottle, and/or argo
        for data_type in range(1, 4):

            # get the data
            try:
                data = get_data(wmo_box, data_type, config, pa_float_name)

                if data:
                    # Sometimes the data comes in wrapped in as a 3d array, so convert to 2d
                    if data['pres'].__len__() == 1:
                        data['pres'] = data['pres'][0]
                        data['sal'] = data['sal'][0]
                        data['ptmp'] = data['ptmp'][0]
                        data['lat'] = data['lat'][0].reshape(-1, 1)
                        data['long'] = data['long'][0].reshape(-1, 1)
                        data['dates'] = data['dates'][0].reshape(-1, 1)

                    #  check the index of each station to see if it should be loaded
                    data_length = data['lat'][0].__len__()
                    data_indices = np.arange(0, data_length) + starting_index

                    # remember location of last entry
                    starting_index = starting_index + data_length

                    # load each station
                    for i in range(0, data_length):

                        good_indices = np.argwhere(index == data_indices[i])

                        if good_indices.__len__() > 0:
                            # only use non-NaN values
                            not_nan = np.argwhere(np.isnan(data['pres'][:, i]) == 0)

                            # get the non-NaN values
                            pres = data['pres'][not_nan, i]
                            sal = data['sal'][not_nan, i]
                            ptmp = data['ptmp'][not_nan, i]

                            # remove values where pressure exceeds the maximum we want
                            # todo: Why is this not done during data fetching ????
                            too_deep = np.argwhere(pres > max_pres)
                            pres = np.delete(pres, too_deep[:, 0])
                            sal = np.delete(sal, too_deep[:, 0])
                            ptmp = np.delete(ptmp, too_deep[:, 0])
                            new_depth = pres.__len__()
                            how_many_rows = np.max([new_depth, max_depth])

                            # if the new data we are adding is longer than our columns, we need to
                            # fill in NaNs in the other columns
                            if new_depth > max_depth != 0:
                                grid_pres = np.append(grid_pres, np.ones(
                                    (how_many_cols, new_depth - max_depth)) * np.nan, axis=1
                                                      ).reshape((how_many_cols, how_many_rows))
                                grid_ptmp = np.append(grid_ptmp, np.ones(
                                    (how_many_cols, new_depth - max_depth)) * np.nan, axis=1
                                                      ).reshape((how_many_cols, how_many_rows))
                                grid_sal = np.append(grid_sal, np.ones(
                                    (how_many_cols, new_depth - max_depth)) * np.nan, axis=1
                                                     ).reshape((how_many_cols, how_many_rows))

                            # if the new data we are adding is shorter than our columns,
                            # then we need to fill in the rest with NaNs so it's the same length
                            elif new_depth < max_depth:
                                pres = np.append(pres, np.ones((max_depth - new_depth, 1)) * np.nan)
                                ptmp = np.append(ptmp, np.ones((max_depth - new_depth, 1)) * np.nan)
                                sal = np.append(sal, np.ones((max_depth - new_depth, 1)) * np.nan)

                            # if we don't have any data saved yet, create the grid matrix with the
                            # first data set
                            if grid_pres.__len__() == 0:

                                grid_pres = pres.reshape((1, pres.__len__()))
                                grid_ptmp = ptmp.reshape((1, pres.__len__()))
                                grid_sal = sal.reshape((1, pres.__len__()))

                            # if we already have data saved, add the new data to the saved data
                            else:

                                grid_pres = np.append(grid_pres, pres).reshape(
                                    how_many_cols + 1, how_many_rows)
                                grid_ptmp = np.append(grid_ptmp, ptmp).reshape(
                                    how_many_cols + 1, how_many_rows)
                                grid_sal = np.append(grid_sal, sal).reshape(
                                    how_many_cols + 1, how_many_rows)

                            # save the latitude, longitude, and date of the new data
                            grid_lat = np.append(grid_lat, data['lat'][0, i])
                            grid_long = np.append(grid_long, data['long'][0, i])
                            grid_dates = np.append(grid_dates, data['dates'][0, i])

                            # readjust our values so we know what column to add the new data to,
                            # and what shape we should expect the data to be
                            max_depth = grid_pres.shape[1]
                            how_many_cols = grid_pres.shape[0]

            except:
                pass

    # convert longitude to 0 to 360 degrees
    try:
        grid_long = wrap_longitude(grid_long)

        # make sure salinity, pressure, and potential temperature data have all the same NaNs
        sal_nans = np.argwhere(np.isnan(grid_sal))
        for nan in sal_nans:
            grid_pres[nan[0], nan[1]] = np.nan
            grid_ptmp[nan[0], nan[1]] = np.nan

        pres_nans = np.argwhere(np.isnan(grid_pres))
        for nan in pres_nans:
            grid_sal[nan[0], nan[1]] = np.nan
            grid_ptmp[nan[0], nan[1]] = np.nan

        ptmp_nans = np.argwhere(np.isnan(grid_ptmp))
        for nan in ptmp_nans:
            grid_sal[nan[0], nan[1]] = np.nan
            grid_pres[nan[0], nan[1]] = np.nan

        grid_dates = change_dates(grid_dates)

        # transpose data
        grid_sal = grid_sal.T
        grid_pres = grid_pres.T
        grid_ptmp = grid_ptmp.T

    except:
        raise RuntimeError("NO DATA FOUND") from None
    # we have encountered a problem where some data coming in is all NaN
    # these columns need to be removed from the data set
    nans = 0
    for column in range(grid_sal.shape[1]):

        if np.all(np.isnan(grid_sal[:, column - nans])) or \
                np.all(np.isnan(grid_pres[:, column - nans])):

            grid_sal = np.delete(grid_sal, column - nans, 1)
            grid_ptmp = np.delete(grid_ptmp, column - nans, 1)
            grid_pres = np.delete(grid_pres, column - nans, 1)
            grid_lat = np.delete(grid_lat, column - nans)
            grid_long = np.delete(grid_long, column - nans)
            grid_dates = np.delete(grid_dates, column - nans)
            nans += 1

    if nans > 0:
        print("Warning: found ", nans,
              " all NaNs in your dataset. These water columns have been removed")

    return grid_sal, grid_ptmp, grid_pres, grid_lat, grid_long, grid_dates


def frontal_constraint_saf(config, grid_sal, grid_ptmp, grid_pres, grid_lat, grid_long, grid_dates, grid_z, float_lat, float_pres, float_tmp,float_sal):

    """ Function to

        Parameters
        ----------
        grid_sal: practical salinity of historical data
        grid_ptmp: potential temperature of  historical data
        grid_pres: pressure of  historical data
        grid_lat: latitude of  historical data
        grid_long:longitude of  historical data
        grid_dates: dates of  historical data
        grid_z: depth of historical data
        float_lat: Argo profile latitude
        float_pres:Argo profile pressure
        float_tmp:Argo profile absolute temperature
        float_sal:Argo profile salinity
        config: Dictionary containing configuration settings. Used to find locations of folders and file containing data

        Returns
        -------
        Array of historical data after application of the SAF frontal criterion

    """

    typical_profile_around_saf = loadmat(os.path.sep.join([config['CONFIG_DIRECTORY'], config['CONFIG_SAF']]))

    S_meanS = typical_profile_around_saf['S_meanS']
    S_meanN = typical_profile_around_saf['S_meanN']
    S_stdS = typical_profile_around_saf['S_stdS']
    T_stdS = typical_profile_around_saf['T_stdS']
    S_stdN = typical_profile_around_saf['S_stdN']
    T_stdN = typical_profile_around_saf['T_stdN']
    T_meanS = typical_profile_around_saf['T_meanS']
    T_meanN = typical_profile_around_saf['T_meanN']
    Deph = typical_profile_around_saf['Deph']

    #keep only data below 100m depth data :
    s_meanS = S_meanS[2:]
    s_meanN = S_meanN[2:]
    s_stdS = S_stdS[2:]
    t_stdS = T_stdS[2:]
    s_stdN = S_stdN[2:]
    t_stdN = T_stdN[2:]
    t_meanS = T_meanS[2:]
    t_meanN = T_meanN[2:]
    deph = Deph[2:]

    # SAF is calculated for data below -30 degree
    if float_lat < -30:

        # Calculations for Argo float
        isok = np.argwhere(~np.isnan(float_pres) & ~np.isnan(float_tmp))

        if np.argwhere(np.diff(float_pres[isok]) == 0).__len__() == 0 and isok.__len__() > 2 and np.min(float_pres[isok]) < 300 and np.max(float_pres[isok]) > 300:

            #SAF:
            t300 = np.interp(300, float_pres[isok].flatten(), float_tmp[isok].flatten())

            if t300 > 0:
                if t300 > 5:
                    critSAF = 1
                elif t300 < 3:
                    critSAF = -1
                else:
                    critSAF = 0
            else:
                critSAF = 0

            # Second step: the envelope test
            if critSAF == 0:
                sal_int = np.interp(deph.T.flatten(), float_pres[isok].flatten(), float_sal[isok].flatten())
                temp_int = np.interp(deph.T.flatten(), float_pres[isok].flatten(), float_tmp[isok].flatten())

                southinf_int = np.ones((deph.T, 1)) * np.nan
                southsup_int = np.ones((deph.T, 1)) * np.nan
                northinf_int = np.ones((deph.T, 1)) * np.nan
                northsup_int = np.ones((deph.T, 1)) * np.nan

                northinf_int = np.interp(temp_int.flatten(), t_meanN - t_stdN, s_meanN - s_stdN)
                northsup_int = np.interp(temp_int.flatten(), t_meanN + t_stdN, s_meanN + s_stdN)
                southinf_int = np.interp(sal_int.flatten(), s_meanS - s_stdS, t_meanS - t_stdS)
                southsup_int = np.interp(sal_int.flatten(), s_meanS + s_stdS, t_meanS + t_stdS)

                isok2 = np.argwhere(~np.isnan(sal_int) & ~np.isnan(southinf_int) & ~np.isnan(southsup_int)
                                    & ~np.isnan(northinf_int) & ~np.isnan(northsup_int)
                                    & deph.T > 150 & deph.T < 1700)

                if isok2.__len__() > 0:
                    ptSouth = np.argwhere((temp_int[isok2] > southinf_int[isok2] & temp_int[isok2] < southsup_int[isok2]))

                    ptNorth = np.argwhere((sal_int[isok2] > northinf_int[isok2] & sal_int[isok2] < northsup_int[isok2]))

                    isSouth = 0
                    isNorth = 0

                    if len(ptSouth) == len(isok2):
                        isSouth = 1

                    if len(ptNorth) == len(isok2):
                        isNorth = 1

                    if isSouth & isNorth:#np.logical_and(isSouth, isNorth):
                        critSAF = 0

                    elif (isSouth):
                        critSAF = -1

                    elif (isNorth):
                        critSAF = 1
            else:
                critSAF = 0
        else:
            critSAF = 0


        # Historical data
        grid_critSAF = [0] * len(grid_long)

        if critSAF != 0:
            for i in grid_long:
                isok = np.argwhere(~np.isnan(grid_pres[:, i]) & ~np.isnan(grid_sal[:, i])
                                 & ~np.isnan(grid_ptmp[:, i]) & [np.diff(grid_pres[:, i]), 0] != 0)

                if np.argwhere(np.diff(grid_pres[isok, i]) == 0).__len__() == 0 and isok.__len__() > 2 and np.min(grid_pres[isok, i]) < 300 and np.max(grid_pres[isok, i]) > 300 :

                    grid_SA = gsw.conversions.SA_from_SP(grid_sal, grid_pres, grid_long, grid_lat)
                    grid_CT = gsw.conversions.CT_from_pt(grid_SA, grid_ptmp)
                    grid_tmp = gsw.conversions.t_from_CT(grid_SA, grid_CT, grid_pres)

                    t300 = np.interp(300, grid_pres[isok, i].flatten(), grid_tmp[isok, i].flatten())

                    if t300 > 5:
                        grid_critSAF[i] = 1
                    elif t300 < 3:
                        grid_critSAF[i] = -1
                    else:
                        grid_critSAF[i] = 0

        #Test
        # check here use logic_and? !!!
        isSelect = np.argwhere(grid_critSAF == critSAF)

        best_hist_sal = grid_sal[:, isSelect]
        best_hist_ptmp = grid_ptmp[:, isSelect]
        best_hist_pres = grid_pres[:, isSelect]
        best_hist_lat = grid_lat[isSelect]
        best_hist_lon = grid_long[isSelect]
        best_hist_dates = grid_dates[isSelect]
        best_hist_z = grid_z[isSelect]
    # if Lat<-30
    else:
        best_hist_sal = grid_sal
        best_hist_ptmp = grid_ptmp
        best_hist_pres = grid_pres
        best_hist_lat = grid_lat
        best_hist_lon = grid_long
        best_hist_dates = grid_dates
        best_hist_z = grid_z

    return best_hist_sal, best_hist_ptmp, best_hist_pres, best_hist_lat, best_hist_lon, best_hist_dates, best_hist_z