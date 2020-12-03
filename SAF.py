# If we are near the Subantarctic Front we need to figure out if
                # the profile is north or south of it. Then we should remove data not on
                # the same side the profile is on
import numpy as np
import gsw
from scipy.io import loadmat
import scipy.interpolate as interpolate
import os
def frontal_constraint_saf(config,grid_sal, grid_ptmp, grid_pres, grid_lat, grid_long, grid_dates, grid_z, float_lat, float_z,float_tmp,float_sal):

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
        float_z:Argo profile pressure
        float_tmp:Argo profile absolute temperature
        float_sal:Argo profile salinity
        po_config_data:

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
    s_stdS  = S_stdS[2:]
    t_stdS  = T_stdS[2:]
    s_stdN  = S_stdN[2:]
    t_stdN  = T_stdN[2:]
    t_meanS = T_meanS[2:]
    t_meanN = T_meanN[2:]
    deph     = Deph[2:]

    # SAF is calculated for data below -30 degree
    if float_lat < -30:

        # Calculations for Argo float
        isok = np.argwhere(~np.isnan(float_z) and ~np.isnan(float_tmp))

        # finding differences in depth levels
        float_depth_diff =  np.argwhere(np.diff(float_z[isok]) == 0)

        if (float_depth_diff.__len__() == 0 & isok.__len__() > 2 & np.min(float_z[isok]) < 300 & np.max(float_z[isok]) > 300):

            #SAF:
            T300 = interpolate.interp1d(float_z[isok], float_tmp[isok], 300)

            if T300.__len__() > 0:
                if T300 > 5:
                    critSAF = 1
                elif T300 < 3:
                    critSAF = -1
                else:
                    critSAF = 0
            else:
                critSAF = 0

            # Second step: the envelope test
            if critSAF == 0:
                sal_int = interpolate.interp1d(float_z[isok], float_sal[isok], deph.T)
                temp_int = interpolate.interp1d(float_z[isok], float_tmp[isok], deph.T)

                southinf_int = np.ones ((deph.T,1)) * np.nan
                southsup_int = np.ones ((deph.T,1)) * np.nan
                northinf_int = np.ones ((deph.T,1)) * np.nan
                northsup_int = np.ones ((deph.T,1)) * np.nan

                northinf_int = interpolate.interp1d(t_meanN - t_stdN, s_meanN - s_stdN, temp_int)
                northsup_int = interpolate.interp1d(t_meanN + t_stdN, s_meanN + s_stdN, temp_int)
                southinf_int = interpolate.interp1d(s_meanS - s_stdS, t_meanS - t_stdS, sal_int)
                southsup_int = interpolate.interp1d(s_meanS + s_stdS, t_meanS + t_stdS, sal_int)

                isok2 = np.argwhere(~np.isnan(sal_int) and ~np.isnan(southinf_int) and
                                    ~np.isnan(southsup_int) and ~np.isnan(northinf_int) and
                                    ~np.isnan(northsup_int) and deph.T > 150 and deph.T < 1700)

                if isok2.__len__() > 0:
                    ptSouth = np.argwhere((temp_int[isok2] > southinf_int[isok2] and temp_int[isok2] < southsup_int[isok2]))

                    ptNorth = np.argwhere((sal_int[isok2] > northinf_int[isok2] and sal_int[isok2] < northsup_int[isok2]))

                    isSouth = 0
                    isNorth = 0

                    if len(ptSouth) == len(isok2):
                        isSouth = 1

                    if len(ptNorth) == len(isok2):
                        isNorth = 1

                    if isSouth and isNorth:#np.logical_and(isSouth, isNorth):
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
        grid_critSAF=np.zeros((grid_long,1))

        if critSAF != 0:
            for i in grid_long:
                isok=np.argwhere(~np.isnan(grid_pres[:,i]) and ~np.isnan(grid_sal[:,i]) and ~np.isnan(grid_ptmp[:,i]) and [np.diff(grid_pres[:,i]), 0] != 0)

                grid_press_diff = np.argwhere(np.diff(grid_pres[isok,i]) == 0)
                if (grid_press_diff.__len__() == 0 and isok.__len__() > 2 and np.min(grid_pres[isok,i]) < 300 and np.max(grid_pres[isok,i]) > 300):
                    grid_SA=gsw.conversions.SA_from_SP(grid_sal,grid_pres,grid_long,grid_lat)
                    grid_CT=gsw.conversions.CT_from_pt(grid_SA, grid_ptmp)
                    grid_tmp=gsw.conversions.t_from_CT(grid_SA, grid_CT,grid_pres)

               	    T300=interpolate.interp1d(grid_pres[isok,i],grid_tmp[isok,i],300)
		            if T300 > 5:
			            grid_critSAF[i]=1
                    elif T300 < 3:
			            grid_critSAF[i]=-1
                    else:
                    grid_critSAF[i]=0

        #Test
        # check here use logic_and? !!!
        isSelect=np.argwhere(grid_critSAF==critSAF)

        best_hist_sal = grid_sal[:,isSelect]
        best_hist_ptmp = grid_ptmp[:,isSelect]
        best_hist_pres = grid_pres[:,isSelect]
        best_hist_lat = grid_lat[isSelect]
        best_hist_lon = grid_long[isSelect]
        best_hist_dates	= grid_dates[isSelect]
        best_hist_z		= grid_z[isSelect]

    else: #if Lat<-30
        best_hist_sal = grid_sal
        best_hist_ptmp = grid_ptmp
        best_hist_pres = grid_pres
        best_hist_lat = grid_lat
        best_hist_lon = grid_long
        best_hist_dates = grid_dates
        best_hist_z = grid_z

    return best_hist_sal, best_hist_ptmp, best_hist_pres, best_hist_lat, best_hist_lon, best_hist_dates, best_hist_z
