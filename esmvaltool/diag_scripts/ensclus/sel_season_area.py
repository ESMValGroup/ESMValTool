"""Selecting a season (DJF,DJFM,NDJFM,JJA)."""

import numpy as np


def sel_season(var, dates, season):
    """Selecting a season (DJF,DJFM,NDJFM,JJA).

    USAGE: var_season, dates_season = sel_season(var, dates, season)
    """
    # -----------------------------------------------------------------------
    # print('Selecting only {0} data'.format(season))
    dmonth = np.array([date.month for date in dates])
    if season == 'DJF':       # ONLY DEC-JAN-FEB
        imon = [1, 2]
        emon = [12]
        mask = (dmonth == 12) | (dmonth == 1) | (dmonth == 2)
    elif season == 'DJFM':    # ONLY DEC-JAN-FEB-MAR
        imon = [1, 2, 3]
        emon = [12]
        mask = (dmonth == 12) | (dmonth == 1) | (dmonth == 2) | (dmonth == 3)
    elif season == 'NDJFM':   # ONLY NOV-DEC-JAN-FEB-MAR
        imon = [1, 2, 3]
        emon = [11, 12]
        mask = (dmonth == 11) | (dmonth == 12) | (dmonth == 1) |\
               (dmonth == 2) | (dmonth == 3)
    elif season == 'JJA':   # ONLY JUN-JUL-AUG
        imon = []
        emon = []
        mask = (dmonth == 6) | (dmonth == 7) | (dmonth == 8)
    else:
        print('season is not one of the following: DJF, DJFM, NDJFM, JJA')
    var_season = var[mask, :, :]
    dates_season = dates[mask]

    dmonth = np.array([date.month for date in dates_season])
    dyear = np.array([date.year for date in dates_season])

    imask = list(mon not in imon for mon in dmonth) | (dyear != dyear[0])
    emask = list(mon not in emon for mon in dmonth) | (dyear != dyear[-1])

    var_season = var_season[imask & emask, :, :]
    dates_season = dates_season[imask & emask]

    return var_season, dates_season


# ____________Selecting only [lat_s-lat_n, lon_w-lon_e] box region
def sel_area(lat, lon, var, area):
    """Selecting the area of interest.

    USAGE: var_area, lat_area, lon_area =sel_area(lat,lon,var,area)
    area can be 'EAT', 'PNA', 'NH', 'EU'
    """
    if area == 'EAT':
        # printarea = 'Euro-Atlantic'
        lat_n = 87.5
        lat_s = 30.0
        lon_w = -80.0    # 280
        lon_e = 40.0     # 40
        # lat and lon are extracted from the netcdf file, assumed to be 1D
        # If 0<lon<360, convert to -180<lon<180
        if lon.min() >= 0:
            lon_new = lon - 180
            var_roll = np.roll(var, int(len(lon) / 2), axis=2)
        else:
            var_roll = var
            lon_new = lon

    elif area == 'PNA':
        # printarea = 'Pacific North American'
        lat_n = 87.5
        lat_s = 30.0
        lon_w = 140.0
        lon_e = 280.0
        # lat and lon are extracted from the netcdf file, assumed to be 1D
        # If -180<lon<180, convert to 0<lon<360
        if lon.min() < 0:
            lon_new = lon + 180
            var_roll = np.roll(var, int(len(lon) / 2), axis=2)
        else:
            var_roll = var
            lon_new = lon

    elif area == 'NH':
        # printarea = 'Northern Hemisphere'
        lat_n = 90.0
        lat_s = 0.0
        # lat and lon are extracted from the netcdf file, assumed to be 1D
        # If 0<lon<360, convert to -180<lon<180
        if lon.min() >= 0:
            lon_new = lon - 180
            var_roll = np.roll(var, int(len(lon) / 2), axis=2)
        else:
            var_roll = var
            lon_new = lon
        lon_w = lon_new.min()
        lon_e = lon_new.max()

    elif area == 'EU':
        # printarea = 'Europe'
        lat_n = 72.0
        lat_s = 27.0
        lon_w = -22.0
        lon_e = 45.0
        # lat and lon are extracted from the netcdf file, assumed to be 1D
        # If 0<lon<360, convert to -180<lon<180
        if lon.min() >= 0:
            lon_new = lon - 180
            var_roll = np.roll(var, int(len(lon) / 2), axis=2)
        else:
            var_roll = var
            lon_new = lon

    # -------------------------------------------------------------------
    # print('__________________________________________________________')
    # print('Selecting the area of interest: {0}'.format(printarea))
    # --------------------------------------------------------------------
    # -------------------------Selecting only an area

    latidx = (lat >= lat_s) & (lat <= lat_n)
    lonidx = (lon_new >= lon_w) & (lon_new <= lon_e)

    var_area = var_roll[:, latidx][..., lonidx]
    # print('Grid dimension of the selected area --->
    #         {0}'.format(var_area[0].shape))

    return var_area, lat[latidx], lon_new[lonidx]
