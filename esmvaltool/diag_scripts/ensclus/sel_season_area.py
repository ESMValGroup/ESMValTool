"""Selecting a season (DJF,DJFM,NDJFM,JJA)."""

import numpy as np
import pandas as pd
from netCDF4 import datetime


def sel_season(var, dates, season):
    """Selecting a season (DJF,DJFM,NDJFM,JJA)."""
    # -----------------------------------------------------------------------
    # print('Selecting only {0} data'.format(season))
    dates_pdh = pd.to_datetime(dates)
    if season == 'DJF':       # ONLY DEC-JAN-FEB
        mon = [12, 1, 2]
        mask = (dates_pdh.month == 12) | (dates_pdh.month == 1) |\
               (dates_pdh.month == 2)
    elif season == 'DJFM':    # ONLY DEC-JAN-FEB-MAR
        mon = [12, 1, 2, 3]
        mask = (dates_pdh.month == 12) | (dates_pdh.month == 1) |\
               (dates_pdh.month == 2) | (dates_pdh.month == 3)
    elif season == 'NDJFM':   # ONLY NOV-DEC-JAN-FEB-MAR
        mon = [11, 12, 1, 2, 3]
        mask = (dates_pdh.month == 11) | (dates_pdh.month == 12) |\
               (dates_pdh.month == 1) | (dates_pdh.month == 2) |\
               (dates_pdh.month == 3)
    elif season == 'JJA':   # ONLY JUN-JUL-AUG
        mon = [6, 7, 8]
        mask = (dates_pdh.month == 6) | (dates_pdh.month == 7) |\
               (dates_pdh.month == 8)
    else:
        print('season is not one of the following: DJF, DJFM, NDJFM, JJA')
    var_season = var[mask, :, :]
    dates_season = dates[mask]

    if (12 in mon) or (1 in mon):
        # REMOVING THE FIRST MONTHS (for the first year)
        # because there is no previuos december
        start = int(np.where(dates_season ==
                             datetime(dates_pdh.year[0], mon[0],
                                      dates_pdh.day[0], dates_pdh.hour[0],
                                      dates_pdh.minute[0]))[0])
        # REMOVING THE LAST MONTHS (for the last year)
        # because there is no following january
        end = int(np.where(dates_season ==
                           datetime(dates_pdh.year[-1],
                                    mon[0], dates_pdh.day[0],
                                    dates_pdh.hour[0],
                                    dates_pdh.minute[0]))[0])

        var_season = var_season[start:end, :, :]
        dates_season = dates_season[start:end]

    return var_season, dates_season


# ____________Selecting only [lat_s-lat_n, lon_w-lon_e] box region
def sel_area(lat, lon, var, area):
    """Selecting the area of interest.

    USAGE: var_area, lat_area, lon_area =sel_area(lat,lon,var,area)
    area can be 'EAT', 'PNA', 'NH'
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
        lon_w = lon.min()
        lon_e = lon.max()
        var_roll = var
        lon_new = lon

    elif area == 'Eu':
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
