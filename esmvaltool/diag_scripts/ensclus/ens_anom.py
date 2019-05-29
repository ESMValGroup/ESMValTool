"""Computation of ensemble anomalies based on a desired value."""

import os
import numpy as np
from scipy import stats

# User-defined packages
from read_netcdf import read_iris, save_n_2d_fields
from sel_season_area import sel_area, sel_season


def ens_anom(filenames, dir_output, name_outputs, varname, numens, season,
             area, extreme):
    """Ensemble anomalies.

    Computation of the ensemble anomalies based on the desired value
    from the input variable (it can be the percentile, mean, maximum, standard
    deviation or trend)
    OUTPUT: NetCDF files of ensemble mean of climatology, selected value and
    anomaly maps.
    """
    print('The name of the output files will be <variable>_{0}.txt'
          .format(name_outputs))
    print('Number of ensemble members: {0}'.format(numens))

    outfiles = []
    # Reading the netCDF file of 3Dfield, for all the ensemble members
    var_ens = []
    for ens in range(numens):
        ifile = filenames[ens]
        # print('ENSEMBLE MEMBER %s' %ens)
        var, varunits, lat, lon, dates, _ = read_iris(ifile)

        # Convertion from kg m-2 s-1 to mm/day
        if varunits == 'kg m-2 s-1':
            var = var * 86400  # there are 86400 seconds in a day
            varunits = 'mm/day'

        # Selecting a season (DJF,DJFM,NDJFM,JJA)
        var_season, _ = sel_season(var, dates, season)

        # Selecting only [latS-latN, lonW-lonE] box region
        var_area, lat_area, lon_area = sel_area(lat, lon, var_season, area)

        var_ens.append(var_area)

    if varunits == 'kg m-2 s-1':
        print('\nPrecipitation rate units were converted from kg m-2 s-1 '
              'to mm/day')

    print('The variable is {0} ({1})'.format(varname, varunits))
    print('Original var shape: (time x lat x lon)={0}'.format(var.shape))
    print('var shape after selecting season {0} and area {1}: '
          '(time x lat x lon)={2}'.format(season, area, var_area.shape))

    if extreme == 'mean':
        # Compute the time mean over the entire period, for each ens member
        varextreme_ens = [np.nanmean(var_ens[i], axis=0)
                          for i in range(numens)]

    elif len(extreme.split("_")) == 2:
        # Compute the chosen percentile over the period, for each ens member
        quant = int(extreme.partition("th")[0])
        varextreme_ens = [np.nanpercentile(var_ens[i], quant, axis=0)
                          for i in range(numens)]

    elif extreme == 'maximum':
        # Compute the maximum value over the period, for each ensemble member
        varextreme_ens = [np.nanmax(var_ens[i], axis=0) for i in range(numens)]

    elif extreme == 'std':
        # Compute the standard deviation over the period, for each ens member
        varextreme_ens = [np.nanstd(var_ens[i], axis=0) for i in range(numens)]

    elif extreme == 'trend':
        # Compute the linear trend over the period, for each ensemble member
        trendmap = np.empty((var_ens[0].shape[1], var_ens[0].shape[2]))
        trendmap_ens = []
        for i in range(numens):
            for jla in range(var_ens[0].shape[1]):
                for jlo in range(var_ens[0].shape[2]):
                    slope, _, _, _, _ = \
                        stats.linregress(range(var_ens[0].shape[0]),
                                         var_ens[i][:, jla, jlo])
                    trendmap[jla, jlo] = slope
            trendmap_ens.append(trendmap.copy())
        varextreme_ens = trendmap_ens

    varextreme_ens_np = np.array(varextreme_ens)
    print('Anomalies are computed with respect to the {0}'.format(extreme))

    # Compute and save the anomalies with respect to the ensemble
    ens_anomalies = varextreme_ens_np - np.nanmean(varextreme_ens_np, axis=0)
    varsave = 'ens_anomalies'
    ofile = os.path.join(dir_output, 'ens_anomalies_{0}.nc'
                         .format(name_outputs))
    # print(ofile)
    print('ens_anomalies shape: (numens x lat x lon)={0}'
          .format(ens_anomalies.shape))
    save_n_2d_fields(lat_area, lon_area, ens_anomalies, varsave,
                     varunits, ofile)
    outfiles.append(ofile)
    # Compute and save the climatology
    vartimemean_ens = [np.mean(var_ens[i], axis=0) for i in range(numens)]
    ens_climatologies = np.array(vartimemean_ens)
    varsave = 'ens_climatologies'
    ofile = os.path.join(dir_output, 'ens_climatologies_{0}.nc'
                         .format(name_outputs))
    save_n_2d_fields(lat_area, lon_area, ens_climatologies, varsave,
                     varunits, ofile)
    outfiles.append(ofile)
    ens_extreme = varextreme_ens_np
    varsave = 'ens_extreme'
    ofile = os.path.join(dir_output, 'ens_extreme_{0}.nc'.format(name_outputs))
    save_n_2d_fields(lat_area, lon_area, ens_extreme, varsave,
                     varunits, ofile)
    outfiles.append(ofile)

    return outfiles
