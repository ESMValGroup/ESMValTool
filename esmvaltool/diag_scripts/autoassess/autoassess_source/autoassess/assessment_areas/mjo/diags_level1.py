import numpy as np
import iris
import os
import mjo_utils as mu
import mjo_plots as mjp
import iris.quickplot as qplt
from iris.time import PartialDateTime


def diagnos_level1(data, outplotdir, runid, tmp_dir):

    with iris.FUTURE.context(cell_datetime_objects=True):
        start_date = data.coord('time').cell(0).point.strftime('%Y-%m-%d')
        end_date = data.coord('time').cell(-1).point.strftime('%Y-%m-%d')

    varname = data.name()

    if varname == 'x_wind':
        assert len(data.coord('pressure').points) == 1
        pressure_level = data.coord('pressure').points[0]
        if pressure_level == 850:
            varname = 'x_wind_850hPa'
        elif pressure_level == 200:
            varname = 'x_wind_200hPa'

    assert varname in [
        'toa_outgoing_longwave_flux',
        'precipitation_flux',
        'x_wind_850hPa',
        'x_wind_200hPa'
    ]

    print 'Starting Level 1 diagnostics for %s...' % varname

    # Define contour levels for plots
    if varname == 'toa_outgoing_longwave_flux':
        datamean_clevels = range(180, 280, 10)
        datamean_colorReverse = True
        datavar_clevels = range(200, 2800, 200)
        datavar_filt_clevels = range(200, 650, 50)
    elif varname == 'precipitation_flux':
        datamean_clevels = range(0, 16, 1)
        datamean_colorReverse = False
        datavar_clevels = range(0, 500, 50)
        datavar_filt_clevels = range(4, 44, 4)
    elif varname == 'x_wind_850hPa':
        datamean_clevels = range(-16, 18, 2)
        datamean_colorReverse = False
        datavar_clevels = range(0, 110, 10)
        datavar_filt_clevels = range(3, 13, 1)
    elif varname == 'x_wind_200hPa':
        datamean_clevels = range(-20, 44, 4)
        datamean_colorReverse = False
        datavar_clevels = range(0, 220, 20)
        datavar_filt_clevels = range(10, 70, 10)

    # mean plots
    summerMean = mu.SummerExtract(data)
    summerMean = summerMean.collapsed('time', iris.analysis.MEAN)

    winterMean = mu.WinterExtract(data)
    winterMean = winterMean.collapsed('time', iris.analysis.MEAN)

    ####################################################
    print '1. Plot Summer/winter mean'
    title = runid + " " + varname + " mean summer"
    title += '\n' + start_date + ' - ' + end_date
    forename = runid + "_" + varname + "_Mean_summer"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(summerMean, title, datamean_clevels, figname, colorReverse=datamean_colorReverse)
    # save the plotted field to netcdf file
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(summerMean, ncname)


    title = runid + " " + varname + " mean winter"
    title += '\n' + start_date + ' - ' + end_date
    forename = runid + "_" + varname + "_Mean_winter"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(winterMean, title, datamean_clevels, figname, colorReverse=datamean_colorReverse)
    # save the plotted field to netcdf file
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(winterMean, ncname)
    ####################################################


    # variance plots
    # Compute anomalies by removing a climatological annual cycle as in NCL
    harmonics = 8
    # Climatology
    clim = mu.clmDayTLL(data)
    # Smooth climatology
    clim_sm = mu.smthClmDayTLL(clim, harmonics)
    # Compute anomalies
    data = mu.calcDayAnomTLL(data, clim_sm)

    # Compute variances for both seasons
    summerVar = mu.SummerExtract(data)
    summerVar = summerVar.collapsed('time', iris.analysis.VARIANCE)

    winterVar = mu.WinterExtract(data)
    winterVar = winterVar.collapsed('time', iris.analysis.VARIANCE)

    print '2. Plot Summer/winter variance'
    title = runid + " " + varname + " Var summer"
    title += '\n' + start_date + ' - ' + end_date
    forename = runid + "_" + varname + "_Var_summer"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(summerVar, title, datavar_clevels, figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(summerVar, ncname)

    title = runid + " " + varname + " Var winter"
    title += '\n' + start_date + ' - ' + end_date
    forename = runid + "_" + varname + "_Var_winter"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(winterVar, title, datavar_clevels, figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(winterVar, ncname)
    ####################################################

    # filtered variance plots
    print 'Filtering the series...'
    # Code below is a exact replica of the Fortran version
    # using a set of supplied filter coefficients
    datafilt = mu.Filter(data)


    print '3. Plot Summer/winter Filtered variance'
    summerFiltVar = mu.SummerExtract(datafilt)
    summerFiltVar = summerFiltVar.collapsed('time', iris.analysis.VARIANCE)

    winterFiltVar = mu.WinterExtract(datafilt)
    winterFiltVar = winterFiltVar.collapsed('time', iris.analysis.VARIANCE)

    original_varname = varname
    varname = varname + '_Filt'
    title = runid + " " + varname + " Var summer"
    title += '\n' + start_date + ' - ' + end_date
    forename = runid + "_" + varname + "_Var_summer"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(summerFiltVar, title, datavar_filt_clevels, figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(summerFiltVar, ncname)


    title = runid + " " + varname + " Var winter"
    title += '\n' + start_date + ' - ' + end_date
    forename = runid + "_" + varname + "_Var_winter"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(winterFiltVar, title, datavar_filt_clevels, figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(winterFiltVar, ncname)

    ####################################################
    print '4. Plot ratio of filtered variance to unfiltered variance'
    summerRatio = summerFiltVar / summerVar * 100.
    winterRatio = winterFiltVar / winterVar * 100.

    title = runid + "_" + varname + " VarRatio summer"
    title += '\n' + start_date + ' - ' + end_date
    forename = runid + "_" + varname + "_VarRatio_summer"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(summerRatio, title, range(10, 60, 5), figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(summerRatio, ncname)

    title = runid + "_" + varname + " VarRatio winter"
    title += '\n' + start_date + ' - ' + end_date
    forename = runid + "_" + varname + "_VarRatio_winter"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(winterRatio, title, range(10, 60, 5), figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(winterRatio, ncname)
    ####################################################
    print '*' * 50
    print original_varname + ' Level 1 diagnostics completed.'
    print '*' * 50

    return
