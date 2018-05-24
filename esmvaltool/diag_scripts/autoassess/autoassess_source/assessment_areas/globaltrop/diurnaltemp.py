'''
Module docstring
'''

import os

import matplotlib as mpl
import matplotlib.cm as mpl_cm
import matplotlib.colors as mcol
import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cf_units
import iris
import iris.coord_categorisation
import iris.plot as iplt

from ..auto_assess_deprecated.loaddata import load_run_ss
from method_constraint import MethodConstraint
from plotting import flddiff, segment2list

SEASONS = ('ann', 'djf', 'mam', 'jja', 'son')


def _calc_means(run, proc_code):
    '''
    Routine to calculate annual and seasonal means from daily data. Uses
    categorisation to calculate seasonal means. Can do similar from monthly
    data.

    Arguments:
        datadir - Directory containing data for all of a single model run.
        proc_code - PP processing code used as extra filter in load data module
        dates - dictionary with start and end dates for filtering model data

    Returns:
        means - dictionary of cubes keyed by meaning period ('ann', 'djf', ...)
    '''

    # Extract daily min/max 1.5m temperatures and calculate seasonal
    # and annual means
    # m01s03i236
    temp = load_run_ss(run, 'daily', 'air_temperature', lbproc=proc_code)
    annual = temp.collapsed('time', iris.analysis.MEAN)
    iris.coord_categorisation.add_season(temp, 'time', name='clim_season')
    season = temp.aggregated_by('clim_season', iris.analysis.MEAN)
    djf = season.extract(iris.Constraint(clim_season='djf'))
    mam = season.extract(iris.Constraint(clim_season='mam'))
    jja = season.extract(iris.Constraint(clim_season='jja'))
    son = season.extract(iris.Constraint(clim_season='son'))

    return dict(ann=annual, djf=djf, mam=mam, jja=jja, son=son)


def plot_tmin_tmax(runs):
    '''
    Routine to plot seasonal and annual means of daily minimum and maximum
    temperatures.

    This is a type 2 multi-function routine that returns no metrics

    Arguments:
        runs - list of run dictionaries.  Each dictionary contains
               metadata for a single model run.  The first dictionary
               in this list is the control experiment.
               (see auto_assess.model_run.py for description of
               the contents of this dictionary)

    Returns:
        doesn't return any objects - it only writes image files to the
        current working dir
    '''

    run_cntl = runs[0]
    run_expts = runs[1:]

    try:

        max_cons = MethodConstraint('time', method='maximum')
        min_cons = MethodConstraint('time', method='minimum')
        cons = (max_cons, min_cons)
        clim = 'HadGHCND (1982 - 2001)'
        climdir = os.path.join(run_cntl['clim_root'], 'HadGHCND')

        for run_expt in run_expts:

            cntl = run_cntl.id
            expt = run_expt.id
            ids = dict(expt=expt.upper(), cntl=cntl.upper(), clim=clim)
            vs_runids = '{0}_v_{1}'.format(expt, cntl)
            f_file = '{}_1.5m_daily_temperature'.format(vs_runids)

            for season in SEASONS:
                clim_file = 'tmin_tmax_{}_1982_2001.pp'.format(season)
                season_file = 'tmin_tmax_{0}_{1}.nc'.format(season,
                                                            run_expt.period)

                clim_file = os.path.join(climdir, clim_file)

                (max_clim, min_clim) = iris.load_cubes(clim_file, cons)
                # XXX Bug workaround:
                # Coordinate names are unicode when loaded from netcdf, and thus
                # incompatible to the cube coordinates from PP files
                max_clim.coord(u'longitude').rename('longitude')
                max_clim.coord(u'latitude').rename('latitude')
                min_clim.coord(u'longitude').rename('longitude')
                min_clim.coord(u'latitude').rename('latitude')

                # extract meaned daily min and max 1.5m temperatures
                min_expt = _calc_means(run_expt, 4096)[season]
                max_expt = _calc_means(run_expt, 8192)[season]
                min_cntl = _calc_means(run_cntl, 4096)[season]
                max_cntl = _calc_means(run_cntl, 8192)[season]

                # Deal with deficencies in climatology metadata
                # Should this be done with a callback function?
                # Will I need callbacks for all climatology data sets!
                max_clim.units = cf_units.Unit('celsius')
                min_clim.units = cf_units.Unit('celsius')
                max_clim.convert_units('K')
                min_clim.convert_units('K')

                outfile = '{0}_{{}}_{1}.png'.format(f_file, season)

                # Calculate temperature ranges
                rng_expt = max_expt - min_expt
                rng_cntl = max_cntl - min_cntl
                rng_clim = max_clim - min_clim

                # Plot data
                _plot_temps(min_expt, min_cntl, min_clim, ids, 'minimum',
                            season, outfile)
                _plot_temps(max_expt, max_cntl, max_clim, ids, 'maximum',
                            season, outfile)
                _plot_temps(rng_expt, rng_cntl, rng_clim, ids, 'range',
                            season, outfile)
    except (IOError, iris.exceptions.ConstraintMismatchError) as e:
        print e
        print 'No min/max 1.5 temperature diagnostics found, so skipping ...'


def _make_title(plot, field, expt1, expt2=None):
    if expt2 is None:
        return '{0}) {1}\n{2}'.format(plot, field, expt1)
    else:
        return '{0}) {1}\n{2} minus {3}'.format(plot, field, expt1, expt2)


def _plot_temps(expt, cntl, clim, ids, stat, season, outfile):
    if stat == 'range':
        trange = True
    else:
        trange = False
    field = '1.5m temperature {0} for {1}'.format(stat, season)
    fig = plt.figure()
    atitle = _make_title('a', field, expt1=ids['expt'])
    btitle = _make_title('b', field, expt1=ids['expt'], expt2=ids['cntl'])
    ctitle = _make_title('c', field, expt1=ids['cntl'], expt2=ids['clim'])
    dtitle = _make_title('d', field, expt1=ids['expt'], expt2=ids['clim'])
    _plot_temp(expt, trange=trange, subplot=221, title=atitle)
    _plot_temp(expt, diff=cntl, trange=trange, subplot=222, title=btitle)
    _plot_temp(cntl, diff=clim, trange=trange, subplot=223, title=ctitle)
    _plot_temp(expt, diff=clim, trange=trange, subplot=224, title=dtitle)
    fig.savefig(outfile.format(stat))
    plt.close()


def _plot_temp(cube, diff=None, trange=False, subplot=111, title=None):
    if diff:
        cmap = mpl_cm.get_cmap('brewer_RdBu_11')
        levels = np.arange(-8, 9, 1)
        plotcube = flddiff(cube, diff)
    else:
        cmap = mpl_cm.get_cmap('jet')
        if trange:
            levels = np.arange(0, 25, 1)
        else:
            levels = np.arange(190, 315, 5)
        plotcube = cube.copy()
    colormap = segment2list(cmap, levels.size)
    normalisation = mcol.BoundaryNorm(levels, levels.size-1)

    axes = plt.subplot(subplot, projection=ccrs.PlateCarree())
    mesh = iplt.pcolormesh(plotcube, cmap=colormap, norm=normalisation)
    cbar = plt.colorbar(mesh, orientation='horizontal',
                        extend='both', drawedges=True)
    cbar.ax.tick_params(labelsize='xx-small')
    axes.set_global()
    axes.add_feature(cfeature.COASTLINE)
    glines = axes.gridlines(draw_labels=True)
    glines.xlines = None
    glines.ylines = None
    glines.xlabels_top = False
    glines.ylabels_right = False
    glines.xformatter = LONGITUDE_FORMATTER
    glines.yformatter = LATITUDE_FORMATTER
    glines.xlabel_style = {'size': 'xx-small'}
    glines.ylabel_style = {'size': 'xx-small'}
    if title:
        axes.set_title(title, fontdict={'fontsize': 'small'})
