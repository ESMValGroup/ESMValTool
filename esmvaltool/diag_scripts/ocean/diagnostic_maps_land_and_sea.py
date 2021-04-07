"""
Maps diagnostics
================

Diagnostic to produce images of a map with land and sea model.
These plost show latitude vs longitude and the cube value is used as the colour
scale.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_map:
      extract_levels:
        levels:  [100., ]
         scheme: linear_extrap
      climate_statistics:
        operator: mean

Note that this recipe may not function on machines with no access to the
internet, as cartopy may try to download the shapefiles. The solution to
this issue is the put the relevant cartopy shapefiles on a disk visible to your
machine, then link that path to ESMValTool via the `auxiliary_data_dir`
variable. The cartopy masking files can be downloaded from::

  https://www.naturalearthdata.com/downloads/

Here, cartopy uses the 1:10, physical coastlines and land files::

      110m_coastline.dbf  110m_coastline.shp  110m_coastline.shx
      110m_land.dbf  110m_land.shp  110m_land.shx

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import logging
import os
import sys
from itertools import product
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib import colors, colorbar

import iris
import iris.quickplot as qplt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import numpy as np
from esmvalcore.preprocessor._regrid import regrid

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


plot_pairs= {'pp':{'land': 'gpp', 'sea': 'intpp'},
            }
fx_mips = ['Ofx', 'fx', 'Lfx']


def longnameify(name):
    if isinstance(name, list):
        return ' '.join([longnameify(n) for n in name])
    if name == 'pp': return 'Primary Production'
    if name == 'intpp': return 'Integrated Primary Production'
    if name == 'gpp': return 'Gross Primary Production'
    if name == 'npp': return 'Net Primary Production'
    return name


def regrid_to_1x1(cube, scheme = 'linear'):
    """
    regrid a cube to a common 1x1 grid.
    """
    # regrid to a common grid:
    return regrid(cube, '1x1', scheme)


def cube_interesction(cube):
    central_longitude = 0. #W #-160.+3.5
    #central_latitude = 0.
    cube = regrid_to_1x1(cube)
    cube = cube.intersection(longitude=(central_longitude-180., central_longitude+180.))
    return cube

def detrended_contourf(cube):
    """
    make the detrended contouf plot
    """
    cmap = 'RdBu_r'
    drange = diagtools.get_cube_range_diff([cube, ])
    dlinspace = np.linspace(drange[0], drange[1], 22, endpoint=True)
    ticks = [t for t in np.linspace(dlinspace.min(),dlinspace.max(), 7)]
    cube = cube_interesction(cube)
    try:
        qplot = qplt.contourf(cube, dlinspace, cmap=cmap) # linewidth=0, rasterized=True,
        qplot.colorbar.set_ticks(ticks)
    except:
        print('Unable to plot cube:', cube)
        return False
    return True

def make_lon_cartopy_safe(lon):
    """
    CMOR states longitude should be between 0 and 360.
    but, Cartopy preferse -180->180.

    """
    lon[np.where(lon>=180.)] = lon[np.where(lon>=180.)]-360.
    return lon

def trended_contourf(cube, cmap='YlOrRd', zrange=[], drawcbar=True):
    """
    make the detrended contouf plot
    """
    cube = cube_interesction(cube)

    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    if lons.ndim ==2:
        ax = plt.axes(projection=ccrs.PlateCarree())

        # if lons.max()>180.:
        #     lons = make_lon_cartopy_safe(lons)

        # cube.data = np.ma.masked_where(lons>179.9, cube.data)
        print(cube.data)
        #assert 0
        #plt.contourf(lons, lats, cube.data, 18, cmap=cmap)
        im = plt.pcolormesh(lons, lats, cube.data, cmap=cmap)
        if drawcbar:
            plt.colorbar(orientation='horizontal')
        if zrange:
            im.set_clim(vmin=zrange[0], vmax=zrange[1])

        plt.gca().coastlines()
        if np.ma.is_masked(cube.data[0,0]):
            ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', edgecolor='face', facecolor='w'))
    else:
        try:
            qplt.contourf(cube, 18, cmap=cmap) # linewidth=0, rasterized=True,
        except:
            print('Unable to plot cube:', cube)
            return False
    return True


def single_pane_land_sea_pane(cfg,
        metadatas,
        fig,
        ax,
        land_cube=None,
        sea_cube=None,
        land_cmap = 'viridis',
        sea_cmap = 'Blues',
        ):
    """
    Generic way to plot single cube
    """
    #land = qplt.contourf(land_cube, 18, cmap=land_cmap) # linewidth=0, rasterized=True,
    #sea = qplt.contourf(sea_cube, 18, cmap=sea_cmap) # linewidth=0, rasterized=True,

    land = iris.plot.contourf(land_cube, 18, cmap=land_cmap) # linewidth=0, rasterized=True
    #cb = plt.colorbar(axp,ax=[ax],location='left')
    #landcbar = plt.colorbar(location='right')
    sea = iris.plot.contourf(sea_cube, 18, cmap=sea_cmap) # linewidth=0, rasterized=True,
    #seacbar = plt.colorbar(location='left')
    plt.gca().coastlines()

    return fig, ax, land, sea

def single_pane_land_sea_plot(
        cfg,
        metadatas,
        land_cube=None,
        sea_cube=None,
        plot_pair={},
        unique_keys = []
        ):
    """
    Try to plot both a land cube and a sea cube together.
    """
    sea_cube = cube_interesction(sea_cube)
    land_cube = cube_interesction(land_cube)
    fig = plt.figure()

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    suffix = '_'.join(unique_keys)+image_extention
    path = diagtools.folder([cfg['plot_dir'], 'single_land_sea_plots'])+ suffix

    #if os.path.exists(path): return

    fig.set_size_inches(10, 6)
    ax = plt.subplot(111, projection=ccrs.PlateCarree())
    fig, ax, land, sea = single_pane_land_sea_pane(cfg,
            metadatas,
            fig,
            ax,
            land_cube=land_cube,
            sea_cube=sea_cube,
            land_cmap = 'Greens', #'viridis',
            sea_cmap = 'Blues', # 'ocean_r'
            )
    divider = make_axes_locatable(ax)
    #pad_fraction = 0.5
    #aspect=20.
    #width = axes_size.AxesY(ax, aspect=1./aspect)
    #pad = axes_size.Fraction(pad_fraction, width)
   # cax = divider.append_axes("right", size=width, pad=pad)
    #caxL = divider.append_axes("left", size=width, pad=pad)
    #caxR = divider.append_axes("right", size=width, pad=pad)

    land_label = ', '.join([longnameify(land_cube.var_name), str(land_cube.units)])
    sea_label = ', '.join([longnameify(sea_cube.var_name), str(sea_cube.units)])

    landcbar = plt.colorbar(land, ax=[ax, ], location='left', label=land_label, shrink=0.55)
    seacbar = plt.colorbar(sea, ax=[ax, ], location='right', label=sea_label, shrink=0.55)

    # Add title to plot
    # if detrend:
    #     title = ' '.join([metadata['dataset'], key, '- detrended'])
    # else:
    #     title = ' '.join([metadata['dataset'], key, '- trend intact'])
    plt.title(longnameify(unique_keys))

    plt.gca().set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0, right = 0.80, left = 0.25, 
            hspace = 0, wspace = 0)
    plt.margins(0,0)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path) #, bbox_inches = 'tight',)
    plt.close()

def trended_pcolormesh(cube, ax, cmap='viridis', zrange=[], drawcbar=True):
    """
    make the detrended contouf plot
    """
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    if lons.ndim ==1:
        #lats,lons = np.meshgrid(lats.copy(),lons.copy())
        lons,lats = np.meshgrid(lons.copy(),lats.copy())

    if lons.ndim ==2:

        if lons.max()>180.:
            lons = make_lon_cartopy_safe(lons)
        cube.data = np.ma.masked_where(lons>179.9, cube.data)
        im = plt.pcolormesh(lons, lats, cube.data, cmap=cmap)

        if drawcbar:
            plt.colorbar(orientation='horizontal')

        if zrange:
            im.set_clim(vmin=zrange[0], vmax=zrange[1])

        plt.gca().coastlines()
        if np.ma.is_masked(cube.data[0,0]):
            ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', edgecolor='face', facecolor='w'))
    else:
        assert 0
        if lons.max()>180.:
            lons = make_lon_cartopy_safe(lons.copy())
        cube.data = np.ma.masked_where(lons>179.9, cube.data)
        im = plt.pcolormesh(lons, lats, cube.data, cmap=cmap)

        if drawcbar=='vertical':
            plt.colorbar(orientation='vertical')

        elif drawcbar in [True, 'horizontal']:
            plt.colorbar(orientation='horizontal')

        if zrange:
            im.set_clim(vmin=zrange[0], vmax=zrange[1])

        plt.gca().coastlines()
        if np.ma.is_masked(cube.data[0,0]):
            ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', edgecolor='face', facecolor='w'))

    return ax


def split_variable_groups(variable_group, debug=True ):
    """
    Split variable group into variable and experiment.
    """
    if debug:
        print('splitting:', variable_group)
    splits = variable_group.split('_')
    if len(splits) == 3:
        variable, exp, threshold = variable_group.split('_')
    elif len(splits) == 2:
        variable, exp = variable_group.split('_')
        threshold = ''
    else:
        print('What?',variable_group)
        assert 0
    #if variable == 'tas':
    #    variable = 'Surface Temperature'
    #if exp in [ 'fx', 'Ofx']:
    #    pass
    #else:
    #    exp = exp.upper()
    #    exp = ' '.join([exp[:3], exp[3], exp[4]+'.'+exp[5]])
    #if threshold == '15':
    #    threshold = '1.5'
    #if threshold: threshold += u'\N{DEGREE SIGN}'
    return variable, exp, threshold


def calc_ensemble_mean(cube_list):
    """
    Calculate the ensemble mean of a list of cubes.

    """
    if not isinstance(cube_list, list):
        print("calc_ensemble_mean: cube_list,",cube_list," is not a list, it's a ", type(cube_list))
        assert 0
    if len(cube_list) ==0:
        print("calc_ensemble_mean: cube_list,",cube_list," has no contents. ")
        assert 0
    if len(cube_list) ==1:
        return cube_list[0]

    cube_data = cube_list[0].data
    for c in cube_list[1:]:
        cube_data += c.data

    cube_data = cube_data/float(len(cube_list))
    ensemble_mean = cube_list[0]
    ensemble_mean.data = cube_data
    return ensemble_mean


def make_map_plots(
        cfg,
        metadata,
        cube,
        key,
        detrend,
        cmap='YlOrRd'
):
    """
    Make a simple map plot for an individual model.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary
    filename: str
        the preprocessed model file.
    detrend: bool
        key to show wherether the mean was subtracted.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if detrend:
        suffix = '_'.join([metadata['dataset'], key, 'detrended'])+image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Detrended', 'single_plots'])+ suffix
    else:
        suffix = '_'.join([metadata['dataset'], key])+image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Trend_intact', 'single_plots'])+suffix

    if os.path.exists(path): return

    # Making plots
    if detrend:
        success = detrended_contourf(cube)
    else:
        success = trended_contourf(cube, cmap=cmap)
    if not success:
        print('Failed to make figure:', path)
        print('key:', key)
        plt.close()
        return

    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    # Add title to plot
    if detrend:
        title = ' '.join([metadata['dataset'], key, '- detrended'])
    else:
        title = ' '.join([metadata['dataset'], key, '- trend intact'])

    plt.title(title)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()


def add_map_subplot(subplot, cube, nspace, title='', cmap=''):
    """
    Add a map subplot to the current pyplot figure.

    Parameters
    ----------
    subplot: int
        The matplotlib.pyplot subplot number. (ie 221)
    cube: iris.cube.Cube
        the iris cube to be plotted.
    nspace: numpy.array
        An array of the ticks of the colour part.
    title: str
        A string to set as the subplot title.
    cmap: str
        A string to describe the matplotlib colour map.

    """
    plt.subplot(subplot)
    cube = cube_interesction(cube)

    qplot = qplt.contourf(cube, nspace, linewidth=0,
                          cmap=plt.cm.get_cmap(cmap))
    qplot.colorbar.set_ticks([nspace.min(),
                              (nspace.max() + nspace.min()) / 2.,
                              nspace.max()])

    try: plt.gca().coastlines()
    except: pass
    plt.title(title)



def weighted_mean(cube, fx_fn):
    """
    Calculate the weighted mean.
    """
    fx_cube = iris.load_cube(fx_fn)
    grid_areas = ''
    fx_vars = ['areacella', 'areacello']
    fx_key = ''
    for fx_var in fx_vars:
        if fx_fn.find(fx_var) > -1:
            fx_key = fx_var

    #fx_cube = fx_cube.extract_strict(iris.Constraint(name=fx_key))

    return cube.collapsed(['latitude', 'longitude'],
                          iris.analysis.MEAN,
                          weights=fx_cube.data)


def make_four_pane_map_plot(
        cfg,
        cube_ssp,
        cube_hist,
        cube_anomaly,
        cube_detrended,
        key,
        plot_type,
):
    """
    Make a four pane map plot for an individual.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary
    filename: str
        the preprocessed model file.
    detrend: bool
        key to show wherether the mean was subtracted.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    suffix = '_'.join([plot_type, key]) + image_extention
    suffix = suffix.replace(' ', '')
    path = diagtools.folder([cfg['plot_dir'], 'four_pane_plots', plot_type]) + suffix

    fig = plt.figure()
    fig.set_size_inches(9, 6)

    # Create the cubes
    cube221 = cube_ssp
    cube222 = cube_hist
    cube223 = cube_anomaly
    cube224 = cube_detrended

    # create the z axis for plots 2, 3, 4.
    zrange1 = diagtools.get_cube_range([cube221, cube222])
    zrange3 = diagtools.get_cube_range_diff([cube223, ])
    zrange4 = diagtools.get_cube_range_diff([cube224, ])

    linspace1 = np.linspace(zrange1[0], zrange1[1], 12, endpoint=True)
    linspace3 = np.linspace(zrange3[0], zrange3[1], 12, endpoint=True)
    linspace4 = np.linspace(zrange4[0], zrange4[1], 12, endpoint=True)

    # Add the sub plots to the figure.
    add_map_subplot(221, cube221, linspace1, cmap='viridis',
                    title='Projection')
    add_map_subplot(222, cube222, linspace1, cmap='viridis',
                    title='Historical mean (1850-1900)')
    add_map_subplot(223, cube223, linspace3, cmap='RdBu_r',
                    title='Projection minus historical')
    add_map_subplot(224, cube224, linspace4, cmap='RdBu_r',
                    title='Detrended projection minus historical')

    # Add overall title
    title = key.replace('_', ' ')
    fig.suptitle(title, fontsize=14)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()


def make_ensemble_map_plots(
        cfg,
        cube,
        variable_group,
        detrend,
        key = '',
        zrange=[],
        cmap='YlOrRd'
):
    """
    Make a simple map plot for each variable_group.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if detrend:
        suffix = '_'.join(['variable_group_ensembles', variable_group, key, 'Detrended']) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Detrended', 'variable_group_ensembles']) + suffix
    else:
        suffix = '_'.join(['variable_group_ensembles', variable_group, key,'Trend_intact']) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Trend_intact', 'variable_group_ensembles']) + suffix
    path = path.replace('__', '_')
    if os.path.exists(path): return

    # Making plots
    if detrend:
        success = detrended_contourf(cube, cmap=cmap, zrange=zrange)
    else:
        success = trended_contourf(cube, cmap=cmap, zrange=zrange)
    if not success:
        print('Unable to make figure:', path)
        print('variable_group:', variable_group)
        plt.close()
        return

    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    # tas_ssp585_6
    variable, exp, threshold = split_variable_groups(variable_group)

    # Add title to plot
    if detrend:
        title = ' '.join([variable, '- ensemble mean of', exp, 'after', threshold, 'warming - detrended', key])
    else:
        title = ' '.join([variable, '- ensemble mean of', exp, 'after', threshold, 'warming - trend intact', key])
    title = title.replace('  ', ' ')
    plt.title(title)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()


def make_ensemble_map_plots_diff(
        cfg,
        cube,
        variable_group,
        detrend,
        temp,
        cmap='YlOrRd'
):
    """
    Make a simple difference map plot for each variable_group.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if detrend:
        suffix = '_'.join(['variable_group_ensembles_diff', variable_group, 'Detrended', temp]) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Detrended', 'variable_group_ensembles_diff']) + suffix
    else:
        suffix = '_'.join(['variable_group_ensembles_diff', variable_group, 'Trend_intact', temp]) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Trend_intact', 'variable_group_ensembles_diff']) + suffix
    if os.path.exists(path): return

    # Making plots
    if detrend:
        success = detrended_contourf(cube)
    else:
        success = trended_contourf(cube, cmap=cmap)
    if not success:
        print('Unable to make figure:', path)
        print('variable_group:', variable_group)
        plt.close()
        return

    #plt.gca().coastlines()

    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    # tas_ssp585_6
    variable, exp, threshold = split_variable_groups(variable_group)

    # Add title to plot
    if detrend:
        title = ' '.join([variable, '- ', exp, 'after', threshold, 'warming - difference against ssp126- detrended'])
    else:
        title = ' '.join([variable, '- ', exp, 'after', threshold, 'warming - diff against ssp126 trend intact'])
    plt.title(title)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()



def make_threshold_ensemble_map_plots(
        cfg,
        cube,
        variable,
        threshold,
        detrend,
        cmap='YlOrRd'
):
    """
    Make a simple map plot for an individual model.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if detrend:
        suffix = '_'.join(['variable_group_ensembles', 'Detrended', threshold, variable]) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Detrended', 'threshold_ensemble']) + suffix
    else:
        suffix = '_'.join(['variable_group_ensembles', 'Trend_intact', threshold, variable]) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Trend_intact', 'threshold_ensemble']) + suffix

    # Making plots for each layer
    #try: qplt.contourf(cube, 12, linewidth=0, rasterized=True, cmap=cmap)
    #except:
    #    print('Trying to make figure:', path)
    #    print('variable_group:', variable, threshold)
    #    print('Unable to plot cube:', cube)
    #    plt.close()
    #    return
    #    #qplt.contourf(cube, 12, linewidth=0, rasterized=True, cmap=cmap)

    # Making plots
    if detrend:
        success = detrended_contourf(cube)
    else:
        success = trended_contourf(cube,cmap=cmap)

    if not success:
        print('Unable to make figure:', path)
        print('variable_group:',threshold, variable)
        plt.close()
        return

    # plt.gca().coastlines()
    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    # Add title to plot
    if detrend:
        title = ' '.join(['Ensemble mean after', threshold, 'warming - detrended'])
    else:
        title = ' '.join(['Ensemble mean after', threshold, 'warming - trend intact'])

    plt.title(title)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()




def make_variable_group_mean_figure(cfg, variable_group_means, variable_group_details, threshold,
        variable='', units=''):
    """
    The idea here is to make a 6 pane figure showing:
    the historic mean (1850-1900) top left,
    Then a single pane for each scenario
    """

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    suffix = '_'.join(['variable_group_ensembles', 'Trend_intact', threshold, variable]) + image_extention
    path = diagtools.folder([cfg['plot_dir'], 'variable_group_ensembles'])+suffix

    scenarios ={
        'historical': 331,
        'ssp119':     332,
        'ssp126':     333,
        'ssp245':     334,
        'ssp370':     335,
        'ssp434':     336,
        'ssp585':     337,
        'ssp534-over': 338,}
    scenario_cubes = {}

    for variable_group, cube in variable_group_means.items():
        historical_group = variable_group[:variable_group.find('_')] +'_historical'
        if variable_group == historical_group:
            scenario_cubes['historical'] = cube
            continue

        if variable_group.find('_'+threshold) ==-1:
            continue
        for scen in scenarios.keys():
            if variable_group.find(scen)>-1:
                scenario_cubes[scen] = cube
    zrange = diagtools.get_cube_range([c  for s,c in scenario_cubes.items()])

    #linspace1 = np.linspace(zrange[0], zrange[1], 25, endpoint=True)

    fig = plt.figure()
    fig.set_size_inches(12, 7)

    # Add the sub plots to the figure.
    for scenario, cube in scenario_cubes.items():
        if scenario == 'historical':
            title = 'Historical (1850-1900)'
        elif scenario == 'ssp534-over':
            title = 'SSP534 Overshoot'
        else:
            title = scenario.upper()
#        add_map_subplot(scenarios[scenario], cube, linspace1, cmap='viridis',
#                    title=title)
        ax = plt.subplot(scenarios[scenario], projection=ccrs.PlateCarree())
        ax = trended_pcolormesh(cube, ax, cmap='viridis', zrange=zrange, drawcbar=False)
        plt.title(title)

    # Add overall title.
    suptitle = ' '.join([variable, 'Ensemble mean after', threshold, 'warming'])
    plt.suptitle(suptitle)

    # Adjsuyt white space
    plt.subplots_adjust(wspace=0.05, hspace=0.02, left=0, right=1, bottom=0, top=0.91)

    # add colorbar on side.
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.12, 0.05, 0.7])
    ##img = plt.pcolormesh(np.array([zrange, zrange]), cmap='viridis')
    #img = plt.imshow(np.array([zrange,zrange]), cmap='viridis')
    ##img.set_clim(vmin=zrange[0], vmax=zrange[1])
    #cbar_ax.set_visible(False)
    #plt.colorbar(img, cax=cbar_ax)
    norm = colors.Normalize(vmin=zrange[0],vmax=zrange[1])
    cb1  = colorbar.ColorbarBase(cbar_ax,cmap=plt.cm.get_cmap('viridis'),norm=norm,orientation='vertical')
    cb1.ax.set_ylabel(str(units))

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()


def make_variable_group_mean_anomaly_figure(cfg, hist_cube, variable_group_anomaly_means, threshold,
        variable='', units = ''):
    """
    The idea here is to make an 8 pane figure showing:
    the historic mean (1850-1900) top left,
    Then a single pane for each scenario showingthe anomaly against the historical.
    """

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    suffix = '_'.join(['variable_group_anomalies', 'Trend_intact', threshold, variable]) + image_extention
    path = diagtools.folder([cfg['plot_dir'], 'variable_group_anomalies'])+suffix

    scenarios ={
        'historical': 321,
        'ssp119':     322,
        'ssp126':     323,
        'ssp245':     324,
        'ssp370':     325,
        #'ssp434':     326,
        'ssp585':     326,}
        #'ssp534-over': 338,}
    scenario_cubes = {}

    for variable_group, cube in variable_group_anomaly_means.items():
        historical_group = variable_group[:variable_group.find('_')] +'_historical'
        if variable_group == historical_group:
            scenario_cubes['historical'] = cube
            continue

        if variable_group.find('_'+threshold) ==-1:
            continue
        for scen in scenarios.keys():
            if variable_group.find(scen)>-1:
                scenario_cubes[scen] = cube
    print(scenario_cubes)
    if not len(scenario_cubes.keys()): return

    zrange = diagtools.get_cube_range_diff([c  for s,c in scenario_cubes.items()])

    scenario_cubes['historical'] = hist_cube

    fig = plt.figure()
    fig.set_size_inches(12, 8)

    # Add the sub plots to the figure.
    for scenario, cube in scenario_cubes.items():
        if scenario == 'historical':
            title = 'Historical (1850-1900)'
            ax = plt.subplot(scenarios[scenario], projection=ccrs.PlateCarree())
            ax = trended_pcolormesh(cube, ax, cmap='viridis', drawcbar=True)
            plt.title(title)
            continue

        elif scenario == 'ssp534-over':
            title = 'SSP534 Overshoot'
        else:
            title = scenario.upper()
        ax = plt.subplot(scenarios[scenario], projection=ccrs.PlateCarree())
        ax = trended_pcolormesh(cube, ax, cmap='PRGn', zrange=zrange, drawcbar=False)
        plt.title(title)

    # Add overall title.
    suptitle = ' '.join([variable, 'Ensemble mean anomaly against hist. after', threshold, 'warming'])
    plt.suptitle(suptitle)

    # Adjsuyt white space
    plt.subplots_adjust(wspace=0.05, hspace=0.02, left=0, right=1, bottom=0, top=0.91)

    # add colorbar on side.
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.12, 0.05, 0.7])
    norm = colors.Normalize(vmin=zrange[0],vmax=zrange[1])
    cb1  = colorbar.ColorbarBase(cbar_ax,cmap=plt.cm.get_cmap('PRGn'), norm=norm, orientation='vertical')
    cb1.ax.set_ylabel(str(units))

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()



def make_variable_group_mean_anomaly_ssp126_figure(cfg, variable_group_anomaly_means, threshold,
       variable='', units = ''):

    """
    The idea here is to make an 8 pane figure showing:
    the historic mean (1850-1900) top left,
    Then a single pane for each scenario showingthe anomaly against the historical.
    """

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    suffix = '_'.join(['variable_group_anomalies_ssp126', 'Trend_intact', threshold, variable]) + image_extention
    path = diagtools.folder([cfg['plot_dir'], 'variable_group_anomalies_ssp126'])+suffix

    scenarios ={
    #    'historical': 331,
    #    'ssp119':     332,
        'ssp126':     221,
        'ssp245':     222,
        'ssp370':     223,
    #    'ssp434':     336,
        'ssp585':     224,
    #    'ssp534-over': 338,
    }
    scenario_cubes = {}

    for variable_group, cube in variable_group_anomaly_means.items():
        print(variable_group, threshold)
        if variable_group.find('_'+threshold) ==-1:
            continue
        for scen in scenarios.keys():
            print(scen)
            if variable_group.find(scen)>-1:
                scenario_cubes[scen] = cube
    print(scenario_cubes.keys())
    if 'ssp126' not in scenario_cubes.keys(): return

    ssp126_cube = scenario_cubes['ssp126']
    diff_cubes = {}
    for scen in scenarios.keys():
        if scen == 'ssp126': continue
        diff_cubes[scen] = scenario_cubes[scen] - ssp126_cube

    zrange = diagtools.get_cube_range_diff([c  for s,c in diff_cubes.items()])

    diff_cubes['ssp126'] = ssp126_cube

    fig = plt.figure()
    fig.set_size_inches(12, 7)

    # Add the sub plots to the figure.
    for scenario, cube in diff_cubes.items():
        if scenario == 'ssp126':
            title = 'SSP126'
            ax = plt.subplot(scenarios[scenario], projection=ccrs.PlateCarree())
            ax = trended_pcolormesh(cube, ax, cmap='viridis', drawcbar=True)
            plt.title(title)
            continue
        elif scenario == 'ssp534-over':
            title = 'SSP534 Overshoot'
        else:
            title = scenario.upper()
        ax = plt.subplot(scenarios[scenario], projection=ccrs.PlateCarree())
        ax = trended_pcolormesh(cube, ax, cmap='PRGn', zrange=zrange, drawcbar=False)
        plt.title(title+ ' minus SSP126')

    # Add overall title.
    suptitle = ' '.join([variable, 'Ensemble mean anomaly against SSP126 after', threshold, 'warming'])
    plt.suptitle(suptitle)

    # Adjsuyt white space
    plt.subplots_adjust(wspace=0.05, hspace=0.02, left=0, right=1, bottom=0, top=0.91)

    # add colorbar on side.
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.12, 0.05, 0.7])
    norm = colors.Normalize(vmin=zrange[0],vmax=zrange[1])
    cb1  = colorbar.ColorbarBase(cbar_ax,cmap=plt.cm.get_cmap('PRGn'), norm=norm, orientation='vertical')
    cb1.ax.set_ylabel(str(units))

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()




def make_gwt_map_four_plots(cfg, ):
    """
    Make plots

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    do_single_plots: bool
        Make individual plots for each dataset.

    """
    metadatas = diagtools.get_input_files(cfg)
    do_single_plots=True
    do_variable_group_plots=True
    do_threshold_plots=True

    #
    files_dict = {}
    short_names = set()
    standard_names = {}
    ensembles = set()
    exps = set()
    variable_groups = set()
    variables = set()
    mips = set()
    thresholds = {}
    #fx_fn = ''
    fx_fns = {}
    all_cubes = {}
    key_index = ['dataset', 'mip', 'exp', 'ensemble', 'short_name', 'variable_group']

    for fn, details in sorted(metadatas.items()):
        #print(fn, details.keys())
        short_name = details['short_name']
        if short_name not in ['areacella', 'areacello']:
            short_names.add(details['short_name'])
        standard_names[details['short_name']] = details['standard_name']
        ensembles.add(details['ensemble'])
        exps.add(details['exp'])
        variable_groups.add(details['variable_group'])
        unique_key = (details['variable_group'], details['ensemble'])
        mips.add(details['mip'])

        if details['mip'] in fx_mips:
            #fx_fn = fn
            fx_fns[details['short_name']] = fn
            print('Found FX MIP', fn)
        try:
            files_dict[unique_key].append(fn)
        except:
            files_dict[unique_key] = [fn, ]

        cube = iris.load_cube(fn)
        cube = diagtools.bgc_units(cube, details['short_name'])

        all_cubes[tuple([details[key] for key in key_index])]=  cube

    if len(fx_fns)=='':
        print('unable to find fx files:', fx_fns)
        assert 0
    print(files_dict.keys())

    # Make a plot for a single land-sea pair.
    cube_pairs = {}
    for (dataset, mip, exp, ensemble, short_name, variable_group), cube in all_cubes.items():
        for var_name, plot_pair in plot_pairs.items():
            if not do_single_plots:
                continue
            if mip in fx_mips:
                continue
            # avoid double plotting.
            if short_name != plot_pair['sea']:
                continue
            variable1, exp1, threshold = split_variable_groups(variable_group)
            print('index:', (dataset, mip, exp, ensemble, short_name, variable_group))
            print('split_variable_groups', variable1, exp1, threshold)

            sea_cube = cube
            if len(threshold):
                land_variable_group = '_'.join([plot_pair['land'], exp1, threshold])
            else:
                land_variable_group = '_'.join([plot_pair['land'], exp1])

            landmip = mip.replace('O', 'L')
            print('land stuff:', landmip, land_variable_group, plot_pair['land'])
            land_index = (dataset, landmip, exp, ensemble, plot_pair['land'], land_variable_group)
            print('land index', land_index)

            try: land_cube = all_cubes[land_index]
            except:
                for index in all_cubes.keys():
                   if plot_pair['land'] not in index: continue
                   print(index)
                assert 0
#            land_cube = all_cubes.get(land_index, None)
            #if land_cube is None: continue
            cube_pairs[(dataset, mip, exp, ensemble, threshold)] = {'sea': sea_cube, 'land': land_cube}
            
            single_pane_land_sea_plot(
                cfg,
                metadatas,
                land_cube=land_cube,
                sea_cube=sea_cube,
                plot_pair=plot_pair,
                unique_keys = [var_name, dataset, exp, ensemble, threshold]
            )

    assert 0

    # lets look at  minus the historical
    # hist_cubes = {variable_group:{} for variable_group in variable_groups}
    ssp_cubes = {variable_group:{} for variable_group in variable_groups}
    # anomaly_cubes = {variable_group:{} for variable_group in variable_groups}
    # detrended_anomaly_cubes = {variable_group:{} for variable_group in variable_groups}
    # viable_keys = {}

    hist_cubes_fns = {}
    # Calculate the anomaly for each ensemble/threshold combination
    for ensemble in sorted(ensembles):
        for variable_group in sorted(variable_groups):

                if not do_single_plots:
                    continue
                variable, exp, threshold = split_variable_groups(variable_group)
                # avoid double plotting.
                if variable != plot_pair['sea']:
                    continue

                # if (variable_group, ensemble) not in viable_keys:
                #     print("Not in viable keys:", ensemble, variable_group)
                #     continue
                sea_cube = ssp_cubes[variable_group][ensemble]
                land_variable_group = '_'.join([plot_pair['land'], exp, threshold])

                land_cube = ssp_cubes[land_variable_group][ensemble]

                single_pane_land_sea_plot(
                    cfg,
                    metadata,
                    land_cube=land_cube,
                    sea_cube=sea_cube,
                    pair_name=var_name,
                )
        hist_cubes_fns = {}
        # Calculate the anomaly for each ensemble/threshold combination
        for ensemble in sorted(ensembles):
            for variable_group in sorted(variable_groups):
            # guess historical group name:
                historical_group = variable_group[:variable_group.find('_')] +'_historical'

                if variable_group == historical_group:
                    continue
            assert 0
            fx_group = variable_group[:variable_group.find('_')] +'_fx'
            if variable_group == fx_group:
                continue

            print('Plotting:', ensemble, variable_group)
            variable, exp, threshold = split_variable_groups(variable_group)
            variables.add(variable)

            if (variable_group, ensemble) not in files_dict:
                continue
            fn = files_dict[(variable_group, ensemble)][0]
            fn_hist = files_dict[(historical_group, ensemble)][0]

            details = metadatas[fn]
            cube = iris.load_cube( fn)
            cube = diagtools.bgc_units(cube, details['short_name'])

            # Check is historial cube is loaded already
            if fn_hist in hist_cubes_fns:
                cube_hist = hist_cubes_fns[fn_hist]
            else:
                cube_hist =  iris.load_cube( fn_hist)
                cube_hist = diagtools.bgc_units(cube_hist, details['short_name'])

            print('fx files:', fx_fn, details['mip'])
            if details['mip'] in ['Omon', 'Oyr']:
                fx_fn = fx_fns['areacello']
            elif  details['mip'] in ['Lmon', 'Lyr', 'Amon', 'Ayr']:
                fx_fn = fx_fns['areacella']

            cube_anomaly = cube - cube_hist
            detrended_cube_anomaly = cube_anomaly - weighted_mean(cube_anomaly, fx_fn)

            ssp_cubes[variable_group][ensemble] = cube
            hist_cubes[variable_group][ensemble] = cube_hist
            anomaly_cubes[variable_group][ensemble] = cube_anomaly
            detrended_anomaly_cubes[variable_group][ensemble] = detrended_cube_anomaly
            viable_keys[(variable_group, ensemble)] = True
            try:
                thresholds[threshold].append([variable_group, ensemble])
            except:
                thresholds[threshold] = [[variable_group, ensemble], ]

    # single pane land-sea plots pairs:
    for var_name, plot_pair in plot_pairs.items():
        for ensemble in sorted(ensembles):
            for variable_group in sorted(variable_groups):
                if not do_single_plots:
                    continue
                variable, exp, threshold = split_variable_groups(variable_group)

                # avoid double plotting.
                if variable != plot_pair['sea']:
                    continue

                # if (variable_group, ensemble) not in viable_keys:
                #     print("Not in viable keys:", ensemble, variable_group)
                #     continue
                sea_cube = ssp_cubes[variable_group][ensemble]
                land_variable_group = '_'.join([plot_pair['land'], exp, threshold])

                land_cube = ssp_cubes[land_variable_group][ensemble]

                single_pane_land_sea_plot(
                    cfg,
                    metadata,
                    land_cube=land_cube,
                    sea_cube=sea_cube,
                    plot_pair=plot_pair,
                )

                # historical_group = variable_group[:variable_group.find('_')] +'_historical'
                # if variable_group == historical_group:
                #     continue
                #
                # cube_ssp = ssp_cubes[variable_group][ensemble]
                # cube_hist = hist_cubes[variable_group][ensemble]
                # cube_anomaly = anomaly_cubes[variable_group][ensemble]
                # cube_detrend = detrended_anomaly_cubes[variable_group][ensemble]
                #
                # key = ' '.join([standard_names[variable].title(),
                #                 exp,
                #                 'at'   , threshold,
                #                 '('+ ensemble+')'])
                # make_four_pane_map_plot(
                #     cfg,
                #     cube_ssp,
                #     cube_hist,
                #     cube_anomaly,
                #     cube_detrend,
                #     key,
                #     'single_plots',
                #     )
    assert 0


    # Ensemble mean for each variable_group:
    for variable_group in sorted(variable_groups):
        # check to make plots.
        if not do_variable_group_plots:
            continue
        if variable_group.split('_')[-1] == 'fx':
            continue
        # guess historical group name:
        historical_group = variable_group[:variable_group.find('_')] + '_historical'
        if variable_group == historical_group:
            continue

        # Calculate ensemble means
        output_cubes = []
        for cubes in [ssp_cubes, hist_cubes, anomaly_cubes, detrended_anomaly_cubes]:
            cube_list = []
            for vari, cube in sorted(cubes[variable_group].items()):
                print(variable_group, vari)
                cube_list.append(cube)
            if cube_list == []:
                continue
            ensemble_mean = calc_ensemble_mean(cube_list)
            output_cubes.append(ensemble_mean)

        variable, exp, threshold = split_variable_groups(variable_group)
        key = ' '.join([standard_names[variable].title(),
                        'in', exp,
                        'at', threshold ,
                        'ensemble mean'
                        ])

        if len(output_cubes) != 4:
            print('Problem with ', variable_group, 'four plots.', historical_group)
            assert 0

        make_four_pane_map_plot(
            cfg,
            output_cubes[0],
            output_cubes[1],
            output_cubes[2],
            output_cubes[3],
            key,
            'variable_group_ensembles',
            )

    # Ensemble mean for each threshold:
    for threshold, paths in sorted(thresholds.items()):
        if not do_threshold_plots: continue
        output_cubes = []
        for cubes in [ssp_cubes, hist_cubes, anomaly_cubes, detrended_anomaly_cubes]:
            cube_list = []
            for (variable_group, ensemble) in viable_keys.keys():
                var_g, exp, thres = split_variable_groups(variable_group)
                if thres != threshold:
                    continue
                cube_list.append(cubes[variable_group][ensemble])
            ensemble_mean = calc_ensemble_mean(cube_list)
            output_cubes.append(ensemble_mean)

        key = ' '.join([short_name, threshold])
        key = ' '.join([standard_names[variable].title(),
                        'at', threshold,
                        'ensemble mean'
                        ])

        make_four_pane_map_plot(
            cfg,
            output_cubes[0],
            output_cubes[1],
            output_cubes[2],
            output_cubes[3],
            key,
            'threshold_ensemble',
            )


def make_gwt_map_plots(cfg, detrend = True,):
    """
    Make plots

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    detrend: bool
        Subtract the mean from the figure.
    do_single_plots: bool
        Make individual plots for each dataset.

    """
    metadatas = diagtools.get_input_files(cfg)

    do_single_plots=True
    do_single_anomaly_plots = True
    do_variable_group_plots=True
    do_variable_group_diff = True
    do_threshold_plots=True
    do_variable_group_mean_plots=True
    do_variable_group_anomaly_mean_plots=True

    #print('\n', cfg.keys())

    files_dict = {}
    short_names = set()
    ensembles = set()
    exps = set()
    variable_groups = set()
    # note that variable groups need to have the format: field_scenario_threshold
    # ie: fgco2_ssp126_1.5, tas_ssp585_5.0
    variables = set()
    thresholds = {}
    units = {}
    fx_fns = {}


    for fn, details in sorted(metadatas.items()):
        #print(fn, details.keys())
        short_names.add(details['short_name'])
        ensembles.add(details['ensemble'])
        exps.add(details['exp'])
        variable_groups.add(details['variable_group'])
        units[details['short_name']] = details['units']
        unique_key = (details['variable_group'], details['ensemble'])
        mip = details['mip']

        if details['mip'] in ['fx', 'Ofx']:
            fx_fns[details['short_name']] = fn
            print('Found FX MIP', fn, details['short_name'])

        try:
            files_dict[unique_key].append(fn)
        except:
            files_dict[unique_key] = [fn, ]

    if fx_fns is None:
        print('unable to find fx file.')
        assert 0

    print(files_dict.keys())

    # lets look at  minus the historical
    variable_group_cubes = {variable_group:{} for variable_group in variable_groups}
    variable_group_details = {variable_group:{} for variable_group in variable_groups}

    # Load and set up each ensemble member:
    for ensemble in sorted(ensembles):
        for variable_group in sorted(variable_groups):
            # guess historical group name:
            historical_group = variable_group[:variable_group.find('_')] +'_historical'
            if variable_group == historical_group:
                variable, exp = variable_group.split('_')
                threshold = None

            elif variable_group.split('_')[-1] == 'fx':
                variable, exp = variable_group.split('_')
                threshold = None
            else:
                variable, exp, threshold = split_variable_groups(variable_group)
                variables.add(variable)

            if (variable_group, ensemble) not in files_dict:
                continue
            fn = files_dict[(variable_group, ensemble)][0]
            fn_hist = files_dict[(historical_group, ensemble)][0]

            details = metadatas[fn]
            cube = iris.load_cube( fn)
            if variable_group.split('_')[-1] != 'fx':
                cube = diagtools.bgc_units(cube, details['short_name'])

            print('Loaded Cube:',ensemble,variable_group, [variable, exp, threshold])
            variable_group_cubes[variable_group][ensemble] = cube
            variable_group_details[variable_group][ensemble] = details

            key = variable_group.replace('_',' ') + ' '+ensemble
            if do_single_plots:
                make_map_plots(cfg, details, cube, key, detrend)

    print(variable_group_cubes.keys())

    # Calculate the anomaly for each ensemble/threshold combination
    anomaly_cubes = {variable_group:{} for variable_group in variable_groups}
    for ensemble in sorted(ensembles):
        for variable_group in sorted(variable_groups):
            # guess historical group name:
            historical_group = variable_group[:variable_group.find('_')] +'_historical'
            if variable_group == historical_group:
                continue

            print('Plotting anomaly:', ensemble, variable_group)

            if variable_group.split('_')[-1] == 'fx':
                continue

            variable, exp, threshold = split_variable_groups(variable_group)
            variables.add(variable)

            if (variable_group, ensemble) not in files_dict:
                continue

            # get details.
            details = variable_group_details[variable_group][ensemble]

            # get cube.
            cube = variable_group_cubes[variable_group][ensemble]

            # get historial cube.
            cube_hist = variable_group_cubes[historical_group][ensemble]
            cube_anomaly = cube - cube_hist

            if detrend:
                cube_anomaly = cube_anomaly - weighted_mean(cube_anomaly, fx_fn)

            anomaly_cubes[variable_group][ensemble] = cube_anomaly
            try:
                thresholds[threshold].append([variable_group, ensemble])
            except:
                thresholds[threshold] = [[variable_group, ensemble], ]
            key = variable_group.replace('_',' ') + ' '+ensemble +' anomaly'

            # Produce a plot of the anomaly.
            if do_single_anomaly_plots:
                make_map_plots(cfg, details, cube_anomaly, key, detrend)

    ####
    # Create variable_group_means.
    #    This is a mean of all ensemble memnbers of a given scenario at a given threshold.
    variable_group_means = {}
    for variable_group in sorted(variable_groups):
        cube_list = []
        for vari, cube in sorted(variable_group_cubes[variable_group].items()):
            print(variable_group, vari)
            cube_list.append(cube)
        if not len(cube_list):
            continue
        if variable_group.split('_')[-1] == 'fx':
            continue

        variable_group_mean = calc_ensemble_mean(cube_list)
        variable_group_means[variable_group] = variable_group_mean
        if do_variable_group_mean_plots:
            make_ensemble_map_plots(cfg, variable_group_mean, variable_group, detrend, key='', cmap='viridis')

    ####
    # Create variable_group_anomaLty_means.
    variable_group_anomaly_means={}
    for variable_group in sorted(variable_groups):
        cube_list = []
        if variable_group not in anomaly_cubes:
            continue
        for vari, cube in sorted(anomaly_cubes[variable_group].items()):
            print(variable_group, vari)
            cube_list.append(cube)
        if not len(cube_list):continue
        cube = calc_ensemble_mean(cube_list)
        variable_group_anomaly_means[variable_group] = cube

        if do_variable_group_anomaly_mean_plots:
            zrange = diagtools.get_cube_range_diff([cube, ])
            make_ensemble_map_plots(cfg, cube, variable_group, detrend, key='anomaly',cmap='PRGn', zrange=zrange)

    # All variable_group_means in a single figure:
    for threshold in ['1.5', '2.0', '3.0', '4.0', '5.0']:
        for variable in sorted(variables):
            make_variable_group_mean_figure(cfg, variable_group_means, variable_group_details, threshold,
                variable=variable, units = units[variable])

    # All variable_group_means anomalies in a single figure:
    for threshold in ['1.5', '2.0', '3.0', '4.0', '5.0']:
        for variable in sorted(variables):
            hist_cube = variable_group_means[variable +'_historical']
            make_variable_group_mean_anomaly_figure(cfg, hist_cube, variable_group_anomaly_means, threshold,
                variable=variable, units = units[variable])

    # 4 Pane figure showing anomaly against ssp126.
    for threshold in ['1.5', '2.0', ]:#'3.0', '4.0', '5.0']:
        for variable in sorted(variables):
            make_variable_group_mean_anomaly_ssp126_figure(cfg, variable_group_anomaly_means, threshold,
                variable=variable, units = units[variable])



    # Plot Ensemble mean for each variable_group, for each sceranrio and threshold,
    # include hist and exclude fx data.
    for variable_group, variable_group_mean in variable_group_means.items():
        # check to make plots.
        if not do_variable_group_plots:
            continue

        # guess historical group name:
        #historical_group = variable_group[:variable_group.find('_')] + '_historical'
        #if variable_group == historical_group:
        #    continue

        fx_group = variable_group[:variable_group.find('_')] +'_fx'
        if variable_group == fx_group:
            continue

        make_ensemble_map_plots(cfg, variable_group_mean, variable_group, detrend)

    # # difference bwtween each scenario and diff_against.
    # diff_against='ssp126'
    # for temp in ['2.0', '4.0', ]: # ['1.5', '2.0', '3.0', '4.0', '5.0']:
    #     for variable_group, variable_group_mean in variable_group_means.items():
    #         # check to make plots.
    #         if not do_variable_group_diff:
    #             continue
    #
    #         if variable_group.find('_'+temp) ==-1: continue
    #
    #     ssp126_group = variable_group[:variable_group.find('_')] + '_'+diff_against+'_'+temp
    #     print('diff_group:',ssp126_group, diff_against, temp)
    #     if variable_group == ssp126_group:
    #         continue
    #
    #     if ssp126_group not in variable_group_means.keys(): continue
    #     print('variable_group diff:', variable_group, ssp126_group, temp)
    #     diff_cube = variable_group_mean - variable_group_means[ssp126_group]
    #     make_ensemble_map_plots_diff(cfg, diff_cube, variable_group, detrend, temp)

    # Ensemble mean for each threshold:
    for variable in sorted(variables):
        if not do_threshold_plots: continue
        for threshold, paths in sorted(thresholds.items()):
            cube_list = []
            for [variable_group, ensemble] in sorted(paths):
                var_g, exp, threshold = split_variable_groups(variable_group)
                if var_g != variable:
                    continue
                cube_list.append(anomaly_cubes[variable_group][ensemble])

            if cube_list == []:
                continue
            ensemble_mean = calc_ensemble_mean(cube_list)
            make_threshold_ensemble_map_plots(cfg, ensemble_mean, variable, threshold, detrend, )


def main(cfg):
    """
    Load the config file, and send it to the plot makers.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    cartopy.config['data_dir'] = cfg['auxiliary_data_dir']

    make_gwt_map_four_plots(cfg)

    #for detrend in [False,]: # True ]:
    #    make_gwt_map_plots(cfg, detrend)


    # for index, metadata_filename in enumerate(cfg['input_files']):
    #     logger.info(
    #         'metadata filename:\t%s',
    #         metadata_filename,
    #     )
    #
    #     metadatas = diagtools.get_input_files(cfg, index=index)
    #     #thresholds = diagtools.load_thresholds(cfg, metadatas)
    #
    #
    #     # for filename in sorted(metadatas.keys()):
    #     #
    #     #     logger.info('-----------------')
    #     #     logger.info(
    #     #         'model filenames:\t%s',
    #     #         filename,
    #     #     )
    #     #
    #     #     ######
    #     #     # Contour maps of individual model
    #     #     if thresholds:
    #     #         make_map_contour(cfg, metadatas[filename], filename)
    #     #
    #     #     ######
        #     # Maps of individual model
        #     make_map_plots(cfg, metadatas[filename], filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
