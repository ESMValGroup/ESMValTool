"""
Maps diagnostics
================

Diagnostic to produce images of a map with coastlines from a cube.
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
import itertools
import logging
import os
import sys
import numpy as np

import cartopy
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

from esmvalcore.preprocessor._time import climate_statistics,regrid_time

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def make_mean_of_cube_list(cube_list, operation='mean'):
    """
    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.
    """
    # Fix empty times
    full_times = {}
    times = []
#   for cube in cube_list:
#       # make time coords uniform:
#       #cube.coord('time').long_name='Time axis'
#       #cube.coord('time').attributes={'time_origin': '1950-01-01 00:00:00'}
#       #times.append(cube.coord('time').points)
#       cube = regrid_time(cube, 'mon')
#        for time in cube.coord('time').points:
#            print(cube.name, time, cube.coord('time').units)
#            try:
#                full_times[time] += 1
#            except:
#                full_times[time] = 1
#
#    for t, v in sorted(full_times.items()):
#        if v != len(cube_list):
#            print('FAIL', t, v, '!=', len(cube_list),'\nfull times:',  full_times)
#            assert 0

    cube_mean=cube_list[0]
    #try: iris.coord_categorisation.add_year(cube_mean, 'time')
    #except: pass
    #try: iris.coord_categorisation.add_month(cube_mean, 'time')
    #except: pass
    Cube_mean = regrid_time(cube_mean, 'mon')
    try: cube_mean.remove_coord('year')
    except: pass
    try: cube_mean.remove_coord('day_of_year')
    except: pass

    #cube.remove_coord('Year')
    try: model_name = cube_mean.metadata[4]['source_id']
    except: model_name = ''
    print(model_name,  cube_mean.coord('time'))

    for i, cube in enumerate(cube_list[1:]):
        try: cube.remove_coord('year')
        except: pass
        try: cube_mean.remove_coord('day_of_year')
        except: pass
        cube = regrid_time(cube, 'mon')
        try: model_name = cube_mean.metadata[4]['source_id']
        except: model_name = ''
        print(i, model_name, cube.coord('time'))
        cube_mean+=cube
        #print(cube_mean.coord('time'), cube.coord('time'))
    cube_mean = cube_mean/ float(len(cube_list))
    return cube_mean



def make_mean_of_cube_list_notime(cube_list):
    """
    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.
    """
    cube_mean=cube_list[0].copy()
    print('cacl mean:', cube_mean.units)
    try: cube_mean.remove_coord('year')
    except: pass

    for i, cube in enumerate(cube_list[1:]):
        try: cube.remove_coord('year')
        except: pass
        print('cacl mean:',i, cube.units)

        cube_mean+=cube

    cube_mean = cube_mean/ float(len(cube_list))
    return cube_mean


def make_stddev_of_cube_list(mean_timed_cube):
    """
    This takes the standard deviation of the means.
    What Ana wanted was the mean of the standard deviations (along the time axis).

    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.
    """
    print('std_dev calc mean cube:', mean_cube.units, mean_cube.data.max(),  mean_cube.data.min())
    if mean_cube is None:
        mean_cube = make_mean_of_cube_list_notime(cube_list)

    if count_cube is None:
        count_cube = make_count_of_cube_list_notime(cube_list)

    for i, cube in enumerate(cube_list):
        cube = cube.copy()
        try: cube.remove_coord('year')
        except: pass
        print('std_dev calc', i, cube.units, cube.data.max(),  cube.data.min())
        cube.data = np.ma.power((cube.data - mean_cube.data), 2.)
        cube_var_list.append(cube)

    cube_var = make_mean_of_cube_list_notime(cube_var_list)

    cube_var.data = np.ma.sqrt(cube_var.data/count_cube.data)

    return cube_var

def make_stddev_of_cube_list_notime(cube_list, mean_cube=None, count_cube=None):
    """
    This takes the standard deviation of the means.
    What Ana wanted was the mean of the standard deviations (along the time axis).

    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.
    """
    cube_var_list=[]
    print('std_dev calc mean cube:', mean_cube.units, mean_cube.data.max(),  mean_cube.data.min())
    if mean_cube is None:
        mean_cube = make_mean_of_cube_list_notime(cube_list)

    if count_cube is None:
        count_cube = make_count_of_cube_list_notime(cube_list)

    for i, cube in enumerate(cube_list):
        cube = cube.copy()
        try: cube.remove_coord('year')
        except: pass
        print('std_dev calc', i, cube.units, cube.data.max(),  cube.data.min())
        cube.data = np.ma.power((cube.data - mean_cube.data), 2.)
        cube_var_list.append(cube)

    cube_var = make_mean_of_cube_list_notime(cube_var_list)

    cube_var.data = np.ma.sqrt(cube_var.data/count_cube.data)

    return cube_var


def make_count_of_cube_list(cube_list,):
    """
    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.

    Input cube has time, but output does not.
    """
    # Fix empty times
    cube_count = cube_list[0][0].copy()

    #cube_count.data = (~np.ma.masked_invalid(cube_count.data)).mask.astype(float)
    cube_count.data = np.ma.masked_invalid(cube_count.data)
    cube_count.data = np.ma.clip(cube_count.data, 1.,1.)

    if len(cube_list)==1:
        return cube_count

    for i, cube in enumerate(cube_list[1:]):
        try: cube.remove_coord('year')
        except: pass
        data = np.ma.masked_invalid(cube.data[0]).copy()
        data = np.ma.clip(data, 1., 1.)
        print('calc count:', cube_count.data.min(), cube_count.data.max(), data.min(),data.max())

        cube_count.data +=  data
        #cube_count.data+=(~np.ma.masked_invalid(cube.data).mask).astype(float)

    return cube_count

def make_count_of_cube_list_notime(cube_list,):
    """
    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.
    """
    # Fix empty times
    cube_count = cube_list[0].copy()

    #cube_count.data = (~np.ma.masked_invalid(cube_count.data)).mask.astype(float)
    cube_count.data = np.ma.masked_invalid(cube_count.data)
    cube_count.data = np.ma.clip(cube_count.data, 1.,1.)

    if len(cube_list)==1:
        return cube_count

    for i, cube in enumerate(cube_list[1:]):
        try: cube.remove_coord('year')
        except: pass
        data = np.ma.masked_invalid(cube.data).copy()
        data = np.ma.clip(data, 1., 1.)
        print('calc count:', cube_count.data.min(), cube_count.data.max(), data.min(),data.max())

        cube_count.data +=  data
        #cube_count.data+=(~np.ma.masked_invalid(cube.data).mask).astype(float)

    return cube_count


def make_map_plots(
        cfg,
        metadata,
        filename,
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

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        layer = str(layer)

        qplt.contourf(cube_layer, 25, linewidth=0, rasterized=True)

        try:
            plt.gca().coastlines()
        except AttributeError:
            logger.warning('Not able to add coastlines')

        # Add title to plot
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        if layer:
            title = ' '.join([
                title, '(', layer,
                str(cube_layer.coords('depth')[0].units), ')'
            ])
        plt.title(title)

        # Determine image filename:
        if multi_model:
            path = diagtools.folder(
                cfg['plot_dir']) + os.path.basename(filename).replace(
                    '.nc', '_map_' + str(layer_index) + image_extention)
        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='map_' + str(layer_index) + image_extention,
            )

        # Saving files:
        if cfg['write_plots']:

            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def make_map_contour(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple contour map plot for an individual model.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary
    filename: str
        the preprocessed model file.

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention and threshold.thresholds.
    image_extention = diagtools.get_image_format(cfg)

    # Load threshold/thresholds.
    plot_details = {}
    colours = []
    thresholds = diagtools.load_thresholds(cfg, metadata)

    for itr, thres in enumerate(thresholds):
        if len(thresholds) > 1:
            colour = plt.cm.jet(float(itr) / float(len(thresholds) - 1.))
        else:
            colour = plt.cm.jet(0)
        label = str(thres) + ' ' + str(cube.units)
        colours.append(colour)
        plot_details[thres] = {'c': colour,
                               'lw': 1,
                               'ls': '-',
                               'label': label}

    linewidths = [1 for thres in thresholds]
    linestyles = ['-' for thres in thresholds]
    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        layer = str(layer)
        qplt.contour(cube_layer,
                     thresholds,
                     colors=colours,
                     linewidths=linewidths,
                     linestyles=linestyles,
                     rasterized=True)

        try:
            plt.gca().coastlines()
        except AttributeError:
            logger.warning('Not able to add coastlines')
        try:
            plt.gca().add_feature(cartopy.feature.LAND,
                                  zorder=10,
                                  facecolor=[0.8, 0.8, 0.8])
        except AttributeError:
            logger.warning('Not able to add coastlines')
        # Add legend
        diagtools.add_legend_outside_right(plot_details,
                                           plt.gca(),
                                           column_width=0.02,
                                           loc='below')

        # Add title to plot
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        depth_units = str(cube_layer.coords('depth')[0].units)
        if layer:
            title = '{} ({} {})'.format(title, layer, depth_units)
        plt.title(title)

        # Determine image filename:
        if multi_model:
            path = os.path.join(diagtools.folder(cfg['plot_dir']),
                                os.path.basename(filename))
            path = path.replace('.nc', '_contour_map_' + str(layer_index))
            path = path + image_extention
        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='_contour_map_' + str(layer_index) + image_extention,
            )

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def multi_model_contours(
        cfg,
        metadata,
):
    """
    Make a contour map showing several models.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary.

    """
    ####
    # Load the data for each layer as a separate cube
    model_cubes = {}
    layers = {}
    for filename in sorted(metadata):
        cube = iris.load_cube(filename)
        cube = diagtools.bgc_units(cube, metadata[filename]['short_name'])

        cubes = diagtools.make_cube_layer_dict(cube)
        model_cubes[filename] = cubes
        for layer in cubes:
            layers[layer] = True

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Load threshold/thresholds.
    thresholds = diagtools.load_thresholds(cfg, metadata)

    # Make a plot for each layer and each threshold
    for layer, threshold in itertools.product(layers, thresholds):

        title = ''
        z_units = ''
        plot_details = {}
        cmap = plt.cm.get_cmap('jet')
        land_drawn = False

        # Plot each file in the group
        for index, filename in enumerate(sorted(metadata)):

            if len(metadata) > 1:
                color = cmap(index / (len(metadata) - 1.))
            else:
                color = 'blue'
            linewidth = 1.
            linestyle = '-'

            # Determine line style for Observations
            if metadata[filename]['project'] in diagtools.get_obs_projects():
                color = 'black'
                linewidth = 1.7
                linestyle = '-'

            # Determine line style for MultiModel statistics:
            if 'MultiModel' in metadata[filename]['dataset']:
                color = 'black'
                linestyle = ':'
                linewidth = 1.4

            cube = model_cubes[filename][layer]
            qplt.contour(cube,
                         [threshold, ],
                         colors=[color, ],
                         linewidths=linewidth,
                         linestyles=linestyle,
                         rasterized=True)
            plot_details[filename] = {
                'c': color,
                'ls': linestyle,
                'lw': linewidth,
                'label': metadata[filename]['dataset']
            }

            if not land_drawn:
                try:
                    plt.gca().coastlines()
                except AttributeError:
                    logger.warning('Not able to add coastlines')
                plt.gca().add_feature(cartopy.feature.LAND,
                                      zorder=10,
                                      facecolor=[0.8, 0.8, 0.8])
                land_drawn = True

            title = metadata[filename]['long_name']
            if layer != '':
                z_units = model_cubes[filename][layer].coords('depth')[0].units
            units = str(model_cubes[filename][layer].units)

        # Add title, threshold, legend to plots
        title = ' '.join([title, str(threshold), units])
        if layer:
            title = ' '.join([title, '(', str(layer), str(z_units), ')'])
        plt.title(title)
        plt.legend(loc='best')

        # Saving files:
        if cfg['write_plots']:
            path = diagtools.get_image_path(
                cfg,
                metadata[filename],
                prefix='MultipleModels_',
                suffix='_'.join(['_contour_map_',
                                 str(threshold),
                                 str(layer) + image_extention]),
                metadata_id_list=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
                ],
            )

        # Resize and add legend outside thew axes.
        plt.gcf().set_size_inches(9., 6.)
        diagtools.add_legend_outside_right(
            plot_details, plt.gca(), column_width=0.15)

        logger.info('Saving plots to %s', path)
        plt.savefig(path)
        plt.close()


def map_plot(cfg, metadata, cube, unique_keys = []):

    image_extention = diagtools.get_image_format(cfg)

    path = diagtools.folder([cfg['plot_dir'], 'single_map_plot']) + '_'.join(unique_keys)
    path+='_map'+image_extention

    if os.path.exists(path):return

    fig = plt.figure()

    # Making plots for each layer
    if 'counts' in unique_keys:
         qplt.contourf(cube, 11, vmin=0.,vmax=50., )
         plt.clim(0., 50.,)
         #plt.update()
    else:
         qplt.contourf(cube, 11,)

    plt.gca().coastlines()

    # Add title to plot
    title = ' '.join(unique_keys)
    plt.title(title)

    logger.info('Saving plots to %s', path)
    plt.savefig(path)

    plt.close()


# def make_file_mean(cfg, metadata, fn, operator='mean', ):
#
#     unique_keys = [metadata[k] for k in ['dataset', 'short_name',]]
#     unique_keys.append( operator)
#
#     work_dir = diagtools.folder([cfg['work_dir'], 'variable_group_means'])
#     path = work_dir+'_'.join(list(unique_keys))+'.nc'
#     print('make_file_mean, path:', path)
#     if os.path.exists(path):return path
#
#     cube = iris.load_cube( fn)
#     cube = diagtools.bgc_units(cube, metadata['short_name'])
#
#     cube = climate_statistics(cube, operator=operator)
#     print('Saving output:', path)
#     iris.save(cube, path)
#
#     map_plot(cfg, metadata, cube, unique_keys = unique_keys)
#
#     return path


def calculate_ensemble_stats(cfg, metadatas, key_short_name, key_scenario, key_start_year):
    """
    Calculate the mean, std_dev & number of contributing multi_model_statistics
    for each point in the 1x1 grid.
    """
    unique_keys = [key_short_name, key_scenario, str(key_start_year)]

    work_dir = diagtools.folder([cfg['work_dir'], 'variable_group_stats'])
    mean_timed_path = work_dir+'_'.join((unique_keys))+'_mean_timed.nc'
    mean_path = work_dir+'_'.join((unique_keys))+'_mean.nc'
    stddev_path = work_dir+'_'.join((unique_keys))+'_stddev.nc'
    count_path = work_dir+'_'.join((unique_keys))+'_count.nc'

    if False not in [os.path.exists(pth) for pth in [mean_timed_path, mean_path, stddev_path, count_path]]:
        make_three_pane_plot(cfg, unique_keys, mean_path, stddev_path, count_path)
        return

    files_list = []
    cubes = {}
    for filename, metadata in sorted(metadatas.items()):
        short_name = metadata['short_name']
        if short_name != key_short_name: continue

        scenario = metadata['exp']
        if scenario != key_scenario: continue

        start_year = metadata['start_year']
        if start_year != key_start_year: continue

        files_list.append(filename)
        cube = iris.load_cube(filename)

        # take mean along time axis.
        print('loading', unique_keys, metadata['dataset'], metadata['ensemble'])

        #cube = climate_statistics(cube, operator='mean')
        #cube = diagtools.bgc_units(cube, metadata['short_name'])
        cubes[filename] = cube

    if not len(files_list):
        return

    cube_list = [cube for fn, cube in cubes.items()]

    # mean_timed:
    # Timed included full time axis.
    if os.path.exists(mean_timed_path):
        mean_cube_timed = iris.load_cube(mean_timed_path)
    else:
        mean_cube_timed = make_mean_of_cube_list(cube_list)
        iris.save(mean_cube_timed, mean_timed_path)

    # mean:
    # Timed included full time axis. Mean is mean of that in the time direction.
    if os.path.exists(mean_path):
        mean_cube = iris.load_cube(mean_path)
    else:
        mean_cube = climate_statistics(mean_cube_timed, operator='mean')
        iris.save(mean_cube, mean_path)
    map_plot(cfg, {}, mean_cube, unique_keys =  [key_short_name, key_scenario, str(key_start_year), 'mean'])
    #mean_cube = diagtools.bgc_units(mean_cube, key_short_name)

    # Count cube
    if os.path.exists(count_path):
        count_cube = iris.load_cube(count_path)
    else:
        count_cube = make_count_of_cube_list(cube_list)
        iris.save(count_cube, count_path)
    map_plot(cfg, {}, count_cube, unique_keys =  [key_short_name, key_scenario, str(key_start_year), 'counts'])


    # std dev cube:
    if os.path.exists(stddev_path):
        pass
    else:
        stddev_cube = climate_statistics(mean_cube_timed, operator='std_dev')
        iris.save(stddev_cube, stddev_path)
    map_plot(cfg, {}, stddev_cube, unique_keys =  [key_short_name, key_scenario, str(key_start_year), 'stddev'])

    # make a plot:
    make_three_pane_plot(cfg, unique_keys, mean_path, stddev_path, count_path)




def calculate_ensemble_stats_wrong(cfg, metadatas, key_short_name, key_scenario, key_start_year):
    """
    Calculate the mean, std_dev & number of contributing multi_model_statistics
    for each point in the 1x1 grid.
    """
    assert 0
    # this method is wrong. it takes the time average first!

    unique_keys = [key_short_name, key_scenario, str(key_start_year)]

    work_dir = diagtools.folder([cfg['work_dir'], 'variable_group_stats'])
    mean_path = work_dir+'_'.join((unique_keys))+'_mean.nc'
    stddev_path = work_dir+'_'.join((unique_keys))+'_stddev.nc'
    count_path = work_dir+'_'.join((unique_keys))+'_count.nc'

    if False not in [os.path.exists(pth) for pth in [mean_path, stddev_path, count_path]]:
        make_three_pane_plot(cfg, unique_keys, mean_path, stddev_path, count_path)
        return

    files_list = []
    cubes = {}
    for filename, metadata in sorted(metadatas.items()):
        short_name = metadata['short_name']
        if short_name != key_short_name: continue

        scenario = metadata['exp']
        if scenario != key_scenario: continue

        start_year = metadata['start_year']
        if start_year != key_start_year: continue

        files_list.append(filename)
        cube = iris.load_cube(filename)

        # take mean along time axis.
        print('loading', unique_keys, metadata['dataset'], metadata['ensemble'])

        cube = climate_statistics(cube, operator='mean')
        #cube = diagtools.bgc_units(cube, metadata['short_name'])
        cubes[filename] = cube

    if not len(files_list):
        return

    cube_list = [cube for fn, cube in cubes.items()]

    # mean:
    if os.path.exists(mean_path):
        mean_cube = iris.load_cube(mean_path)
    else:
        mean_cube = make_mean_of_cube_list_notime(cube_list)
        iris.save(mean_cube, mean_path)
    map_plot(cfg, {}, mean_cube, unique_keys =  [key_short_name, key_scenario, str(key_start_year), 'mean'])
    #mean_cube = diagtools.bgc_units(mean_cube, key_short_name)

    # Count cube
    if os.path.exists(count_path):
        count_cube = iris.load_cube(count_path)
    else:
        count_cube = make_count_of_cube_list_notime(cube_list)
        iris.save(count_cube, count_path)
    map_plot(cfg, {}, count_cube, unique_keys =  [key_short_name, key_scenario, str(key_start_year), 'counts'])


    # std dev cube:
    if os.path.exists(stddev_path):
        pass
    else:
        stddev_cube = make_stddev_of_cube_list_notime(cube_list, mean_cube=mean_cube, count_cube=count_cube)
        iris.save(stddev_cube, stddev_path)
    map_plot(cfg, {}, stddev_cube, unique_keys =  [key_short_name, key_scenario, str(key_start_year), 'stddev'])

    # make a plot:
    make_three_pane_plot(cfg, unique_keys, mean_path, stddev_path, count_path)


def make_three_pane_plot(cfg, unique_keys, mean_path, stddev_path, count_path):
    """
    Make a three pane plot showing mean, std, count.
    """
    #make a plot
    image_extention = diagtools.get_image_format(cfg)
    path = diagtools.folder([cfg['plot_dir'],'three_pane_plot'])
    path += '_'.join(unique_keys)+image_extention

    if os.path.exists(path):return

    mean_cube= iris.load_cube(mean_path)
    stddev_cube = iris.load_cube(stddev_path)
    count_cube= iris.load_cube(count_path)

    fig = plt.figure()
    fig.set_size_inches(10., 8.)

    for sbp, cube, title in zip(
            [221,222,223],
            [mean_cube,stddev_cube, count_cube],
            ['Mean', 'Standard deviation', 'Count'],
            ):
        fig.add_subplot(sbp)
        if title == 'Count':
            qplt.contourf(cube, 11, vmin=0.,vmax=50.)
            plt.clim(0., 50.)
            #plt.update()
        else:
            qplt.contourf(cube, 11,)
        plt.title(title)
        plt.gca().coastlines()


    # Add title to plot
    title = ' '.join(unique_keys)
    plt.suptitle(title)
    logger.info('Saving plots to %s', path)
    plt.savefig(path)
    plt.close()


def main(cfg):
    """
    Load the config file, and send it to the plot makers.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    # cartopy.config['data_dir'] = cfg['auxiliary_data_dir']

    short_names = {}
    datasets = {}
    ensembles = {}
    start_years = {}
    scenarios = {}

    metadatas = diagtools.get_input_files(cfg, )
    for fn, metadata in metadatas.items():
        short_names[metadata['short_name']] = True
        datasets[metadata['dataset']] = True
        ensembles[metadata['ensemble']] = True
        start_years[metadata['start_year']] = True
        scenarios[metadata['exp']] = True

    for short_name, scenario, start_year in itertools.product(short_names, scenarios, start_years):
        if short_name == 'areacello': continue
        print('main', short_name, scenario, start_year)
        calculate_ensemble_stats(cfg, metadatas, short_name, scenario, start_year)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
