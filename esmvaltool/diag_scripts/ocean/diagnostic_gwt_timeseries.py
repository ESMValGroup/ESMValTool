"""
GWT Time series diagnostics
=======================

Diagnostic to produce figures of the time development of the GWT field from
cubes. These plost show time on the x-axis and cube value (ie temperature) on
the y-axis.

Two types of plots are produced: individual model timeseries plots and
multi model time series plots. The inidivual plots show the results from a
single cube, even if this is a mutli-model mean made by the _multimodel.py
preproccessor. The multi model time series plots show several models
on the same axes, where each model is represented by a different line colour.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has a time component, no depth component, and no
latitude or longitude coordinates.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_timeseries_1:# For Global Volume Averaged
      volume_statistics:
        operator: mean


An approproate preprocessor for a 3D+time field at the surface would be::

    prep_timeseries_2: # For Global surface Averaged
      extract_levels:
        levels:  [0., ]
        scheme: linear_extrap
      area_statistics:
        operator: mean


An approproate preprocessor for a 2D+time field would be::

    prep_timeseries_2: # For Global surface Averaged
      area_statistics:
        operator: mean


This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""

import logging
import os
import datetime
import iris
import matplotlib.pyplot as plt
import numpy as np
from itertools import product

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def timeplot(cube, **kwargs):
    """
    Create a time series plot from the cube.

    Note that this function simple does the plotting, it does not save the
    image or do any of the complex work. This function also takes and of the
    key word arguments accepted by the matplotlib.pyplot.plot function.
    These arguments are typically, color, linewidth, linestyle, etc...

    If there's only one datapoint in the cube, it is plotted as a
    horizontal line.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube

    """
    cubedata = np.ma.array(cube.data)
    if len(cubedata.compressed()) == 1:
        plt.axhline(cubedata.compressed(), **kwargs)
        return

    times = diagtools.cube_time_to_float(cube)
    plt.plot(times, cubedata, **kwargs)


def moving_average(cube, window):
    """
    Calculate a moving average.

    The window is a string which is a number and a measuremet of time.
    For instance, the following are acceptable window strings:

    * ``5 days``
    * ``12 years``
    * ``1 month``
    * ``5 yr``

    Also note the the value used is the total width of the window.
    For instance, if the window provided was '10 years', the the moving
    average returned would be the average of all values within 5 years
    of the central value.

    In the case of edge conditions, at the start an end of the data, they
    only include the average of the data available. Ie the first value
    in the moving average of a ``10 year`` window will only include the average
    of the five subsequent years.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube
    window: str
        A description of the window to use for the

    Returns
    ----------
    iris.cube.Cube:
        A cube with the movinage average set as the data points.

    """
    window = window.split()
    window_len = int(window[0]) / 2.
    win_units = str(window[1])

    if win_units not in [
            'years', 'yrs',
            'year', 'yr'
    ]:
        raise ValueError("Moving average window units not recognised: " +
                         "{}".format(win_units))

    times = cube.coord('time').units.num2date(cube.coord('time').points)

    cal_dt = diagtools.guess_calendar_datetime(cube)

    output = []

    times = np.array([
        cal_dt(time_itr.year, time_itr.month, time_itr.day, time_itr.hour,
               time_itr.minute) for time_itr in times
    ])

    for time_itr in times:
        tmin = cal_dt(time_itr.year - window_len, time_itr.month,
                      time_itr.day, time_itr.hour, time_itr.minute)

        tmax = cal_dt(time_itr.year + window_len, time_itr.month,
                      time_itr.day, time_itr.hour, time_itr.minute)

        #if time_itr.year - times.min().year < window_len:

        arr = np.ma.masked_where((times < tmin) + (times > tmax), cube.data)
        output.append(arr.mean())
    cube.data = np.array(output)
    return cube


def calculate_anomaly(cube, anomaly):
    """
    Calculate the anomaly using a specified time range.

    The anomaly window is a list which includes a starting year and and end
    year to indicate the start and end of the time period in which to calculate
    the anomaly.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube
    anomaly: list
        A start year and end year to calculate an anomaly.

    Returns
    ----------
    iris.cube.Cube:
        A cube with the anomaly calculated.
    """
    start_year = int(np.array(anomaly).min())
    end_year = int(np.array(anomaly).max())

    end_day = 31
    time_units = cube.coord('time').units
    if time_units.calendar == '360_day':
        end_day = 30

    start_date = datetime.datetime(int(start_year), 1, 1)
    end_date = datetime.datetime(int(end_year), 12, end_day)

    t_1 = time_units.date2num(start_date)
    t_2 = time_units.date2num(end_date)
    constraint = iris.Constraint(
        time=lambda t: t_1 < time_units.date2num(t.point) < t_2)

    new_cube = cube.extract(constraint)
    if new_cube is None:
        return None
    mean = new_cube.data.mean()
    cube.data = cube.data - mean
    return cube


def get_threshold_exceedance_date(cube, threshold):
    """
    Calculate the threshold exceedance date.
    """
    loc = np.where(cube.data > threshold)[0]
    print('exceedance indices:', loc)
    if not len(loc): return None
    times = cube.coord('time').units.num2date(
        cube.coord('time').points)
    time = times[loc[0]]
    return time


def make_time_series_plots(
        cfg,
        metadata,
        filename,
        thresholds = [1.5, 2., 3., 4., 5., 6.]
):
    """
    Make a simple time series plot for an indivudual model 1D cube.

    This tool loads the cube from the file, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.
    Returns:
    --------
        dict

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

    exceedance_dates = {}
    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        layer = str(layer)
        print(metadata['ensemble'])
        if metadata['ensemble'] in ['r3i1p1f2', 'r4i1p1f2','r5i1p1f2', 'r8i1p1f2',]: assert 0

        if 'anomaly' in cfg:
            cube_layer = calculate_anomaly(cube_layer, cfg['anomaly'])
            if cube_layer is None:
                return

        if 'moving_average' in cfg:
            cube_layer = moving_average(cube_layer, cfg['moving_average'])

        if multi_model:
            timeplot(cube_layer, label=metadata['dataset'], ls=':')
        else:
            timeplot(cube_layer, label=metadata['dataset'])

        # Add title, legend to plots
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        if 'anomaly' in cfg:
            title = ' '.join([title, 'anomaly'])
        if layer != '':
            if cube_layer.coords('depth'):
                z_units = cube_layer.coord('depth').units
            else:
                z_units = ''
            title = ' '.join([title, '(', layer, str(z_units), ')'])
        plt.title(title)
        #plt.legend(loc='best')

        ylabel = str(cube_layer.units)
        if 'ylabel' in cfg:
            ylabel = cfg['ylabel']
        plt.ylabel(ylabel)

        # Determine image filename:
        if multi_model:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                prefix='MultiModel',
                suffix='_'.join(['timeseries',
                                 str(layer) + image_extention]),
                metadata_id_list=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
                ],
            )

        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='timeseries_' + str(layer_index) + image_extention,
            )

        # Add thresholds part
        for threshold in thresholds:
            if threshold > np.max(plt.ylim()): continue
            plt.axhline(threshold, color='k', ls='--', lw=0.5)
        plt.axhline(0., color='k', ls='--', lw=0.5)

        # Add thresholds legend:
        txt = [metadata['dataset'], metadata['exp'], metadata['ensemble']]

        print('----------\nSearching for thresholds', txt)
        for threshold in thresholds:
            time = get_threshold_exceedance_date(cube, threshold)
            print('threshold:', threshold, time)
            if not time:
                continue
            txt.append('>'+str(threshold)+': '+str(time.year))
            plt.axvline(time.year, color='k', ls='--', lw=0.5,)

            exd_key = (metadata['exp'], metadata['ensemble'], threshold)
            print(metadata['ensemble'])
            if metadata['ensemble'] in ['r3i1p1f2', 'r4i1p1f2','r5i1p1f2', 'r8i1p1f2',]: assert 0
            exceedance_dates[exd_key] = time.year

        # Add top left legend
        ntxt = plt.text(.05, .95, '\n'.join(txt),
            transform=plt.gca().transAxes, ha="left", va="top")
        ntxt.set_bbox(dict(facecolor='grey', alpha=0.7, edgecolor='black'))

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()
    return exceedance_dates


def multi_model_time_series(
        cfg,
        metadata,
        thresholds = [1.5, 2., 3., 4., 5., 6.],
):
    """
    Make a time series plot showing several preprocesssed datasets.

    This tool loads several cubes from the files, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.

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

    # Make a plot for each layer
    for layer in layers:
        title = ''
        z_units = ''
        plot_details = {}
        cmap = plt.cm.get_cmap('viridis')

        # Plot each file in the group
        for index, filename in enumerate(sorted(metadata)):
            if len(metadata) > 1:
                color = cmap(index / (len(metadata) - 1.))
            else:
                color = 'blue'
            cube = model_cubes[filename][layer]
            if 'anomaly' in cfg:
                cube = calculate_anomaly(cube, cfg['anomaly'])
                if cube is None:
                    logger.warning('Not enough time to calculate anomaly: %s',
                                   metadata[filename]['dataset'])
                    continue

            # Take a moving average, if needed.
            if 'moving_average' in cfg:
                cube = moving_average(cube, cfg['moving_average'])

            if metadata[filename]['dataset'].lower().find('multimodel') > -1:
                logger.info('plotting - Multi: %s',
                            metadata[filename]['dataset'])
                timeplot(
                    cube,
                    c='black',
                    ls='--',
                    lw=2.,
                )
                plot_details[filename] = {
                    'c': 'black',
                    'ls': '--',
                    'lw': 2.,
                    'label': metadata[filename]['dataset']
                }
            else:
                timeplot(
                    cube,
                    c=color,
                    ls='-',
                    lw=2.,
                )
                plot_details[filename] = {
                    'c': color,
                    'ls': '-',
                    'lw': 2.,
                    'label': metadata[filename]['dataset']
                }

            title = metadata[filename]['long_name']
            ylabel = str(model_cubes[filename][layer].units)
            if layer != '':
                if model_cubes[filename][layer].coords('depth'):
                    z_units = model_cubes[filename][layer].coord('depth').units
                else:
                    z_units = ''

        # Add thresholds part
        for threshold in thresholds:
            if threshold > np.max(plt.ylim()):
                continue
            plt.axhline(threshold, color='k', ls='--', lw=0.5)
        plt.axhline(0., color='k', ls='--', lw=0.5)

        # Add title, legend to plots
        if 'anomaly' in cfg:
            title = ' '.join([title, 'anomaly'])

        if layer:
            title = ' '.join([title, '(', str(layer), str(z_units), ')'])

        # check to see if the title is mentionned in the recipe.
        # If so, it overwrites the dafult title.
        if 'title' in cfg:
            title = cfg['title']

        if 'ylabel' in cfg:
            ylabel = cfg['ylabel']

        plt.title(title)
        plt.legend(loc='best')
        plt.ylabel(ylabel)

        # Saving files:
        if cfg['write_plots']:
            path = diagtools.get_image_path(
                cfg,
                metadata[filename],
                prefix='MultipleModels',
                suffix='_'.join(['timeseries',
                                 str(layer) + image_extention]),
                metadata_id_list=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
                ],
            )

        # Resize and add legend outside thew axes.
        if len(plot_details) < 25:
            plt.gcf().set_size_inches(9., 6.)
            diagtools.add_legend_outside_right(
                plot_details, plt.gca(), column_width=0.15)
        if len(plot_details) > 25:
            plt.gcf().set_size_inches(11., 6.)
            diagtools.add_legend_outside_right(
                plot_details, plt.gca(), column_width=0.18)

        logger.info('Saving plots to %s', path)
        plt.savefig(path)
        plt.close()


def print_exceedance_dates(cfg, exceedance_dates, window = 10, short_name = 'tas',mip='Amon', preprocessor='prep_1', grid = 'gn'):
    """
    prints the exceedance_dates in a format ready to go into a recipe.
    exceednace key: (metadata['exp'], metadata['ensemble'], threshold)
    Like this:
      tas_ssp119_15:
        short_name: tas
        preprocessor: prep_1
        additional_datasets:
          - {dataset: UKESM1-0-LL, project: CMIP6, mip: Amon, exp: [historical, ssp119], ensemble: r1i1p1f2, start_year: 2014, end_year: 2034, grid: gn}
          - {dataset: UKESM1-0-LL, project: CMIP6, mip: Amon, exp: [historical, ssp119], ensemble: r3i1p1f2, start_year: 2013, end_year: 2033, grid: gn}
          - {dataset: UKESM1-0-LL, project: CMIP6, mip: Amon, exp: [historical, ssp119], ensemble: r4i1p1f2, start_year: 2017, end_year: 2037, grid: gn}
    """
    # Define the contents
    exps = set()
    ensembles = set()
    thresholds = set()

    for (exp, ens, thresh) in exceedance_dates.keys():
        if exp == 'historical':
            continue
        # print(exp, exp.find('-')+1, exp[exp.find('-')+1:])
        # exp = exp[exp.find('-'):]
        exps.add(exp)
        ensembles.add(ens)
        thresholds.add(thresh)
    print(exceedance_dates)

    txt=''
    # For each combination of short_name, threshold:
    for exp, thresh in product(sorted(exps), sorted(thresholds)):
        ssp = exp[exp.find('-')+1:]
        lines = []
        lines.append('\n') #
        lines.append('      '+ '_'.join([short_name, ssp, str(thresh)])) #  tas_ssp119_15:
        lines.append('        short_name: '+ short_name)
        lines.append('        preprocessor: '+ preprocessor)
        lines.append('        additional_datasets:')

        # matches = []
        for ens in ensembles:
            print(exp, thresh, ens)
            try:
                exceedance_date = float(exceedance_dates[(exp, ens, thresh)])
            except:
                continue

            start_year = str(exceedance_date - window)
            end_year = str(exceedance_date + window)

            lines.append('         - {'
                         'dataset: UKESM1-0-LL, '
                         'project: CMIP6, '
                         'mip: '+mip+', '
                         'exp: [historical, '+ssp+'], '
                         'ensemble: '+ens +', '
                         'start_year: '+start_year+', '
                         'end_year: '+end_year+', '
                         'grid:'+grid+', }'
                         )
        if len(lines) == 5:
            continue
        txt += '\n'.join(lines)

    txt += '\n'
    print(txt)
    fn = cfg['work_dir']+'/new_recipe.yml'
    print('Saved to: ', fn)
    out = open(fn, 'w')
    out.write(txt)
    out.close()


def main(cfg):
    """
    Load the config file and some metadata, then make plots.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    exceedance_dates = {}
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('metadata filename:\t%s', metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)

        #######
        # Multi model time series
        multi_model_time_series(
            cfg,
            metadatas,
        )

        for filename in sorted(metadatas):

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            ######
            # Time series of individual model
            print(metadatas[filename]['ensemble'])
            if metadatas[filename]['ensemble'] in ['r3i1p1f2', 'r4i1p1f2','r5i1p1f2', 'r8i1p1f2',]: assert 0
            exceedance_date = make_time_series_plots(cfg, metadatas[filename], filename)
            exceedance_dates.update(exceedance_date)
    print(exceedance_dates)
    print_exceedance_dates(cfg, exceedance_dates)
    logger.info('Success')



if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
