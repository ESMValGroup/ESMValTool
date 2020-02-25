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
import cf_units
import glob

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
    if 'time' not in [coord.name() for coord in cube.coords()]:
        return

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


def calculate_anomaly(cube, anomaly, calc_average=False):
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
    calc_average: bool
        Flag to return the average of the anomaly period, instead of subtracting
        it from the cube.
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
    if calc_average:
       return mean
    cube.data = cube.data - mean
    return cube


def get_threshold_exceedance_date(cube, threshold):
    """
    Calculate the threshold exceedance date.
    """
    loc = np.where(cube.data > threshold)[0]
    if not len(loc): return None
    times = cube.coord('time').units.num2date(
        cube.coord('time').points)
    time = times[loc[0]]
    return time


def calculate_total(cfg, metadata, cube):
    """
    Calcualte the global total in the cube
    """
    if metadata['short_name'] in ['npp', 'rh']:
        area = 'areacella'
    elif metadata['short_name'] in ['npp', 'rh']:
        area = 'areacella'
    else:
        print (metadata['short_name'], 'not recognised')
        assert 0

    return



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
    #cube = diagtools.bgc_units(cube, metadata['short_name'])

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
        cube = calculate_total(cfg, metadata, cube)
        #cube = diagtools.bgc_units(cube, metadata[filename]['short_name'])

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

    txt='      # Procedurally generated recipe contents:'
    # Add the historical ensembles:
    lines = []
    lines.append('\n') #
    lines.append('      '+ '_'.join([short_name, 'historical'])+':') #  tas_ssp119_15:
    lines.append('        short_name: '+ short_name)
    lines.append('        preprocessor: '+ preprocessor)
    lines.append('        additional_datasets:')
    for ens in sorted(ensembles):
        lines.append('         - {'
                     'dataset: UKESM1-0-LL, '
                     'project: CMIP6, '
                     'mip: ' + mip + ', '
                     'exp: historical, '
                     'ensemble: ' + ens + ', '
                     'start_year: 1850, '
                     'end_year: 1900, '
                     'grid: ' + grid + '}'
                     )
    txt += '\n'.join(lines)

    # For each combination of short_name, threshold:
    for exp, thresh in product(sorted(exps), sorted(thresholds)):
        ssp = exp[exp.find('-')+1:]
        lines = []
        lines.append('\n') #
        lines.append('      '+ '_'.join([short_name, ssp, str(thresh)])+':') #  tas_ssp119_15:
        lines.append('        short_name: '+ short_name)
        lines.append('        preprocessor: '+ preprocessor)
        lines.append('        additional_datasets:')

        # matches = []
        for ens in sorted(ensembles):
            print(exp, thresh, ens)
            try:
                exceedance_date = float(exceedance_dates[(exp, ens, thresh)])
            except:
                continue

            start_year = str(int(exceedance_date - window))
            end_year = str(int(exceedance_date + window))

            # What if end year is after 2100?
            if int(end_year)> 2099:
                continue

            lines.append('         - {'
                         'dataset: UKESM1-0-LL, '
                         'project: CMIP6, '
                         'mip: ' + mip + ', '
                         'exp: [historical, ' + ssp + '], '
                         'ensemble: ' + ens + ', '
                         'start_year: ' + start_year + ', '
                         'end_year: ' + end_year + ', '
                         'grid: ' + grid + '}'
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

def marine_gt(data_dict, short, gt):
    """
    Calculate global from the data dictionary.
    """
    areas = []
    for (short_name, exp, ensemble), cube in sorted(data_dict.items()):
        if short_name == 'areacello':
            areas.append(cube)
    if len(areas) != 1:
        assert 0
    areas = areas[0]
    for (short_name, exp, ensemble), cube in sorted(data_dict.items()):
        if short_name != short:
            continue
        cubegt = cube.copy()
        cubegt.data = cube.data * areas.data * 1.E-12 * (360*24*60*60)
        cubegt.units = cf_units.Unit('Pg yr^-1') #cube.units * areas.units

        data_dict[(gt, exp, ensemble)] = cubegt
    return data_dict


def fgco2gt(data_dict): return marine_gt(data_dict, short='fgco2', gt='fgco2gt')
def intppgt(data_dict): return marine_gt(data_dict, short='intpp', gt='intppgt')
def epc100gt(data_dict): return marine_gt(data_dict, short='epc100', gt='epc100gt')
def intdicgt(data_dict): return marine_gt(data_dict, short='intdic', gt='intdicgt')
def intpocgt(data_dict): return marine_gt(data_dict, short='intpoc', gt='intpocgt')
def fricgt(data_dict): return marine_gt(data_dict, short='fric', gt='fricgt')
def frocgt(data_dict): return marine_gt(data_dict, short='froc', gt='frocgt')


def nppgt(data_dict, short='npp', gt='nppgt'):
    """
    Calculate nppgt from the data dictionary.
    """
    areas = []
    for (short_name, exp, ensemble), cube in sorted(data_dict.items()):
        if short_name == 'areacella':
            areas.append(cube)
    if len(areas) != 1:
        print(areas)
        assert 0
    areas = areas[0]
    for (short_name, exp, ensemble), cube in sorted(data_dict.items()):
        #rint('nppgt:', short_name, exp, ensemble, [short,gt])
        if short_name != short:
        #   print('nppgt:', short_name,'!=', short)
            continue
        if (gt, exp, ensemble) in data_dict.keys():
            print('nppgt:', gt, 'already calculated')
            continue
        print('nppgt:', short_name, exp, ensemble, [short,gt])
        cubegt = cube.copy()
        cubegt.data = cube.data * areas.data * 1.E-12 * (360*24*60*60)
        cubegt.units = cf_units.Unit('Pg yr^-1') #cube.units * areas.units
        print('nppgt:', (gt, exp, ensemble), cubegt.data.mean())
        data_dict[(gt, exp, ensemble)] = cubegt
    return data_dict


def rhgt(data_dict):
    """
    Calculate rhgt from the data dictionary.
    """
    return nppgt(data_dict, short='rh', gt='rhgt')


def exchange(data_dict):
    """
    Calculate exchange from the data dictionary.
    """
    data_dict = rhgt(data_dict)
    data_dict = nppgt(data_dict)

    exps = {}
    ensembles = {}
    for (short_name, exp, ensemble)  in sorted(data_dict.keys()):
        exps[exp] = True
        ensembles[ensemble] = True

    for exp, ensemble in product(exps, ensembles):
        if ('nppgt', exp, ensemble) not in data_dict: continue
        cube = data_dict[('nppgt', exp, ensemble)].copy()
        cube2 = data_dict[('rhgt', exp, ensemble)]
        cube.data = cube.data - cube2.data
        data_dict[('exchange', exp, ensemble)] = cube
    return data_dict


def tas_norm(data_dict):
    """
    Calculate tas_norm from the data dictionary.
    """
    exps = {}
    ensembles = {}
    baselines={}
    for (short_name, exp, ensemble), cube  in data_dict.items():
        exps[exp] = True
        ensembles[ensemble] = True
        if short_name != 'tas': continue
        if exp != 'historical': continue
        baselines[(short_name, ensemble)] = calculate_anomaly(cube, [1850, 1900], calc_average=True)

    for exp, ensemble in product(exps, ensembles):
        if not ('tas', exp, ensemble) in data_dict.keys(): continue
        cube = data_dict[('tas', exp, ensemble)].copy()
        cube.data = cube.data - baselines[('tas', ensemble)]
        data_dict[('tas_norm', exp, ensemble)] = cube
    return data_dict


def norm_co2(data_dict, short='nppgt'):
    """
    Weight a value according to the ratio of the forcing co2 for each year
    against the average co2 forcing in 1850-1900.
    """
    print(data_dict.keys())
    print(data_dict[('co2', 'historical', 'r1i1p1f2' )]['time'][:50],
          data_dict[('co2', 'historical', 'r1i1p1f2' )]['co2'][:50])
    baseline = np.mean(data_dict[('co2', 'historical', 'r1i1p1f2' )]['co2'][:50])
    new_data_dict = {}
    for (short_name, exp, ensemble), cube  in data_dict.items():
        if short_name != short: continue
        #if exp != 'historical': continue
        cube = data_dict[(short, exp, ensemble)].copy()
        print('norm_co2:', short_name, exp, ensemble, 'baseline:',baseline)
        out = []
        co2_data= data_dict[('co2', exp, ensemble )]['co2']
        if len(cube.data) != len(co2_data):
            times = cube.coord('time').units.num2date(cube.coord('time').points)
            print('times do not match', (short_name, exp, ensemble), len(cube.data), '!=', len(co2_data))
            for t1, t2 in zip(times, data_dict[('co2', exp, ensemble )]['time']):
                print(short_name, exp, ensemble, short_name+':', t1, 'co2:', t2)
            assert 0
        for d,co2 in zip(cube.data, data_dict[('co2', exp, ensemble)]['co2']):
            out.append(d*baseline/co2)
        cube.data = np.ma.array(out)
        new_data_dict[(short+'_norm', exp, ensemble)] = cube
    data_dict.update(new_data_dict)
    return data_dict


def norm_co2_nppgt(data_dict): return norm_co2(data_dict, short='nppgt')
def norm_co2_rhgt(data_dict): return norm_co2(data_dict, short='rhgt')
def norm_co2_exchange(data_dict): return norm_co2(data_dict, short='exchange')
def norm_co2_fgco2gt(data_dict): return norm_co2(data_dict, short='fgco2gt')


def load_timeseries(cfg, short_names):
    """
    Load times series as a dict.

    Dict is :
    data_dict[(short_name, exp, ensemble) ] = cube
    assume only one model
    """
    transforms = {
        'fgco2gt': ['fgco2', 'areacello'],
        'nppgt': ['npp', 'areacella'],
        'rhgt': ['rh', 'areacella'],
        'epc100gt': ['epc100', 'areacello'],
        'intdicgt': ['intdic', 'areacello'],
        'intpocgt': ['intpoc', 'areacello'],
        'fricgt': ['fric', 'areacello'],
        'frocgt': ['froc', 'areacello'],
        'exchange': ['rh', 'npp', 'areacella'],
        'tas_norm': ['tas', ],
        'nppgt_norm': ['nppgt', ],
        'rhgt_norm': ['rhgt', ],
        'exchange_norm': ['exchange', ],
        'fgco2gt_norm': ['fgco2gt', ],
        }

    transforms_functions = {
        'fgco2gt': fgco2gt,
        'nppgt': nppgt,
        'rhgt': rhgt,
        'epc100gt': epc100gt,
        'intdicgt': intdicgt,
        'intpocgt': intpocgt,
        'fricgt': fricgt,
        'frocgt': frocgt,
        'exchange': exchange,
        'tas_norm': tas_norm,
        'nppgt_norm':norm_co2_nppgt,
        'rhgt_norm':norm_co2_rhgt,
        'exchange_norm':norm_co2_exchange,
        'fgco2gt_norm':norm_co2_fgco2gt,
        }

    short_names_to_load = short_names.copy()
    for sn in short_names:
        if sn in transforms:
            short_names_to_load.extend(transforms[sn])

    data_dict = {}
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('load_timeseries:\t%s', metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)
        for fn in sorted(metadatas):
            short_name = metadatas[fn]['short_name']
            exp = metadatas[fn]['exp']
            ensemble = metadatas[fn]['ensemble']

            if short_name not in short_names_to_load:
                continue

            cube = iris.load_cube(fn)
            #cube = diagtools.bgc_units(cube, short_name)

            print('load_timeseries:\t%s successfull loaded data:', (short_name, exp, ensemble), 'mean:', cube.data.mean())

            data_dict[(short_name, exp, ensemble)] = cube

    if 'co2' in short_names_to_load:
        data_dict = load_co2_forcing(cfg, data_dict)


    for sn in short_names_to_load:
        if sn in transforms:
            data_dict = transforms_functions[sn](data_dict)
    return data_dict


def load_thresholds(cfg, data_dict, short_names = ['tas', ], thresholds = [1.5, 2., 3., 4., 5.], ):
    """
    Load thresholds  as a dict.

    Dict is :
    data_dict[(short_name, exp, ensemble) ] = {threshold: year}
    """
    thresholds_dict = {}
    baselines = {}
    for (short_name, exp, ensemble), cube in sorted(data_dict.items()):
        if short_name not in short_names:
             continue
        if exp != 'historical': continue
        baselines[(short_name, ensemble)] = calculate_anomaly(cube, [1850, 1900], calc_average=True)

    for (short_name, exp, ensemble), cube in sorted(data_dict.items()):
        if short_name not in short_names:
             continue

        cube2 = moving_average(cube.copy(), '21 years')
        cube2.data = cube2.data - baselines[(short_name, ensemble)]

        thresholds_dict[(short_name, exp, ensemble)] = {}

        for threshold in thresholds:
            time = get_threshold_exceedance_date(cube2, threshold)
            thresholds_dict[(short_name, exp, ensemble)][threshold] = time
    return thresholds_dict

def cube_to_years(cube):
    """
    Convert from time coordinate into years.

    Takes an iris time coordinate and returns a list of floats.
    Parameters
    ----------
    cube: iris.cube.Cube
        the opened dataset as a cube.

    Returns
    -------
    list
        List of floats showing the time coordinate in decimal time.

    """
    if not cube.coords('year'):
        iris.coord_categorisation.add_year(cube, 'time')
    return cube.coord('year').points


def load_co2_forcing(cfg, data_dict):
    """
    Load annual CO2 data from the auxiliary datasets.

    Unlike the rest of data_dcit, it's isn't loaded as a cube, but rather as a
    dict.
    """
    fold = cfg['auxiliary_data_dir']+'/atmos_co2_forcing/'
    files = glob.glob(fold+'*.dat')
    print(files)
    hist_datas = []
    hist_times = []
    for fn in files:
        open_fn = open(fn, 'r')
        key = os.path.basename(fn).replace('_co2.dat', '')
        times = []
        data = []
        for line in open_fn.readlines()[1:]:
            line = line.split(' ')
            for x in range(len(line)):
                if '' in line: line.remove('')
                if '\n' in line: line.remove('\n')
            times.append(float(line[0]))
            data.append(float(line[1]))
        if key == 'historical':
            hist_datas = np.array(data).copy()
            hist_times = np.array(times).copy()
        for ens in ['r1', 'r2',  'r3', 'r4', 'r8']:
            data_dict[('co2', key, ens+'i1p1f2' )] = {'time': times, 'co2':data}
            print('load_co2_forcing:\t%s successfull loaded data:', ('co2', key, ens+'i1p1f2'), 'mean:', np.array(data).mean())
        open_fn.close()

    # Check for historical-ssp scenarios.
    tmp_dict = {}
    for (short_name, exp, ensemble), ssp_cube in data_dict.items():
        if short_name in ['co2', 'areacella', 'areacello',]:
            continue
        if ('co2', 'historical-'+exp, ensemble ) in tmp_dict.keys():
            continue
        if exp == 'historical':
            continue
        ssp_only = exp.replace('historical-', '')
        new_times = []
        new_datas = []
        print((short_name, exp,(ssp_only), ensemble))
        ssp_times = cube_to_years(ssp_cube)
        min_time = np.array(ssp_times).min()
        if min_time > np.array(hist_times).max():
            print(short_name, exp, ensemble, 'no overlap', ('ssp:', min_time, '>', 'hist max:', np.array(hist_times).max()))
            # no overlap
            new_times = data_dict[('co2', ssp_only, ensemble)]['time']
            new_datas = data_dict[('co2', ssp_only, ensemble)]['co2']
        else:
            # Some overlap
            print(short_name, exp, ensemble,'some overlap', (min_time, '<=', np.array(hist_times).max()))
            new_times = list(np.ma.masked_where(hist_times<min_time, hist_times).compressed())
            new_datas = list(np.ma.masked_where(hist_times<min_time, hist_datas).compressed())
            new_times.extend(data_dict[('co2', ssp_only, ensemble)]['time'])
            new_datas.extend(data_dict[('co2', ssp_only, ensemble)]['co2'])

        print('co2', exp, ensemble, len(new_times), len(new_datas))
        tmp_dict[('co2', exp, ensemble )] ={'time': new_times, 'co2':new_datas}
    data_dict.update(tmp_dict)

    # Save the co2 image:
    path = diagtools.folder(cfg['plot_dir'])
    image_extention = diagtools.get_image_format(cfg)
    path += 'co2_forcing' + image_extention
    if not os.path.exists(path):
        exp_colours = {'historical':'black',
                       'ssp119':'green',
                       'ssp126':'dodgerblue',
                       'ssp245':'blue',
                       'ssp370':'purple',
                       'ssp434':'goldenrod',
                       'ssp585': 'red',
                       'ssp534-over':'orange'}
        for key in exp_colours.keys():
            plt.plot(data_dict[('co2', key, 'r1i1p1f2' )]['time'],
                     data_dict[('co2', key, 'r1i1p1f2' )]['co2'],
                     c=exp_colours[key],
                     label=key)
        plt.legend()
        plt.savefig(path)
        plt.close()
    image_extention = diagtools.get_image_format(cfg)
    path += 'co2_forcing_hists' + image_extention
    if not os.path.exists(path):
        exp_colours = {'historical':'black',
                       'historical-ssp119':'green',
                       'historical-ssp126':'dodgerblue',
                       'historical-ssp245':'blue',
                       'historical-ssp370':'purple',
                       'historical-ssp434':'goldenrod',
                       'historical-ssp585': 'red',
                       'historical-ssp534-over':'orange'}
        for key in exp_colours.keys():
            plt.plot(data_dict[('co2', key, 'r1i1p1f2' )]['time'],
                     data_dict[('co2', key, 'r1i1p1f2' )]['co2'],
                     c=exp_colours[key],
                     label=key)
        plt.legend()
        plt.savefig(path)
        plt.close()
    return data_dict


# def load_areas(cfg, short_names=['areacella', 'areacello']):
#     """
#     Load area data
#     """
#     area_dict = {}
#     for index, metadata_filename in enumerate(cfg['input_files']):
#         logger.info('load_areas:\t%s', metadata_filename)
#         metadatas = diagtools.get_input_files(cfg, index=index)
#         for fn in sorted(metadatas):
#             short_name = metadatas[fn]['short_name']
#             exp = metadatas[fn]['exp']
#             ensemble = metadatas[fn]['ensemble']
#             variable_group = metadatas[fn]['variable_group']
#
#             if short_name not in short_names: continue
#             print('loaded area:', (short_name, exp, ensemble) )
#             cube = iris.load_cube(fn)
#             area_dict[(short_name,variable_group)] = cube
#     return area_dict


def get_threshold_point(cube, year):
    """
    get the location of the year provided.
    """
    time_units = cube.coord('time').units
    date = datetime.datetime(int(year), 6, 1)
    t_1 = time_units.date2num(date)
    arg_min = np.argmin( np.abs(cube.coord('time').points - t_1))
    return arg_min


def make_ts_figure(cfg, data_dict, thresholds_dict, x='time', y='npp',markers='thresholds', draw_line=True, do_moving_average=True):
    """
    make a 2D figure.
    x axis and y axis are determined by the short_names provuided in x and y
    vars.
    Markers are placed at certain points when the tas goes above thresholds.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    """
    exps = {}
    ensembles = {}
    for (short_name, exp, ensemble)  in sorted(data_dict.keys()):
         exps[exp] = True
         ensembles[ensemble] = True

    exp_colours = {'historical':'black',
                   'ssp119':'green',
                   'ssp126':'dodgerblue',
                   'ssp245':'blue',
                   'ssp370':'purple',
                   'ssp434':'goldenrod',
                   'ssp585': 'red',
                   'ssp534-over':'orange',
                   'historical-ssp119':'green',
                   'historical-ssp126':'dodgerblue',
                   'historical-ssp245':'blue',
                   'historical-ssp370':'purple',
                   'historical-ssp434':'goldenrod',
                   'historical-ssp585': 'red',
                   'historical-ssp534-over':'orange'}

    marker_styles = {1.5: 'o', 2.:'*', 3.:'^', 4.:'s', 5.:'X'}

    fig = plt.figure()
    x_label,y_label = [], []
    for exp_1, ensemble_1 in product(exps, ensembles):
        x_data, y_data = [], []
        for (short_name, exp, ensemble), cube in sorted(data_dict.items()):
            #print('plotting', short_name, exp, ensemble, ' from', x,y)
            if exp != exp_1: continue
            if ensemble != ensemble_1: continue
            if short_name not in [x,y]: continue

            print('make_ts_figure: found', short_name, exp, ensemble, x,y)
            if x == 'time':
                x_label = 'Year'
            elif x == short_name == 'co2':
                x_data = cube['co2']
                x_label = ' '.join(['atmopheric co2'])
            elif short_name == x:
                x_data = cube.data
                x_label = ' '.join([x, str(cube.units)])

            if y == 'time':
                if isinstance(cube, dict):
                    x_data = cube['time']
                else:
                    y_data = diagtools.cube_time_to_float(cube)
                y_label = 'Year'
            elif y == short_name == 'co2' :
                y_data = cube['co2']
                y_label = ' '.join(['atmopheric co2'])
            elif short_name == y:
                y_data = cube.data
                y_label = ' '.join([y, str(cube.units)])

            if x == 'time' and short_name == y:
                x_label = 'Year'
                if isinstance(cube, iris.cube.Cube):
                    x_data = diagtools.cube_time_to_float(cube)
                else:
                    x_data = cube['time']
            print('make_ts_figure: loaded x data', short_name, exp, ensemble, x, np.mean(x_data))
            print('make_ts_figure: loaded y data', short_name, exp, ensemble, y, np.mean(y_data))



        if 0 in [len(x_data), len(y_data)]: continue

        label = ' '.join([exp_1, ensemble_1])
        if draw_line:
            plt.plot(x_data, y_data, lw=0.5, color=exp_colours[exp_1], )#label=label)

        if markers == 'thresholds':
            try: threshold_times = thresholds_dict[('tas', exp_1, ensemble_1)]
            except:
               threshold_times = {}
            for threshold, time in threshold_times.items():
                if not time: continue

                y_point = get_threshold_point(cube, time.year)
                plt.plot(x_data[y_point],
                         y_data[y_point],
                         marker_styles[threshold],
                         fillstyle='none',
                         color=exp_colours[exp_1])

    exp_colours_leg = {'historical':'black',
                   'ssp119':'green',
                   'ssp126':'dodgerblue',
                   'ssp245':'blue',
                   'ssp370':'purple',
                   'ssp434':'goldenrod',
                   'ssp585': 'red',
                   'ssp534-over':'orange'}
    plot_details = {}
    for exp,color in sorted(exp_colours_leg.items()):
        plot_details[exp] = {
                    'c': color,
                    'ls': '-',
                    'lw': 2.,
                    'label': exp
                }
    for thres,ms in sorted(marker_styles.items()):
        plot_details[str(thres)] = {
                    'c': 'black',
                    'marker': ms,
                    'fillstyle':'none',
                    'label': '>' + str(thres)+u'\u00B0C'
                }

    diagtools.add_legend_outside_right(
                plot_details, plt.gca(), column_width=0.175)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(' '.join([x, 'by', y ]))

    image_extention = diagtools.get_image_format(cfg)

    path = diagtools.folder(cfg['plot_dir'])

    path += '_'.join([x,y,markers,]) + image_extention
    if do_moving_average:
        path = path.replace(image_extention, '_21ma'+image_extention)
    print('saving figure:', path)
    plt.savefig(path)
    plt.close()


def main(cfg):
    """
    Load the config file and some metadata, then make plots.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    #jobss for tomoorrow:
    #    check to make sure that you're using the 'right areacella for land.
    #    do you even need the areacella for air? probably not, right?
    #    change the recipe to add the other ensemble members to the job.
    #    email the figues to other authors.



    # short_names = ['tas', 'tas_norm', 'nppgt', 'fgco2gt', 'rhgt', 'exchange']
    # short_names_x = ['time','tas', 'tas_norm','nppgt', 'fgco2gt', 'rhgt', 'exchange']
    # short_names_y = ['tas', 'tas_norm', 'nppgt',  'fgco2gt', 'rhgt', 'exchange']

    # jobtype = 'marine'
    short_names, short_names_x, short_names_y = [], [], []
    jobtype = 'debug'
    if jobtype == 'debug':
        short_names = [
                       'intpp',
                       'intppgt',
                       ]
        short_names_x = ['time', ]
        short_names_y = ['intpp', 'intppgt',]

    if jobtype == 'marine':
        short_names = ['tas', 'tas_norm', 'co2',
                #    'npp', 'nppgt', 'rhgt', 'exchange',
                #    'nppgt_norm','rhgt_norm','exchange_norm','fgco2gt_norm', 'intppgt_norm',
                       'intpp', 'fgco2', 'epc100', 'intdic', 'intpoc', 'fric', 'froc',
                       'intppgt','fgco2gt', 'epc100gt', 'intdicgt', 'intpocgt', 'fricgt', 'frocgt',
                       ]
        short_names_x = ['time', 'co2', 'tas', 'tas_norm', 'fgco2gt',]
        #'intpp', 'epc100', 'intdic', 'intpoc', 'fric', 'froc'] #'nppgt', 'fgco2gt', 'rhgt', 'exchange']
        #short_names_y = ['nppgt', 'nppgt_norm','rhgt_norm','exchange_norm','fgco2gt_norm', 'co2',]
        short_names_y = ['intpp', 'fgco2', 'epc100', 'intdic', 'intpoc', 'fric', 'froc', 'fgco2gt', 'intppgt','epc100gt', 'intdicgt', 'intpocgt', 'fricgt', 'frocgt',]

    if jobtype == 'full':
        short_names = ['tas', 'tas_norm', 'co2',
                       'npp', 'nppgt', 'rhgt', 'exchange',
                       'nppgt_norm','rhgt_norm','exchange_norm','fgco2gt_norm', 'intppgt_norm',
                       'intpp', 'fgco2', 'epc100', 'intdic', 'intpoc', 'fric', 'froc',
                       'intppgt','fgco2gt', 'epc100gt', 'intdicgt', 'intpocgt', 'fricgt', 'frocgt',
                       ]
        short_names_x = ['time', 'co2', 'tas', 'tas_norm', 'fgco2gt',
                         'intpp', 'epc100', 'intdic', 'intpoc', 'fric', 'froc', 'nppgt', 'fgco2gt', 'rhgt', 'exchange']
        short_names_y = ['nppgt', 'nppgt_norm','rhgt_norm','exchange_norm','fgco2gt_norm', 'co2',
                         'intpp', 'fgco2', 'epc100', 'intdic', 'intpoc', 'fric', 'froc', 'fgco2gt', 'intppgt','epc100gt', 'intdicgt', 'intpocgt', 'fricgt', 'frocgt',]

     # ]'npp', 'nppgt', 'intpp', 'intppgt_norm', 'fgco2gt', 'rhgt', 'exchange', 'nppgt_norm','rhgt_norm','exchange_norm','fgco2gt_norm']

    pairs = []

    for do_ma in [True, ]:#False]:
        data_dict = load_timeseries(cfg, short_names)
        thresholds_dict = load_thresholds(cfg, data_dict)

        for (short_name, exp, ensemble),cube  in sorted(data_dict.items()):
            if do_ma and short_name != 'co2':
                data_dict[(short_name, exp, ensemble)] = moving_average(cube, '21 years')

        print(short_names)
        for x in short_names_x:
            for y in short_names_y:
                if x == y: continue
                if (x,y) in pairs: continue
                print('main:', do_ma, x, y)
                make_ts_figure(cfg, data_dict, thresholds_dict, x=x, y=y,
                               markers='thresholds', do_moving_average=do_ma)
                pairs.append((x,y))
                pairs.append((y,x))

    logger.info('Success')



if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
