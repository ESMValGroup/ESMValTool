"""
Mpas's plotting tool. PAS
============================

Basically the same as time series, just with a few additions.

Diagnostic to produce figures of the time development of a field from
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

import iris
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from esmvalcore.preprocessor._time import extract_time

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

ipcc_colours={
    'historical':'blue',
    'ssp126': 'green',
    'ssp245': 'gold',
    'ssp370': 'orange',
    'ssp585': 'red',}

long_name_dict = {
    'thetao': 'Temperature',
    'tos': 'Surface Temperature',
    'tob': 'Seafloor Temperature',
    'sos': 'Surface Salinity',
    'uo': 'Zonal Velocity',
    'vo': 'Meridional Velocity',
    'ph': 'Surface pH',
    'chl': 'Surface chlorophyll',
    'zos': 'Sea Surface Height',
    'no3': 'Dissolved Nitrate',
    'o2': 'Dissolved Oxygen',
    'intpp': 'Integrated Primary production'}

def timeplot(cube, **kwargs):
    """
    Create a time series plot from the cube.

    Note that this function simple does the plotting, it does not save the
    image or do any of the complex work. This function also takes and of the
    key word arguments accepted by the matplotlib.plt.plot function.
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
    if window  in ['annual', ]:
        window = '1 year'
    window = window.split()
    window_len = int(window[0]) / 2.
    win_units = str(window[1])

    if win_units not in [
            'days', 'day', 'dy', 'months', 'month', 'mn', 'years', 'yrs',
            'year', 'yr'
    ]:
        raise ValueError("Moving average window units not recognised: " +
                         "{}".format(win_units))

    times = cube.coord('time').units.num2date(cube.coord('time').points)

    datetime = diagtools.guess_calendar_datetime(cube)

    output = []

    times = np.array([
        datetime(time_itr.year, time_itr.month, time_itr.day, time_itr.hour,
                 time_itr.minute) for time_itr in times
    ])

    for time_itr in times:
        if win_units in ['years', 'yrs', 'year', 'yr']:
            tmin = datetime(time_itr.year - window_len, time_itr.month,
                            time_itr.day, time_itr.hour, time_itr.minute)
            tmax = datetime(time_itr.year + window_len, time_itr.month,
                            time_itr.day, time_itr.hour, time_itr.minute)

        if win_units in ['months', 'month', 'mn']:
            tmin = datetime(time_itr.year, time_itr.month - window_len,
                            time_itr.day, time_itr.hour, time_itr.minute)
            tmax = datetime(time_itr.year, time_itr.month + window_len,
                            time_itr.day, time_itr.hour, time_itr.minute)

        if win_units in ['days', 'day', 'dy']:
            tmin = datetime(time_itr.year, time_itr.month,
                            time_itr.day - window_len, time_itr.hour,
                            time_itr.minute)
            tmax = datetime(time_itr.year, time_itr.month,
                            time_itr.day + window_len, time_itr.hour,
                            time_itr.minute)

        arr = np.ma.masked_where((times < tmin) + (times > tmax), cube.data)
        output.append(arr.mean())
    cube.data = np.array(output)
    return cube



def multi_model_time_series(
        cfg,
        metadatas,
        ts_dict = {},
        moving_average_str='',
        colour_scheme = 'IPCC',
        fig = None,
        ax = None,
        save = False
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
    model_cubes_paths = {}

    mean_cubes = {}
    mean_cubes_paths = {}

    details = []
    #metadata_list = []'short_name'
    for variable_group, filenames  in ts_dict.items():
        for fn in sorted(filenames):
            print(variable_group, fn)
            #assert 0 
            if metadatas[fn]['mip'] in ['Ofx', 'fx']: continue
            print('loading', fn)

            cube = iris.load_cube(fn)
            cube = diagtools.bgc_units(cube, metadatas[fn]['short_name'])

            model_cubes = add_dict_list(model_cubes, variable_group, cube)
            model_cubes_paths = add_dict_list(model_cubes_paths, variable_group, fn)

        # load means:
        work_dir = diagtools.folder([cfg['work_dir'], 'variable_group_means'])
        path = work_dir+'_'.join([variable_group, 'mean'])+'.nc'

        if os.path.exists(path):
            print('loading mean file:',variable_group, path)

            cube = iris.load_cube(path)
            mean_cubes[variable_group] = cube
        else:
            mean_cubes[variable_group] = diagtools.make_mean_of_cube_list(model_cubes[variable_group])
            print('saving mean file:',variable_group, path)
            iris.save(cube, path)

        metadatas[path] = metadatas[fn].copy()
        mean_cubes_paths[variable_group] = path

    if fig is None or ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        save = True
    else:
        plt.sca(ax)

    title = ''
    z_units = ''
    plot_details = {}
    if colour_scheme in ['viridis', 'jet']:
        cmap = plt.cm.get_cmap(colour_scheme)

    # Plot individual model in the group
    for variable_group, cubes in model_cubes.items():
        for i, cube in enumerate(cubes):
            path = model_cubes_paths[variable_group][i]
            metadata = metadatas[path]
            scenario = metadata['exp']
            dataset = metadata['dataset']

            if colour_scheme in ['viridis', 'jet']:
                if len(metadata) > 1:
                    color = cmap(index / (len(metadata) - 1.))
                else:
                    color = 'blue'
                label = dataset
            if colour_scheme in 'IPCC':
                color = ipcc_colours[scenario]
                label =  scenario
            # Take a moving average, if needed.
            if moving_average_str:
                    cube = moving_average(cube,
                          moving_average_str)

            timeplot(
                cube,
                c=color,
                ls='-',
                lw=0.6,
            )

    # Plot each mean file in the group
    for variable_group, cube in mean_cubes.items():
        path = mean_cubes_paths[variable_group]
        metadata = metadatas[path]
        scenario = metadata['exp']
        dataset = metadata['dataset']

        if colour_scheme in ['viridis', 'jet']:
            if len(metadata) > 1:
                color = cmap(index / (len(metadata) - 1.))
            else:
                color = 'blue'
            label = dataset
        if colour_scheme in 'IPCC':
            color = ipcc_colours[scenario]
            label =  scenario
        # Take a moving average, if needed.
        if moving_average_str:
                cube = moving_average(cube,
                                      moving_average_str)

        timeplot(
            cube,
            c=color,
            # label=metadata[filename]['dataset'])
            ls='-',
            lw=2.,
        )
        plot_details[path] = {
            'c': color,
            'ls': '-',
            'lw': 2.,
            'label': label,
        }

#   # Add title, legend to plots
    plt.title(title)
    plt.legend(loc='best')
    #plt.ylabel(str(model_cubes[filename][layer].units))

    # Saving files:
    if save:
        path = diagtools.folder(cfg['plot_dir']+'/individual_panes')
        path += '_'.join(['multi_model_ts',] )
        path += diagtools.get_image_format(cfg)
 
        # Resize and add legend outside thew axes.
        fig.set_size_inches(9., 6.)
        diagtools.add_legend_outside_right(
            plot_details, plt.gca(), column_width=0.15)

        logger.info('Saving plots to %s', path)
        plt.savefig(path)
        plt.close()
    else:
        return fig, ax



def multi_model_clim_figure(
        cfg,
        metadatas,
        short_name,
        figure_style = 'plot_all_years',
        hist_time_range = [1990., 2000.],
        ssp_time_range = [2040., 2050.],
        fig = None,
        ax = None,
        save=False
    ):
    """
    produce a monthly climatology figure.
    """
    if fig is None or ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        save=True

    labels = []
    for filename, metadata in metadatas.items():
        if short_name != metadata['short_name']:
            continue
        scenario = metadata['exp']
        cube = iris.load_cube(filename)
        cube = diagtools.bgc_units(cube, metadata['short_name'])

        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')

        if not cube.coords('month_number'):
            iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')

        if scenario == 'historical':
            cube = extract_time(cube, hist_time_range[0], 1, 1, hist_time_range[1], 12, 31)
        else:
            cube = extract_time(cube, ssp_time_range[0], 1, 1, ssp_time_range[1], 12, 31)

        months =  cube.coord('month_number').points
        years = cube.coord('year').points

        years_range = sorted({yr:True for yr in cube.coord('year').points}.keys())
        data = cube.data
        label = scenario
        if figure_style == 'plot_all_years':
            for year in years_range:
                t = np.ma.masked_where(years!= year, months)
                d =  np.ma.masked_where(years!= year, data)
                if int(year)%10==0: lw = 1
                else: lw = 0.4
                if label not in labels:
                    plt.plot(t, d, color = ipcc_colours[scenario], lw = lw,label=label )
                    labels.append(label)
                else:
                    plt.plot(t, d, color = ipcc_colours[scenario], lw = lw )

        if figure_style == 'mean_and_range':
            cube_mean = cube.copy().aggregated_by(['month_number', ], iris.analysis.MEAN)
            cube_min = cube.copy().aggregated_by(['month_number', ], iris.analysis.MIN)
            cube_max = cube.copy().aggregated_by(['month_number', ], iris.analysis.MAX)

            if label not in labels:
                plt.plot(cube_mean.coord('month_number').points, cube_mean.data, color = ipcc_colours[scenario], lw = 2., label=label)
                labels.append(label)
            else:
                plt.plot(cube_mean.coord('month_number').points, cube_mean.data, color = ipcc_colours[scenario], lw = 2.)

            plt. fill_between(cube_mean .coord('month_number').points,
                    cube_min.data,
                    cube_max.data,
                    alpha = 0.2,
                    color=ipcc_colours[scenario])

    plt.legend()

    time_str = '_'.join(['-'.join([str(t) for t in hist_time_range]), 'vs',
                         '-'.join([str(t) for t in ssp_time_range])])

    plt.suptitle(' '.join([long_name_dict[short_name], 'in Ascension'
                           ' Island MPA \n Historical', '-'.join([str(t) for t in hist_time_range]),
                           'vs SSP', '-'.join([str(t) for t in ssp_time_range]) ]))

    units = cube.units
    ax.set_xlabel('Months')
    ax.set_ylabel(' '.join([short_name+',', str(units)]))

    if save:
        # save and close.
        path = diagtools.folder(cfg['plot_dir']+'/clim')
        path += '_'.join(['multi_model_clim', short_name, figure_style, time_str])
        path += diagtools.get_image_format(cfg)

        logger.info('Saving plots to %s', path)
        plt.savefig(path)
        plt.close()
    else:
        return fig, ax







def do_gridspec():
    #fig = None
    fig = plt.figure()
    fig.set_size_inches(12., 9.) #, 6.)

    gs = matplotlib.gridspec.GridSpec(ncols=1, nrows=2)
    gs0 =gs[0,0].subgridspec(ncols=4, nrows=2, ) #ght_ratios=[2, 1], hspace=0.)
    subplots = {}
    subplots['timeseries'] = fig.add_subplot(gs0[0:2,0:2])
    subplots['climatology'] = fig.add_subplot(gs0[0, 2])
    subplots['clim_diff'] = fig.add_subplot(gs0[0, 3])
    subplots['profile'] = fig.add_subplot(gs0[1, 2])
    subplots['prof_diff'] = fig.add_subplot(gs0[1, 3])
    #
    # maps:
    gs1 =gs[1,0].subgridspec(ncols=3, nrows=2, )
    subplots['map_hist'] = fig.add_subplot(gs1[0,0])
    subplots['map_obs'] = fig.add_subplot(gs1[1,0])
    subplots['map_ssp125'] = fig.add_subplot(gs1[0,1])
    subplots['map_ssp245'] = fig.add_subplot(gs1[1,1])
    subplots['map_ssp379'] = fig.add_subplot(gs1[0,2])
    subplots['map_ssp585'] = fig.add_subplot(gs1[1,2])
    #plt.show()
    #plt.savefig('tmp.png')

    return fig, subplots


def add_dict_list(dic, key, item):
   try: dic[key].append(item)
   except: dic[key]= [item, ]
   return dic

def main(cfg):
    """
    Load the config file and some metadata, then pass them the plot making
    tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    # Here is plan the parts of the plot:
    # Plot A:
    # top left: time series
    # top right: climatology
    # middle: profile:
    # right: maps

    # Strategy is to get all the individual plots made.
    # individual plots of the ensemble mean
    # Then stitch them together into a big single plot.
    metadatas = diagtools.get_input_files(cfg, )

    time_series_fns = {}
    profile_fns = {}
    maps_fns = {}

    for fn, metadata in metadatas.items():
        #print(os.path.basename(fn),':',metadata['variable_group'])
        variable_group = metadata['variable_group']
        if variable_group.find('_ts_')>-1:
            time_series_fns = add_dict_list(time_series_fns, variable_group, fn)
        if variable_group.find('_profile_')>-1:
            profile_fns = add_dict_list(profile_fns, variable_group, fn)
        if variable_group.find('_map_')>-1:
            maps_fns = add_dict_list(maps_fns, variable_group, fn)
    # Individual plots - standalone
    multi_model_time_series(
            cfg,
            metadatas,
            ts_dict = time_series_fns,

            #moving_average_str='',
            #colour_scheme = 'viridis',
            #fig = fig,
            #ax =  subplots['timeseries'],
    )
    
    #print(time_series_fns)
    #print(profile_fns)
    #print(maps_fns)

    fig, subplots = do_gridspec()
    fig, subplots['timeseries'] = multi_model_time_series(
            cfg,
            metadatas,
            ts_dict = time_series_fns,
            #moving_average_str='',
            #colour_scheme = 'viridis',
            fig = fig,
            ax =  subplots['timeseries'],
    )

    # save and close.
    path = diagtools.folder(cfg['plot_dir']+'/whole_plot')
    path += '_'.join(['multi_model_whole_plot'])
    path += diagtools.get_image_format(cfg)

    logger.info('Saving plots to %s', path)
    plt.savefig(path)
    plt.close()


    assert 0


    assert 0
    #moving_average_str = cfg.get('moving_average', None)
    short_names = {}
    metadatas = diagtools.get_input_files(cfg, )
    for fn, metadata in metadatas.items():
        short_names[metadata['short_name']] = True

    for short_name in short_names.keys():
        hist_time_ranges = [[1850, 2015], [1850, 1900], [1950, 2000], [1990, 2000], [1990, 2000]]
        ssp_time_ranges  = [[2015, 2100], [2050, 2100], [2050, 2100], [2040, 2050], [2090, 2100]]
        for hist_time_range, ssp_time_range in zip(hist_time_ranges, ssp_time_ranges):

            for figure_style in ['plot_all_years', 'mean_and_range']:

                multi_model_clim_figure(
                    cfg,
                    metadatas,
                    short_name,
                    figure_style=figure_style,
                    hist_time_range=hist_time_range,
                    ssp_time_range=ssp_time_range,
                )
    return

    moving_average_strs = ['', 'annual', '5 years', '10 years', '20 years']
    for moving_average_str in moving_average_strs:
        for index, metadata_filename in enumerate(cfg['input_files']):
            logger.info('metadata filename:\t%s', metadata_filename)

            metadatas = diagtools.get_input_files(cfg, index=index)

            #######
            # Multi model time series
            multi_model_time_series(
                cfg,
                metadatas,

                moving_average_str=moving_average_str,
                colour_scheme = 'IPCC',
            )
            continue
            for filename in sorted(metadatas):
                if metadatas[filename]['frequency'] != 'fx':
                    logger.info('-----------------')
                    logger.info(
                        'model filenames:\t%s',
                        filename,
                    )

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
