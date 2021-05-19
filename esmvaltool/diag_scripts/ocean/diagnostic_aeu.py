"""
Time series diagnostics
=======================

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
import iris.quickplot as qplt

import matplotlib.pyplot as plt
logging.getLogger('matplotlib.font_manager').disabled = True
import matplotlib.patches as patches

import cf_units

import numpy as np
import geopy.distance

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvalcore.preprocessor._volume import extract_volume
from esmvalcore.preprocessor._area import extract_region 


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def calculate_crossssection_mean(cube):
    #
    # calculate AEU
    #
    cube.data = np.ma.masked_invalid(cube.data)
    cube.data = np.ma.masked_where(cube.data.mask + (cube.data > 1.E10) + (cube.data <= 0.), cube.data)
    print('calculate_crossssection_mean: raw cube', cube.shape)
    cross_section = calc_cross_section(cube)
    print('calculate_crossssection_mean: X section:', cross_section.shape)

    cube =  cube.collapsed(['time', ], iris.analysis.MEAN)
    print('calculate_crossssection_mean: post-collapse', cube.shape)

    cube.data = cube.data * cross_section * 1.e-6 # convert to 1e6 m3/s (Sv)
    cube.units = cf_units.Unit('Sv')
    return cube


def calc_cross_section(cube):
    """
    Calculate the cross sectional area
    """
    lats = cube.coord('latitude')
    depth = cube.coord('depth')
    thicknesses = np.abs(depth.bounds[:,1] - depth.bounds[:,0])

    cross_section = np.zeros((len(depth.points), len(lats.points)))
    lats.guess_bounds()
    print(lats, depth)
    for l, la in enumerate(lats.points):
        # from@ https://stackoverflow.com/a/43211266/3082690
          
        coords_1 = (lats.bounds[l,0], -23.)
        coords_2 = (lats.bounds[l,1], -23.)
        distance = geopy.distance.distance(coords_1, coords_2).m
        cross_section[:, l] = thicknesses * distance
        print('calc_cross_section:', lats.bounds.shape, lats.bounds[l], 'distance:', distance,'m')

    cross_section = np.ma.array(cross_section)
    print('calc_cross_section', cross_section.shape, 'cs range:', (cross_section.min(), cross_section.max()) )
    return cross_section


def calculate_aeu(cube, cut='loose'):
    #
    # calculate AEU
    #
    # TODO: Add additional sensitivity tests.
    if cut == 'tight':
        cube = extract_volume(cube, 30., 300.)
        cube = cube.intersection(
            latitude=(-1.2, 1.2),
            ignore_bounds=True,
        )
       # cube = extract_region(cube, -25., -23., -1.2, 1.2)

    if cut == 'medium':
        cube = extract_volume(cube, 15., 400.)
        #cube = extract_region(cube, -25., -23., -3., 3.)
        cube = cube.intersection(
            latitude=(-3, 3.),
            ignore_bounds=True,
        )
    if cut == 'loose': pass

    cross_section = calc_cross_section(cube)

    cube.data = np.ma.masked_invalid(cube.data)
    cube.data = np.ma.masked_where(cube.data.mask + (cube.data > 1.E10) + (cube.data <= 0.), cube.data)
    times = cube.coord('time')
    data = cube.data.copy()

    for t,time in enumerate(times.points):
        print('calculat AEU', t, np.sum( data[t]*cross_section))
        data[t] = data[t]*cross_section * 1.e-6 # convert to 1e6 m3/s (Sv)
    cube.data = data
    #print(cube.data.max(),cube.units)
    cube.units = cf_units.Unit('Sv')
    #print(cube.data.max(),cube.units)
    cube =  cube.collapsed(['latitude', 'longitude', 'depth'], iris.analysis.SUM)

    return cube


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
    #if cubedata.shape:
    #    c
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

def make_time_series_sensitivity(
        cfg,
        metadata,
        filename,
):
    """
    Do a sensitivity analysis for the current.

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Making plots for each layer
    #for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
    cuts = ['loose', 'medium', 'tight']
    for cut in cuts:
        ncube = calculate_aeu(cube.copy(), cut=cut)
        timeplot(ncube, label=' '.join([metadata['dataset'], cut]))

    plt.axhline(14, c='k', label='Obs mean, Brandt et al (2021)')

    plt.legend()
    title = ' '.join([metadata['dataset'], 'AEU'] )
    plt.title(title)
    plt.legend(loc='best')
    plt.ylabel(str(cube.units))

    # Determine image filename:
    path = diagtools.get_image_path(
            cfg,
            metadata,
            suffix='sensitivity_timeseries'  + image_extention,
        )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()




def make_time_series_plots(
        cfg,
        metadata,
        filename,
        cut='loose',
        ma=None,
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

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Making plots for each layer
    #for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
    cube = calculate_aeu(cube,cut=cut)
    if ma:
        cube = moving_average(cube, ma)

    timeplot(cube, label=metadata['dataset'])

    title = ' '.join([metadata['dataset'], metadata['long_name'], cut, str(ma)])
    plt.title(title)
    plt.legend(loc='best')
    plt.ylabel(str(cube.units))

    # Determine image filename:
    path = diagtools.get_image_path(
            cfg,
            metadata,
            suffix='timeseries_'+cut+'_'+str(ma)  + image_extention,
        )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()



def make_transects_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple plot of the transect for an indivudual model.

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

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    cube = diagtools.make_depth_safe(cube)
    #cubes = diagtools.make_cube_region_dict(cube)
    cube = calculate_crossssection_mean(cube)

    # Determine y log scale.
    #set_y_logscale = determine_set_y_logscale(cfg, metadata)

    # Make a dict of cubes for each layer.
    # qplt.contourf(cube, 15, linewidth=0, rasterized=True)
    fig = plt.figure()
    ax = fig.add_subplot(111
)
    plt.contourf(cube.coord('latitude').points, 
                 -1.* np.abs(cube.coord('depth').points), 
                 cube.data,
                 15, linewidth=0, rasterized=True)
    plt.colorbar(label=str(cube.units))

    # tight 
    ax.add_patch(
        patches.Rectangle(
            xy=(-1.2, -300.), # lower left
            width=2.4,
            height=270.,
            linewidth=1,
            color='red',
            fill=False
        )
    )
    # medium
    ax.add_patch(
        patches.Rectangle(
            xy=(-3., -400.), # lower left
            width=6.,
            height=385.,
            linewidth=1,
            color='red',
            fill=False
        )
    )

    plt.xlabel('Latitude')
    plt.ylabel('Depth, m')
 
#    if set_y_logscale:
#        plt.axes().set_yscale('log')

#        if region:
#            region_title = region
#        else:
#            region_title = determine_transect_str(cube, region)

    # Add title to plot
    title = ' '.join(
        [metadata['dataset'], metadata['long_name'], ])
    plt.title(title)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    path = diagtools.get_image_path(
        cfg,
        metadata,
        suffix='transect' + image_extention,
    )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def multi_model_time_series(
        cfg,
        metadata,
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
        if metadata[filename]['frequency'] != 'fx':
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

            # Take a moving average, if needed.
            if 'moving_average' in cfg:
                cube = moving_average(model_cubes[filename][layer],
                                      cfg['moving_average'])
            else:
                cube = model_cubes[filename][layer]

            if 'MultiModel' in metadata[filename]['dataset']:
                timeplot(
                    cube,
                    c=color,
                    # label=metadata[filename]['dataset'],
                    ls=':',
                    lw=2.,
                )
                plot_details[filename] = {
                    'c': color,
                    'ls': ':',
                    'lw': 2.,
                    'label': metadata[filename]['dataset']
                }
            else:
                timeplot(
                    cube,
                    c=color,
                    # label=metadata[filename]['dataset'])
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
            if layer != '':
                if model_cubes[filename][layer].coords('depth'):
                    z_units = model_cubes[filename][layer].coord('depth').units
                else:
                    z_units = ''
        # Add title, legend to plots
        if layer:
            title = ' '.join([title, '(', str(layer), str(z_units), ')'])
        plt.title(title)
        plt.legend(loc='best')
        plt.ylabel(str(model_cubes[filename][layer].units))

        # Saving files:
        if cfg['write_plots']:
            path = diagtools.get_image_path(
                cfg,
                metadata[filename],
                prefix='MultipleModels_',
                suffix='_'.join(['timeseries',
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


def main(cfg):
    """
    Load the config file and some metadata, then pass them the plot making
    tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """

    cuts = ['loose', 'medium', 'tight']

    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('metadata filename:\t%s', metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)

        #######
        # Multi model time series
        #multi_model_time_series(
        #    cfg,
        #    metadatas,
        #)

        for filename in sorted(metadatas):
            if metadatas[filename]['frequency'] != 'fx':
                logger.info('-----------------')
                logger.info(
                    'model filenames:\t%s',
                    filename,
                )

                ######
                # Time series of individual model
                make_time_series_sensitivity(cfg, metadatas[filename], filename)
                for cut in cuts: 
                    for ma in [None, '10 years', '1 years']:
                        make_time_series_plots(cfg, metadatas[filename], filename, cut=cut, ma=ma)
                make_transects_plots(cfg, metadatas[filename], filename)
    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
