"""
Profile diagnostics.
====================

Diagnostic to produce figure of the profile over time from a cube.
These plost show cube value (ie temperature) on the x-axis, and depth/height
on the y axis. The colour scale is the time series.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has a time component, and depth component, but no
latitude or longitude coordinates.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_profile:
      extract_volume:
        long1: 0.
        long2:  20.
        lat1:  -30.
        lat2:  30.
        z_min: 0.
        z_max: 3000.
      area_statistics:
        operator: mean


In order to add an observational dataset to the profile plot, the following
arguments are needed in the diagnostic script::

  diagnostics:
    diagnostic_name:
      variables:
        ...
      additional_datasets:
      - {observational dataset description}
      scripts:
        script_name:
          script: ocean/diagnostic_profiles.py
          observational_dataset: {observational dataset description}

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import logging
import os
import sys

import numpy as np
import iris
import iris.coord_categorisation
import iris.exceptions
import iris.quickplot as qplt
import matplotlib
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvalcore.preprocessor._time import extract_time
from esmvalcore.preprocessor._volume import extract_volume

from esmvalcore.preprocessor._regrid import extract_levels

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def determine_profiles_str(cube):
    """
    Determine a string from the cube, to describe the profile.

    Parameters
    ----------
    cube: iris.cube.Cube
        the opened dataset as a cube.

    Returns
    -------
    str
        Returns a string which describes the profile.
    """
    options = ['latitude', 'longitude']
    for option in options:
        coord = cube.coord(option)
        if len(coord.points) > 1:
            continue
        value = coord.points.mean()
        if option == 'latitude':
            return str(value) + ' N'
        if option == 'longitude':
            if value > 180.:
                return str(value - 360.) + ' W'
            return str(value) + ' E'
    return ''


def make_profiles_plots(
        cfg,
        metadata,
        filename,
        obs_metadata={},
        obs_filename='',
):
    """
    Make a profile plot for an individual model.

    The optional observational dataset can also be added.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.
    obs_metadata: dict
        The metadata dictionairy for the observational dataset.
    obs_filename: str
        The preprocessed observational dataset file.

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    try:
        raw_times = diagtools.cube_time_to_float(cube)
    except iris.exceptions.CoordinateNotFoundError:
        return

    # Make annual or Decadal means from:
    if np.max(raw_times) - np.min(raw_times) < 20:
        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
    else:
        cube = diagtools.decadal_average(cube)

    times_float = diagtools.cube_time_to_float(cube)
    time_0 = times_float[0]

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    cmap = plt.cm.get_cmap('jet')

    plot_details = {}
    for time_index, time in enumerate(times_float):
        if times_float[-1] == time_0:
            color = 'black'
        else:
            color = cmap((time - time_0) / (times_float[-1] - time_0))

        qplt.plot(cube[time_index, :], cube[time_index, :].coord('depth'),
                  c=color)

        plot_details[str(time_index)] = {'c': color, 'ls': '-', 'lw': 1,
                                         'label': str(int(time))}

    # Add observational data.
    if obs_filename:
        obs_cube = iris.load_cube(obs_filename)
        obs_cube = diagtools.bgc_units(obs_cube, metadata['short_name'])
        obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)

        obs_key = obs_metadata['dataset']
        qplt.plot(obs_cube, obs_cube.coord('depth'), c='black')

        plot_details[obs_key] = {'c': 'black', 'ls': '-', 'lw': 1,
                                 'label': obs_key}

    # Add title to plot
    title = ' '.join([
        metadata['dataset'],
        metadata['long_name'],
    ])
    plt.title(title)

    # Add Legend outside right.
    diagtools.add_legend_outside_right(plot_details, plt.gca())

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if multi_model:
        path = diagtools.folder(
            cfg['plot_dir']) + os.path.basename(filename).replace(
                '.nc', '_profile' + image_extention)
    else:
        path = diagtools.get_image_path(
            cfg,
            metadata,
            suffix='profile' + image_extention,
        )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def extract_depth_range(cube, target_depths, drange='surface', threshold=1000.):
    """
    Extract either the top 1000m or everything below there.

    drange can be surface or depths.

    Assuming everything is 1D at this stage!
    """
    if drange == 'surface':
        cube = extract_volume(cube.copy(), 0., threshold)
        target_depths = np.ma.masked_where(target_depths > threshold, target_depths).compressed()
    elif drange == 'depths':
        cube = extract_volume(cube.copy(), threshold, target_depths.max())
        target_depths = np.ma.masked_where(target_depths < threshold, target_depths).compressed()

    if threshold not in target_depths:
        target_depths = np.array(sorted(np.append(target_depths, threshold)))

    cube = extract_levels(cube, target_depths, "nearest_horizontal_extrapolate_vertical")
    print("extract_depth_range", target_depths, cube.coord('depth').points)
    return target_depths, cube


def make_multi_model_profiles_plots(
        cfg,
        metadatas,
        short_name,
        obs_metadata={},
        obs_filename='',
        hist_time_range = [1990., 2000.],
        ssp_time_range = [2040., 2050.],
        figure_style = 'difference',
        fig = None,
        gs0 = None,
        draw_legend=True
    ):
    """
    Make a profile plot for an individual model.

    The optional observational dataset can also be added.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.
    obs_metadata: dict
        The metadata dictionairy for the observational dataset.
    obs_filename: str
        The preprocessed observational dataset file.

    """
    # Load cube and set up units
    if fig is None:
        single_pane = True
        fig = plt.figure()
        gs = matplotlib.gridspec.GridSpec(ncols=1, nrows=1) # master
        gs0 =gs[0,0].subgridspec(ncols=1, nrows=2, height_ratios=[2, 1], hspace=0.)
    else:
        single_pane = False

    #gs0 = gs[0].subgridspec(2, 1, hspace=0.35) # scatters
    #gs1 = gs[1].subgridspec(3, 1, hspace=0.06 ) # maps
        #scatters
    ax0=fig.add_subplot(gs0[0,0]) # surface
    ax1=fig.add_subplot(gs0[1,0]) # depths

    cubes = {}
    for filename, metadata in metadatas.items():
        if short_name != metadata['short_name']:
            continue
        cube = iris.load_cube(filename)
        cube = diagtools.bgc_units(cube, metadata['short_name'])
        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')

        raw_times = diagtools.cube_time_to_float(cube)

        times_float = diagtools.cube_time_to_float(cube)
        dataset = metadata['dataset']
        scenario = metadata['exp']
        ensemble = metadata['ensemble']

        if scenario == 'historical':
            cube = extract_time(cube, hist_time_range[0], 1, 1, hist_time_range[1], 1, 1)
        else:
            cube = extract_time(cube, ssp_time_range[0], 1, 1, ssp_time_range[1], 1, 1)

        cube_mean = cube.copy().collapsed('time', iris.analysis.MEAN)
        cube_min = cube.copy().collapsed('time', iris.analysis.MIN)
        cube_max = cube.copy().collapsed('time', iris.analysis.MAX)

        cubes[(dataset, scenario, ensemble, 'mean')]  = cube_mean
        cubes[(dataset, scenario, ensemble, 'min')]  = cube_min
        cubes[(dataset, scenario, ensemble, 'max')]  = cube_max

    for (dataset, scenario, ensemble, metric), cube_mean in cubes.items():
        if metric != 'mean': continue

        cube_min = cubes[(dataset, scenario, ensemble, 'min')]
        cube_max = cubes[(dataset, scenario, ensemble, 'max')]
        cube_hist = cubes[(dataset, 'historical', ensemble, 'mean')]

        color = diagtools.ipcc_colours[scenario]
        plot_details = {}

        for ax, drange in zip([ax0, ax1], ['surface', 'depths']):

            if figure_style == 'difference':
                z_hist, z_cube_hist = extract_depth_range(cube_hist, cube_hist.coord('depth').points, drange=drange)
                z_min, z_cube_min = extract_depth_range(cube_min, cube_hist.coord('depth').points, drange=drange)
                z_max, z_cube_max = extract_depth_range(cube_max, cube_hist.coord('depth').points, drange=drange)
                z_mean, z_cube_mean = extract_depth_range(cube_mean, cube_hist.coord('depth').points, drange=drange)

                ax.plot(z_cube_mean.data - z_cube_hist.data, -1.*z_hist,
                    lw=2,
                    c=color,
                    label= scenario)

                ax.fill_betweenx(-1.*z_hist,
                    z_cube_min.data - z_cube_hist.data,
                    x2=z_cube_max.data - z_cube_hist.data,
                    alpha = 0.2,
                    color=color)

            if figure_style == 'compare':
                z_min, z_cube_min = extract_depth_range(cube_min, cube_mean.coord('depth').points, drange=drange)
                z_max, z_cube_max = extract_depth_range(cube_max, cube_mean.coord('depth').points, drange=drange)
                z_mean, z_cube_mean = extract_depth_range(cube_mean, cube_mean.coord('depth').points, drange=drange)

                ax.plot(z_cube_mean.data, -1.*z_mean,
                    lw=2,
                    c=color,
                    label=scenario)

                ax.fill_betweenx(-1.*z_mean,
                    z_cube_min.data,
                    x2=z_cube_max.data,
                    alpha = 0.2,
                    color=color)

        plot_details[scenario] = {'c': color, 'ls': '-', 'lw': 1,
                                             'label': scenario}

    # Add units:
    units = cubes[(dataset, scenario, ensemble,'mean')].units
    if figure_style == 'difference':
        ax.set_xlabel(r'$\Delta$ ' +str(units))
    if figure_style == 'compare':
        ax.set_xlabel(str(units))

    # Add observational data.
    if obs_filename:
        obs_cube = iris.load_cube(obs_filename)
        obs_cube = diagtools.bgc_units(obs_cube, metadata['short_name'])
        obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)

        obs_key = obs_metadata['dataset']
        qplt.plot(obs_cube, obs_cube.coord('depth'), c='black')

        plot_details[obs_key] = {'c': 'black', 'ls': '-', 'lw': 1,
                                 'label': obs_key}

    time_str = '_'.join(['-'.join([str(t) for t in hist_time_range]), 'vs',
                         '-'.join([str(t) for t in ssp_time_range])])

    # set x axis limits:
    xlims = np.array([ax0.get_xlim(), ax1.get_xlim()])
    ax0.set_xlim([xlims.min(), xlims.max()])
    ax1.set_xlim([xlims.min(), xlims.max()])

    ylims = np.array([ax0.get_ylim(), ax1.get_ylim()])
    ax0.set_ylim([-1000., ylims.max()])
    ax1.set_ylim([ylims.min(), -1001.])

    # hide between pane axes and ticks:
    ax0.spines['bottom'].set_visible(False)
    ax0.xaxis.set_ticks_position('none')
    ax0.xaxis.set_ticklabels([])
    ax1.spines['top'].set_visible(False)

#    ax0.fill_bewteen([xlims.min(), xlims.max()],[-1000., -1000.], [-975., -975.], color= 'k', alpha=0.2)
#    ax1.fill_bewteen([xlims.min(), xlims.max()],[-1100., -1100.],[-1000., -1000.], color= 'k', alpha=0.2)

    # draw a line between figures
    ax0.axhline(-999., ls='--', lw=1.5, c='black')
    ax1.axhline(-1001., ls='--', lw=1.5, c='black')

    if single_pane:
        # Add title to plot
        title = ' '.join([
            short_name,
            figure_style,
            time_str
            ])
        ax0.set_title(title)
    else:
        #ax0.yaxis.set_ticks_position('both')
        #ax1.yaxis.set_ticks_position('both')
        if figure_style == 'difference':
           ax0.yaxis.set_ticks_position('right')
           ax1.yaxis.set_ticks_position('right')
           ax0.yaxis.set_ticklabels([])
           ax1.yaxis.set_ticklabels([])
           ax0.set_title('Difference against historical')

        else:
             ax0.set_title('Mean and Range')


    # Add Legend outside right.
    if draw_legend:
        ax1.legend(loc='lower left')
    #diagtools.add_legend_outside_right(plot_details, plt.gca())

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    path = diagtools.get_image_path(
            cfg,
            metadata,
            prefix='_'.join(['multi_model', short_name, figure_style, str(time_str)]),
            suffix='profile' + image_extention,
        )

    if not single_pane:
        return fig
    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def make_multi_model_profiles_plot_pair(
        cfg,
        metadatas,
        short_name,
        obs_metadata={},
        obs_filename='',
        hist_time_range = [1990., 2000.],
        ssp_time_range = [2040., 2050.],

    ):

    fig = plt.figure()
    fig.set_size_inches(10,6)
    gs = matplotlib.gridspec.GridSpec(ncols=2, nrows=1, wspace=0.01) # master

    gs0 =gs[0, 0].subgridspec(ncols=1, nrows=2, height_ratios=[2, 1], hspace=0.)
    gs1 =gs[0, 1].subgridspec(ncols=1, nrows=2, height_ratios=[2, 1], hspace=0.)

    for figure_style, gs_spot, draw_legend in zip(['compare', 'difference'], [gs0, gs1], [True, False]):
        fig = make_multi_model_profiles_plots(
            cfg,
            metadatas,
            short_name,
            figure_style=figure_style,
            obs_metadata=obs_metadata,
            obs_filename=obs_filename,
            hist_time_range=hist_time_range,
            ssp_time_range=ssp_time_range,
            fig=fig,
            gs0=gs_spot,
            draw_legend=draw_legend,
        )

    time_str = '_'.join(['-'.join([str(t) for t in hist_time_range]), 'vs', 
                         '-'.join([str(t) for t in ssp_time_range])])

    long_name_dict = {
        'thetao': 'Temperature',
        'uo': 'Zonal Velocity',
        'vo': 'Meridional Velocity',
        'no3': 'Dissolved Nitrate',
        'o2': 'Dissolved Oxygen',}

    plt.suptitle(' '.join([long_name_dict[short_name], 'Profile in Ascension'
                           ' Island MPA \n Historical', '-'.join([str(t) for t in hist_time_range]),
                           'vs SSP', '-'.join([str(t) for t in ssp_time_range]) ]))

    # Determine image filename:
    path = diagtools.folder(cfg['plot_dir'])
    path += '_'.join(['multi_model_pair', short_name, str(time_str)])
    path += diagtools.get_image_format(cfg)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def main(cfg):
    """
    Run the diagnostics profile tool.

    Load the config file, find an observational dataset filename,
    pass loaded into the plot making tool.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    metadatas = diagtools.get_input_files(cfg)
    obs_filename = ''
    obs_key = 'observational_dataset'
    obs_metadata = {}

    short_names = {}
    for fn, metadata in metadatas.items():
        short_names[metadata['short_name']] = True

    if obs_key in cfg:
        obs_filename = diagtools.match_model_to_key(obs_key,
                                                    cfg[obs_key],
                                                    metadatas)
        if obs_filename:
            obs_metadata = metadatas[obs_filename]
        else:
            obs_metadata = ''

    hist_time_ranges = [[1850, 2015], [1850, 1900], [1950, 2000], [1990, 2000], [1990, 2000]] 
    ssp_time_ranges  = [[2015, 2100], [2050, 2100], [2050, 2100], [2040, 2050], [2090, 2100]]

    for hist_time_range, ssp_time_range in zip(hist_time_ranges, ssp_time_ranges):
        for short_name in short_names.keys():
            make_multi_model_profiles_plot_pair(
                cfg,
                metadatas,
                short_name,
                obs_metadata=obs_metadata,
                obs_filename=obs_filename,
                hist_time_range=hist_time_range,
                ssp_time_range=ssp_time_range,
            )
            for figure_style in ['compare', 'difference']:
                continue
                make_multi_model_profiles_plots(
                    cfg,
                    metadatas,
                    short_name,
                    figure_style=figure_style,
                    obs_metadata=obs_metadata,
                    obs_filename=obs_filename,
                    hist_time_range=hist_time_range,
                    ssp_time_range=ssp_time_range,
                )
    return

    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('metadata filename:\t%s', metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)

        for filename in sorted(metadatas.keys()):

            if filename == obs_filename:
                continue

            if metadatas[filename]['frequency'] == 'fx':
                continue

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            ######
            # Time series of individual model
            make_profiles_plots(cfg, metadatas[filename], filename,
                                obs_metadata=obs_metadata,
                                obs_filename=obs_filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
