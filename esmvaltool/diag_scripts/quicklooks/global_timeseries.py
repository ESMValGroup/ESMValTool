"""
Global time series diagnostics
=======================

Description
-----------
Time series of global mean
In quicklook mode (cfg['quicklook']['active']: True) this diagnostic plots
the concatinated file.

Author
------
Lisa Bock (DLR, Germany)

Project
-------
CMIP6-DICAD

Configuration options in recipe
-------------------------------
time_int: min and max for time axis
y_min: min of y axis
y_max: max of y axis
multimodel_plot: if True: additional plot with all datasets
                 qicklook mode: all concatinated files
                 no quicklook mode: all dataset given in recipe

"""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools

from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename,
    run_diagnostic)

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(caption):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': 'times',
        'authors': ['bock_ls'],
        'references': ['acknow_project'],
    }
    return record


def compute_diagnostic(filename):
    """Compute global average."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    if cube.ndim > 1:
        grid_areas = iris.analysis.cartography.area_weights(cube)

        # from zonal mean to global mean
        cube = cube.collapsed('latitude', iris.analysis.MEAN,
                              weights=grid_areas)

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
    cubedata = np.ma.array(cube.data)
    if len(cubedata.compressed()) == 1:
        plt.axhline(cubedata.compressed(), **kwargs)
        return

    times = diagtools.cube_time_to_float(cube)
    plt.plot(times, cubedata, **kwargs)


def make_time_series_plots(
        cfg,
        metadata,
        filename,
):
    """
    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file

    """

    # Load and compute cube
    cube = compute_diagnostic(filename)

    logger.info(cube)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    timeplot(cube, label=metadata['dataset'])

    # Add title, legend to plots
    title = ' '.join([metadata['dataset'], metadata['long_name']])
    plt.title(title)
    #plt.legend(loc='best')
    plt.ylabel(str(cube.units))

    # set the y-limits
    if 'y_min' in cfg:
        plt.ylim(bottom=cfg['y_min'])
    if 'y_max' in cfg:
        plt.ylim(top=cfg['y_max'])

    #set time axis limits
    if 'time_int' in cfg:
        plt.xlim(cfg['time_int'])

    # Determine image filename:
    path = diagtools.get_image_path(
        cfg,
        metadata,
        prefix='Model',
        suffix='global_timeseries' + image_extention,
        metadata_id_list=[
            'dataset', 'field', 'short_name'
        ],
    )

    # Saving files:
    if cfg['write_plots']:

        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()

    # Write netcdf file for every plot
    diagname = '_'.join([metadata['dataset'], metadata['short_name'],
                         'global_timeseries'])
    diagnostic_file = get_diagnostic_filename(diagname, cfg)
    logger.info("Saving analysis results to %s", diagnostic_file)
    iris.save(cube, target=diagnostic_file)

    # Provenance
    provenance_record = get_provenance_record(
        "Timeseries of global mean of {} for dataset {}."
        .format(metadata['short_name'], metadata['dataset']))
    provenance_record.update({
        'plot_file': path,
    })

    return (diagnostic_file, provenance_record)

def multi_model_time_series(
        cfg,
        metadata,
):
    """
    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.

    """

    model_cube = {}
    for filename in sorted(metadata):
        # Load and compute cube
        cube = compute_diagnostic(filename)

        logger.info(cube)

        model_cube[filename] = cube

        metadata[filename]['dataset'] = cube.attributes['model_id']

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make the plot
    title = ''
    plot_details = {}
    cmap = plt.cm.get_cmap('viridis')

    # Plot each file in the group
    for index, filename in enumerate(sorted(metadata)):
        if len(metadata) > 1:
            color = cmap(index / (len(metadata) - 1.))
        else:
            color = 'blue'

        cube = model_cube[filename]

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
    # Add title, legend to plots
    plt.title(title)
    plt.legend(loc='best')
    plt.ylabel(str(model_cube[filename].units))

    # set the y-limits
    if 'y_min' in cfg:
        plt.ylim(bottom=cfg['y_min'])
    if 'y_max' in cfg:
        plt.ylim(top=cfg['y_max'])

    # Saving files:
    if cfg['write_plots']:
        path = diagtools.get_image_path(
            cfg,
            metadata[filename],
            prefix='MultiModel',
            suffix='global_timeseries' + image_extention,
            metadata_id_list=[
                'field', 'short_name'
            ],
        )

    # Resize and add legend outside thew axes.
    plt.gcf().set_size_inches(9., 6.)
    diagtools.add_legend_outside_right(
        plot_details, plt.gca(), column_width=0.15)

    logger.info('Saving plots to %s', path)
    plt.savefig(path)
    plt.close()

    # Write netcdf file for every plot
    diagname = '_'.join(['MultiModel',
                         metadata[filename]['short_name'],
                         'global_timeseries'])
    diagnostic_file = get_diagnostic_filename(diagname, cfg)
    logger.info("Saving analysis results to %s", diagnostic_file)
    iris.save(cube, target=diagnostic_file)

    # Provenance
    provenance_record = get_provenance_record(
        "Timeseries of global mean of {}."
        .format(metadata[filename]['short_name']))
    provenance_record.update({
        'plot_file': path,
    })

    return (diagnostic_file, provenance_record)


def main(cfg):
    """
    Load the config file and some metadata, then pass them the plot making
    tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """


    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('index:\t%s', index)
        logger.info('metadata filename:\t%s', metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)

        for filename in sorted(metadatas):

            logger.info('-----------------')
            logger.info(
                'preprocessed model filenames:\t%s',
                filename,
            )

            metadata = metadatas[filename]

            if cfg['quicklook']['active']:
                # if quicklook mode - plotting of concatinating file

                # path to concatinated file
                con_dir = cfg['quicklook']['output_dir'] + "/"
                filename = con_dir + '_'.join([metadata['dataset'],
                                               metadata['short_name'] + '.nc'])
                logger.info('concatinated filename:\t%s', filename)

            # Time series of individual model
            (path, provenance_record) = make_time_series_plots(
                cfg, metadata, filename)

            # Provenance
            if path is not None:
                provenance_record['ancestors'] = filename
                with ProvenanceLogger(cfg) as provenance_logger:
                    provenance_logger.log(path, provenance_record)

        if 'multimodel_plot' in cfg:
            if cfg['multimodel_plot']:

                if cfg['quicklook']['active']:
                    # if quicklook mode - plotting of concatinating file
                    con_dir = cfg['quicklook']['output_dir'] + "/"
                    con_files = [name for name in os.listdir(con_dir)
                                 if name.endswith(
                                     metadata['short_name'] + '.nc')]
                    con_files = [con_dir + name for name in con_files]

                    # if more datasets are given in concatinated files
                    if len(con_files) > 1:
                        meta_datas = {}
                        for filename in con_files:
                            meta_datas[filename] = dict(
                                short_name=metadata['short_name'],
                                long_name=metadata['long_name'])
                        logger.info('meta_datas:\t%s', meta_datas)

                        # Time series plot with all models
                        (path, provenance_record) = multi_model_time_series(
                            cfg, meta_datas)

                        # Provenance
                        if path is not None:
                            provenance_record['ancestors'] = con_files
                            with ProvenanceLogger(cfg) as provenance_logger:
                                provenance_logger.log(path, provenance_record)

                    else:
                        logger.info('Only one concatinated file available')

                else:

                    # if more datasets are given in the recipe
                    if len(metadatas) > 1:
                        # Time series plot with all models
                        (path, provenance_record) = multi_model_time_series(
                            cfg, metadatas)
                        # Provenance
                        if path is not None:
                            provenance_record['ancestors'] = metadatas
                            with ProvenanceLogger(cfg) as provenance_logger:
                                provenance_logger.log(path, provenance_record)


    logger.info('Success')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
