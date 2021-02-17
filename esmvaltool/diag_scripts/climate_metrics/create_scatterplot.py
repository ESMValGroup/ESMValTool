#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to create simple multi-dataset scatterplots.

Description
-----------
Create scatterplot for different datasets of a single variable. This diagnostic
needs preprocessed 1D cubes with single dimension ``dataset``.

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
dataset_style : str, optional
    Dataset style file (located in
    :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
pattern : str, optional
    Pattern to filter list of input data.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots).
y_range : list of float, optional
    Range for the y axis in the plot.

"""

import logging
import os

import iris
import matplotlib.pyplot as plt
import seaborn as sns

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    io,
    iris_helpers,
    plot,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(caption, ancestor_files, **kwargs):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'authors': ['schlund_manuel'],
        'references': ['acknow_project'],
        'ancestors': ancestor_files,
    }
    record.update(kwargs)
    return record


def plot_data(cfg, cube):
    """Create scatterplot for cube."""
    logger.debug("Plotting scatterplot for cube %s",
                 cube.summary(shorten=True))
    (_, axes) = plt.subplots()
    project = cube.attributes.get('project')

    # Plot
    for (idx, dataset_name) in enumerate(cube.coord('dataset').points):
        style = plot.get_dataset_style(dataset_name, cfg.get('dataset_style'))
        y_data = cube.extract(iris.Constraint(dataset=dataset_name)).data
        axes.plot(idx + 1,
                  y_data,
                  marker=style['mark'],
                  linestyle='none',
                  markeredgecolor=style['color'],
                  markerfacecolor=style['facecolor'],
                  label=dataset_name)

    # Plot appearance
    title = cube.long_name
    if project is not None:
        title += f' for {project}'
    axes.set_title(title)
    axes.tick_params(axis='x',
                     which='both',
                     bottom=False,
                     top=False,
                     labelbottom=False)
    axes.set_ylabel(f"{cube.var_name.upper()} / {cube.units}")
    axes.set_ylim(cfg.get('y_range'))
    legend = axes.legend(loc='center left',
                         bbox_to_anchor=[1.05, 0.5],
                         borderaxespad=0.0,
                         ncol=2)

    # Save plot
    plot_path = get_plot_filename(cube.var_name, cfg)
    plt.savefig(plot_path,
                orientation='landscape',
                bbox_inches='tight',
                additional_artists=[legend])
    logger.info("Wrote %s", plot_path)
    plt.close()
    return plot_path


def write_data(cfg, cube):
    """Write netcdf file."""
    cube.attributes.pop('provenance', None)
    netcdf_path = get_diagnostic_filename(cube.var_name, cfg)
    io.iris_save(cube, netcdf_path)
    return netcdf_path


def main(cfg):
    """Run the diagnostic."""
    sns.set(**cfg.get('seaborn_settings', {}))
    input_files = io.get_all_ancestor_files(cfg, pattern=cfg.get('pattern'))
    if len(input_files) != 1:
        raise ValueError(f"Expected exactly 1 file, got {len(input_files)}")
    input_file = input_files[0]
    logger.info("Found input file: %s", input_file)

    # Create plots
    cube = iris.load_cube(input_file)
    try:
        cube.coord('dataset')
    except iris.exceptions.CoordinateNotFoundError as exc:
        logger.error(
            "File '%s' does not contain necessary coordinate 'dataset'",
            input_file)
        raise exc

    # Sort coordinate 'dataset'
    [cube] = iris_helpers.intersect_dataset_coordinates([cube])

    # Create plot and netcdf file
    plot_path = plot_data(cfg, cube)
    netcdf_path = write_data(cfg, cube)

    # Provenance
    project = cube.attributes.get('project')
    caption = "{}{} for multiple datasets.".format(
        cube.long_name, '' if project is None else f' for {project}')
    provenance_record = get_provenance_record(caption, [input_file])
    provenance_record.update({
        'plot_file': plot_path,
        'plot_types': ['scatter'],
    })
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
