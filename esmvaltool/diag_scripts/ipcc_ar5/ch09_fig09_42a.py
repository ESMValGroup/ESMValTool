#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to plot figure 9.42a of IPCC AR5 chapter 9.

Description
-----------
Calculate and plot the equilibrium climate sensitivity (ECS) vs. the global
mean surface temperature (GMSAT) (see IPCC AR5 WG1 ch.9, fig. 9.42a).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
axes_functions : dict, optional
    Keyword arguments for the plot appearance functions.
dataset_style : str, optional
    Dataset style file (located in
    :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
matplotlib_style : str, optional
    Dataset style file (located in
    :mod:`esmvaltool.diag_scripts.shared.plot.styles_python.matplotlib`).
save : dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots).

"""

import logging
import os

import iris
import seaborn as sns

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    extract_variables,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    io,
    plot,
    run_diagnostic,
    variables_available,
)

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(project, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption':
        ('Equilibrium climate sensitivity (ECS) against the global '
         'mean surface temperature of {} models, both for the '
         'period 1961-1990 (larger symbols) and for the '
         'pre-industrial control runs (smaller symbols).'.format(project)),
        'statistics': ['mean'],
        'domains': ['global'],
        'authors': ['schlund_manuel'],
        'references': ['flato13ipcc'],
        'realms': ['atmos'],
        'themes': ['phys'],
        'ancestors':
        ancestor_files,
    }
    return record


def plot_data(cfg, hist_cubes, pi_cubes, ecs_cube):
    """Plot data."""
    x_data = []
    y_data = []
    dataset_names = []
    plot_kwargs = []

    # Collect data
    for dataset in hist_cubes:
        ecs = ecs_cube.extract(iris.Constraint(dataset=dataset))
        if ecs is None:
            raise ValueError(f"No ECS data for '{dataset}' available")

        # Historical data
        x_data.append(ecs.data)
        y_data.append(hist_cubes[dataset].data)
        dataset_names.append(dataset)
        plot_kwargs.append({
            'label': dataset,
            'linestyle': 'none',
            'markersize': 10,
        })

        # PiControl data
        x_data.append(ecs.data)
        y_data.append(pi_cubes[dataset].data)
        dataset_names.append(dataset)
        plot_kwargs.append({
            'label': '_' + dataset,
            'linestyle': 'none',
            'markersize': 6,
        })

    # Plot data
    path = get_plot_filename('ch09_fig09_42a', cfg)
    plot.multi_dataset_scatterplot(
        x_data,
        y_data,
        dataset_names,
        path,
        plot_kwargs=plot_kwargs,
        save_kwargs=cfg.get('save', {}),
        axes_functions=cfg.get('axes_functions', {}),
        dataset_style_file=cfg.get('dataset_style'),
        mpl_style_file=cfg.get('matplotlib_style'),
    )
    return path


def write_data(cfg, hist_cubes, pi_cubes, ecs_cube):
    """Write netcdf file."""
    datasets = []
    data_ecs = []
    data_hist = []
    data_pi = []
    for dataset in list(hist_cubes):
        ecs = ecs_cube.extract(iris.Constraint(dataset=dataset))
        if ecs is None:
            raise ValueError(f"No ECS data for '{dataset}' available")
        datasets.append(dataset)
        data_ecs.append(ecs.data)
        data_hist.append(hist_cubes[dataset].data)
        data_pi.append(pi_cubes[dataset].data)

    # Create cube
    dataset_coord = iris.coords.AuxCoord(datasets, long_name='dataset')
    tas_hist_coord = iris.coords.AuxCoord(data_hist,
                                          attributes={'exp': 'historical'},
                                          **extract_variables(
                                              cfg, as_iris=True)['tas'])
    tas_picontrol_coord = iris.coords.AuxCoord(data_pi,
                                               attributes={'exp': 'piControl'},
                                               **extract_variables(
                                                   cfg, as_iris=True)['tas'])
    cube = iris.cube.Cube(data_ecs,
                          var_name='ecs',
                          long_name='Equilibrium Climate Sensitivity (ECS)',
                          aux_coords_and_dims=[(dataset_coord, 0),
                                               (tas_hist_coord, 0),
                                               (tas_picontrol_coord, 0)])

    # Save file
    path = get_diagnostic_filename('ch09_fig09_42a', cfg)
    io.iris_save(cube, path)
    return path


def main(cfg):
    """Run the diagnostic."""
    sns.set(**cfg.get('seaborn_settings', {}))
    input_data = cfg['input_data'].values()
    project = list(group_metadata(input_data, 'project').keys())
    project = [p for p in project if 'obs' not in p.lower()]
    if len(project) == 1:
        project = project[0]

    # Check if tas is available
    if not variables_available(cfg, ['tas']):
        raise ValueError("This diagnostic needs 'tas' variable")

    # Get ECS data
    ecs_filepath = io.get_ancestor_file(cfg, 'ecs.nc')
    ecs_cube = iris.load_cube(ecs_filepath)

    # Create iris cubes for each dataset
    hist_cubes = {}
    pi_cubes = {}
    for data in input_data:
        name = data['dataset']
        logger.info("Processing %s", name)
        cube = iris.load_cube(data['filename'])

        # Preprocess cubes
        cube.convert_units(cfg.get('tas_units', 'celsius'))
        cube = cube.collapsed(['time'], iris.analysis.MEAN)

        # Save cubes
        if data.get('exp') == 'historical':
            hist_cubes[name] = cube
        elif data.get('exp') == 'piControl':
            pi_cubes[name] = cube
        else:
            pass

    # Plot data
    plot_path = plot_data(cfg, hist_cubes, pi_cubes, ecs_cube)

    # Write netcdf file
    netcdf_path = write_data(cfg, hist_cubes, pi_cubes, ecs_cube)

    # Provenance
    ancestor_files = [d['filename'] for d in input_data]
    ancestor_files.append(ecs_filepath)
    provenance_record = get_provenance_record(project, ancestor_files)
    provenance_record.update({
        'plot_file': plot_path,
        'plot_types': ['scatter'],
    })
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
