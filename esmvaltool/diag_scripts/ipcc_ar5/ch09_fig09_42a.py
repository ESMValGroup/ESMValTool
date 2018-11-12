#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Plot figure 9.42a of IPCC AR5 chapter 9.

###############################################################################
ipcc_ar5/ch09_fig09-42a.py
Author: Manuel Schlund (DLR, Germany)
CRESCENDO project
###############################################################################

Description
-----------
    Calculate and plot the equilibrium climate sensitivity (ECS) vs. the global
    mean surface temperature (GMSAT) for several CMIP5 models (see IPCC AR5 WG1
    ch. 9, fig. 9.42a).

Configuration options
---------------------
    ecs_filename        : Name of the netcdf in which the ECS data is saved.
    output_name         : Name of the output files.
    save                : Keyword arguments for the fig.saveplot() function.
    axes_functions      : Plot appearance functions.

###############################################################################

"""


import logging
import os
from datetime import datetime

import cf_units
import iris

import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))


def plot_data(cfg, datasets):
    """Plot data."""
    if not cfg[n.WRITE_PLOTS]:
        return
    filepath = os.path.join(cfg[n.PLOT_DIR],
                            cfg.get('output_name', 'fig09-42a') + '.' +
                            cfg[n.OUTPUT_FILE_TYPE])
    x_data = []
    y_data = []
    dataset_names = []
    plot_kwargs = []
    names = datasets.get_info_list(n.DATASET, short_name='ecs')
    ecs_data = datasets.get_data_list(short_name='ecs')

    # Historical
    x_data.extend(ecs_data)
    y_data.extend(datasets.get_data_list(short_name='tas', exp=HISTORICAL))
    dataset_names.extend(names)
    for name in names:
        plot_kwargs.append({'label': name, 'linestyle': 'none',
                            'markersize': 10})

    # piControl
    x_data.extend(ecs_data)
    y_data.extend(datasets.get_data_list(short_name='tas', exp=PICONTROL))
    dataset_names.extend(names)
    for name in names:
        plot_kwargs.append({'label': '_' + name, 'linestyle': 'none',
                            'markersize': 6})

    # Plot data
    e.plot.multi_dataset_scatterplot(
        x_data,
        y_data,
        dataset_names,
        filepath,
        plot_kwargs=plot_kwargs,
        save_kwargs=cfg.get('save', {}),
        axes_functions=cfg.get('axes_functions', {}))
    return


def write_data(cfg, datasets, variables):
    """Write netcdf file."""
    if cfg[n.WRITE_PLOTS]:
        data_ecs = datasets.get_data_list(short_name='ecs')
        data_tas_hist = datasets.get_data_list(short_name='tas',
                                               exp=HISTORICAL)
        data_tas_picontrol = datasets.get_data_list(short_name='tas',
                                                    exp=PICONTROL)
        models = datasets.get_info_list(n.DATASET, short_name='ecs')
        dataset_coord = iris.coords.AuxCoord(models, long_name='models')
        tas_hist_coord = iris.coords.AuxCoord(
            data_tas_hist,
            attributes={'experiment': HISTORICAL},
            **variables.iris_dict('tas'))
        tas_picontrol_coord = iris.coords.AuxCoord(
            data_tas_picontrol,
            attributes={'experiment': PICONTROL},
            **variables.iris_dict('tas'))
        attr = {'created_by': 'ESMValTool version {}'.format(cfg[n.VERSION]) +
                              ', diagnostic {}'.format(cfg[n.SCRIPT]),
                'creation_date': datetime.utcnow().isoformat(' ') + 'UTC'}
        cube = iris.cube.Cube(data_ecs, long_name=variables.long_name('ecs'),
                              var_name='ecs', units=variables.units('ecs'),
                              aux_coords_and_dims=[(dataset_coord, 0),
                                                   (tas_hist_coord, 0),
                                                   (tas_picontrol_coord, 0)],
                              attributes=attr)

        # Save file
        filepath = os.path.join(cfg[n.WORK_DIR],
                                cfg.get('output_name', 'fig09_42a') + '.nc')
        iris.save(cube, filepath)
        logger.info("Writing %s", filepath)


###############################################################################
# Setup diagnostic
###############################################################################

# Experiments
PICONTROL = 'piControl'
HISTORICAL = 'historical'
ABRUPT4XCO2 = 'abrupt4xCO2'
DIFF = 'difference of abrupt4xCO2 and piControl'

# Default settings
DEFAULT_TAS_UNITS = 'celsius'


def main(cfg):
    """Run the diagnostic.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.

    """
    ###########################################################################
    # Read recipe data
    ###########################################################################

    # Dataset data containers
    data = e.Datasets(cfg)
    logging.debug("Found datasets in recipe:\n%s", data)

    # Variables
    var = e.Variables(cfg)
    var.modify_var('tas', units=cfg.get('tas_units', DEFAULT_TAS_UNITS))
    logging.debug("Found variables in recipe:\n%s", var)

    # Get ECS data (ignore metadata.yml files)
    input_dirs = [d for d in cfg[n.INPUT_FILES]
                  if not d.endswith(n.METADATA_YAML_FILE)]
    if len(input_dirs) != 1:
        logging.error("Input files directory from ancestors should contain "
                      "exactly one directory (ECS directory)")
    ecs_filepath = os.path.join(input_dirs[0],
                                cfg.get('ecs_filename', 'ecs') + '.nc')

    ###########################################################################
    # Read data
    ###########################################################################

    # Create iris cube for each dataset
    for dataset_path in data:
        cube = iris.load(dataset_path, var.standard_names())[0]

        # Convert units if desired
        cube.convert_units(cfg.get('tas_units', DEFAULT_TAS_UNITS))

        # Total temporal means
        cube = cube.collapsed([n.TIME], iris.analysis.MEAN)
        data.set_data(cube.data, dataset_path)

    # Create iris cube for ECS data
    cube = iris.load_cube(ecs_filepath)
    var.add_vars(ecs={n.SHORT_NAME: cube.var_name,
                      n.LONG_NAME: cube.long_name,
                      n.UNITS: cube.units.format(cf_units.UT_DEFINITION)})
    for (idx, model) in enumerate(cube.coord('datasets').points):
        data.add_dataset('ecs_' + model,
                         data=cube.data[idx],
                         dataset=model,
                         short_name='ecs')

    ###########################################################################
    # Plot data
    ###########################################################################

    plot_data(cfg, data)

    ###########################################################################
    # Write nc file
    ###########################################################################

    write_data(cfg, data, var)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
