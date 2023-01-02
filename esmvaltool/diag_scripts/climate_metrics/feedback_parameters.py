#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to calculate several climate feedback parameters and ECS.

Description
-----------
Calculate climate feedback parameters and ECS similar to Andrews et al. (2012)
using the regression methods proposed by Gregory et al. (2004).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
calculate_mmm : bool, optional (default: True)
    Calculate multi-model means.
only_consider_mmm : bool, optional (default: False)
    Only consider multi-model mean dataset. This automatically sets
    ``calculate_mmm`` to ``True``. For large multi-dimensional datasets, this
    might significantly reduce the computation time if only the multi-model
    mean dataset is relevant.
output_attributes : dict, optional
    Write additional attributes to netcdf files.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots).

"""

import logging
import os
from collections import OrderedDict
from copy import deepcopy
from functools import partial

import cf_units
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    extract_variables,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    io,
    plot,
    run_diagnostic,
    select_metadata,
    sorted_metadata,
    variables_available,
)

logger = logging.getLogger(os.path.basename(__file__))

EXP_4XCO2 = {
    'CMIP5': 'abrupt4xCO2',
    'CMIP6': 'abrupt-4xCO2',
}
FEEDBACK_PARAMETERS = {
    'rtnt': 'Net',
    'rtmt': 'Net (rtmt)',
    'rlntcs': 'LW (cs)',
    'rsntcs': 'SW (cs)',
    'netcre': 'Net CRE',
    'lwcre': 'LW CRE',
    'swcre': 'SW CRE',
}
VAR_NAMES = {
    'rtnt': 'lambda_net',
    'rtmt': 'lambda_net_rtmt',
    'rlntcs': 'lambda_lwcs',
    'rsntcs': 'lambda_swcs',
    'netcre': 'lambda_netcre',
    'lwcre': 'lambda_lwcre',
    'swcre': 'lambda_swcre',
}
UNITS = {
    'rtnt': 'W m-2 K-1',
    'rtmt': 'W m-2 K-1',
    'rlntcs': 'W m-2 K-1',
    'rsntcs': 'W m-2 K-1',
    'netcre': 'W m-2 K-1',
    'lwcre': 'W m-2 K-1',
    'swcre': 'W m-2 K-1',
    'F': 'W m-2',
    'ECS': 'K',
}
LONG_NAMES = {
    'rtnt': 'Net Climate Feedback Parameter',
    'rtmt': 'Net Climate Feedback Parameter (using rtmt)',
    'rlntcs': 'Long Wave (clear sky) Climate Feedback Parameter',
    'rsntcs': 'Short Wave (clear sky) Climate Feedback Parameter',
    'netcre': 'Net Cloud Radiative Effect Climate Feedback Parameter',
    'lwcre': 'Long Wave Cloud Radiative Effect Climate Feedback Parameter',
    'swcre': 'Short Wave Cloud Radiative Effect Climate Feedback Parameter',
    'F': 'Radiative Forcing',
    'ECS': 'Equilibrium Climate Sensitivity',
}
SHORTER_NAMES = {
    'rtnt': 'Net Feedback Parameter',
    'rtmt': 'Net Feedback Parameter (rtmt)',
    'rlntcs': 'LW (clear sky) Feedback Parameter',
    'rsntcs': 'SW (clear sky) Feedback Parameter',
    'netcre': 'Net CRE Feedback Parameter',
    'lwcre': 'LW CRE Feedback Parameter',
    'swcre': 'SW CRE Feedback Parameter',
    'F': 'Radiative Forcing',
    'ECS': 'Equilibrium Climate Sensitivity',
}
NICE_SYMBOLS = {
    'rtnt': r'$\lambda_{net}$',
    'rtmt': r'$\lambda_{net(rtmt)}$',
    'rlntcs': r'$\lambda_{lwcs}$',
    'rsntcs': r'$\lambda_{swcs}$',
    'netcre': r'$\lambda_{netcre}$',
    'lwcre': r'$\lambda_{lwcre}$',
    'swcre': r'$\lambda_{swcre}$',
}
NICE_UNITS = {
    'W m-2 K-1': r'$W m^{-2} K^{-1}$',
    'W m-2': r'$W m^{-2}$',
    'K': r'$K$',
}
RTMT_TEXT = (
    "For datasets {}, 'rtmt' (net top of model radiation) instead of 'rtnt' "
    "(net top of atmosphere radiation) is used due to lack of data. These two "
    "variables might differ.")

# Global variables which will be adapted during runtime
NDIMS = {}
COORDS = {}
RTMT_DATASETS = set()


def calculate_anomaly(data_4x, data_pic):
    """Calculate anomaly cube for a dataset."""
    # Abrupt 4xCO2
    cube_4x = iris.load_cube(data_4x[0]['filename'])
    iris.coord_categorisation.add_year(cube_4x, 'time')
    cube_4x = cube_4x.aggregated_by('year', iris.analysis.MEAN)

    # PiControl
    cube_pic = iris.load_cube(data_pic[0]['filename'])
    iris.coord_categorisation.add_year(cube_pic, 'time')
    cube_pic = cube_pic.aggregated_by('year', iris.analysis.MEAN)

    # Anomaly
    x_data = cube_pic.coord('year').points
    y_data = _get_data_time_last(cube_pic)
    slope = _get_slope(x_data, y_data)
    intercept = _get_intercept(x_data, y_data)
    for _ in range(cube_pic.ndim - 1):
        x_data = np.expand_dims(x_data, -1)
    new_x_data = np.broadcast_to(x_data, cube_pic.shape)
    new_data = slope * new_x_data + intercept
    cube_4x.data -= np.ma.masked_invalid(new_data)
    return cube_4x


def _check_array_shapes(list_of_arrays, var):
    """Check if list of arrays has identical shapes."""
    shapes = {a.shape for a in list_of_arrays}
    if len(shapes) > 1:
        raise ValueError(
            f"Expected cubes with identical shapes for multi-model mean "
            f"calculation of '{var}', got {shapes}")


def _check_cube_dimensions(cube):
    """Check if dimension of :class:`iris.cube.Cube` is valid."""
    var = 'tas' if cube.var_name == 'tas' else 'rad'
    ndim = cube.ndim
    if ndim not in (1, 2, 3):
        raise ValueError(
            f"This diagnostic supports only 1D (time), 2D (time + other "
            f"dimension) or 3D (time + 2 other dimensions) input data, got "
            f"{ndim}D data")
    coord_0_name = cube.coord(dimensions=(0, ), dim_coords=True).name()
    if coord_0_name != 'time':
        raise ValueError(
            f"This diagnostic expects 'time' as dimension 0 for every cube, "
            f"got '{coord_0_name}' for cube\n{cube}")
    if NDIMS.get(var) is None:
        NDIMS[var] = ndim
    else:
        if ndim != NDIMS[var]:
            raise ValueError(
                "This diagnostic supports only '{}' data with consistent "
                "numbers of dimension, got mixed data".format(
                    'radiation' if var is None else var))
    if NDIMS.get('tas', 0) > NDIMS.get('rad', 4):
        raise ValueError(
            f"This diagnostic expects temperature data to have smaller number "
            f"of dimensions than radiation data, got {NDIMS['tas']}D data for "
            f"temperature and {NDIMS['rad']}D data for radiation")
    if ndim < 2:
        return
    coords = [
        coord.name() for coord in cube.coords(dim_coords=True)
        if coord.name() != 'time'
    ]
    if COORDS.get(var) is None:
        COORDS[var] = coords
    else:
        if coords != COORDS[var]:
            raise ValueError(
                "This diagnostic expects identical coordinate names for all "
                "'{}' data, got {} and {}".format(
                    'radiation' if var is None else var, coords, COORDS[var]))


def _create_feedback_file(feedback_cube, dataset_name, cfg, description=None):
    """Save feedback parameter plot vs. remaining dimensions."""
    var = feedback_cube.var_name
    filename = ('{}_vs_{}_{}'.format(VAR_NAMES.get(var, var),
                                     '-'.join(COORDS['rad']), dataset_name))
    attrs = {
        'dataset': dataset_name,
        'radiation_variable': var,
    }
    attrs.update(cfg.get('output_attributes', {}))
    if description is not None:
        attrs['description'] = description
        filename += f"_{description.replace(' ', '_')}"
    feedback_cube.var_name = VAR_NAMES.get(var, var)
    feedback_cube.long_name = LONG_NAMES.get(var, var)
    feedback_cube.units = UNITS.get(var, 'unknown')
    feedback_cube.attributes = attrs

    # Write cube
    netcdf_path = get_diagnostic_filename(filename, cfg)
    io.iris_save(feedback_cube, netcdf_path)

    # Caption
    caption = (
        'Dependence of {} on {} for {}. The calculation follows Andrews et '
        'al., Geophys. Res. Lett., 39 (2012): The {} is defined as the '
        'slope of the linear regression between {}-dependent {} TOA radiance '
        'and the {} surface temperature anomaly{} of the abrupt 4x CO2 '
        'experiment.'.format(
            LONG_NAMES.get(var,
                           var), ' and '.join(COORDS['rad']), dataset_name,
            LONG_NAMES.get(var, var), ' and '.join(COORDS['rad']),
            FEEDBACK_PARAMETERS.get(var, var),
            ('global mean' if NDIMS.get('tas') == 1 else '{}-dependent'.format(
                ' and '.join(COORDS['tas']))),
            '' if description is None else f' for {description}'))
    return (netcdf_path, caption)


def _create_feedback_plot(tas_cube, cube, dataset_name, cfg, description=None):
    """Plot feedback parameter vs. remaining dimensions."""
    var = cube.var_name
    logger.debug("Plotting '%s' vs. %s for '%s'", SHORTER_NAMES.get(var, var),
                 COORDS['rad'], dataset_name)
    x_data = _get_data_time_last(tas_cube)
    y_data = _get_data_time_last(cube)
    coords = [(coord, idx - 1)
              for (idx, coord) in enumerate(cube.coords(dim_coords=True))
              if coord.name() != 'time']
    feedback_cube = iris.cube.Cube(_get_slope(x_data, y_data),
                                   var_name=var,
                                   dim_coords_and_dims=coords,
                                   units='W m-2 K-1')

    # Plot
    if feedback_cube.ndim == 1:
        iplt.plot(feedback_cube)
        plt.xlabel(f"{COORDS['rad'][0]} / "
                   f"{cube.coord(COORDS['rad'][0]).units.origin}")
        plt.ylabel(f"{NICE_SYMBOLS.get(var, var)} / "
                   f"{NICE_UNITS.get(feedback_cube.units.origin, 'unknown')}")
        colorbar = None
    elif feedback_cube.ndim == 2:
        iplt.contourf(feedback_cube, cmap='bwr', levels=_get_levels())
        colorbar = plt.colorbar(orientation='horizontal')
        colorbar.set_label(
            f"{NICE_SYMBOLS.get(var, var)} / "
            f"{NICE_UNITS.get(feedback_cube.units.origin, 'unknown')}")
        ticks = [-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0]
        colorbar.set_ticks(ticks)
        colorbar.set_ticklabels([str(tick) for tick in ticks])
        if COORDS['rad'] == ['latitude', 'longitude']:
            plt.gca().coastlines()
        else:
            plt.xlabel(f"{COORDS['rad'][0]} / "
                       f"{cube.coord(COORDS['rad'][0]).units.origin}")
            plt.ylabel(f"{COORDS['rad'][1]} / "
                       f"{cube.coord(COORDS['rad'][1]).units.origin}")
    else:
        raise ValueError(f"Cube dimension {feedback_cube.ndim} not supported")

    # Appearance
    title = f'{SHORTER_NAMES.get(var, var)} for {dataset_name}'
    filename = ('{}_vs_{}_{}'.format(VAR_NAMES.get(var, var),
                                     '-'.join(COORDS['rad']), dataset_name))
    if description is not None:
        title += f' ({description})'
        filename += f"_{description.replace(' ', '_')}"
    plt.title(title)
    plot_path = get_plot_filename(filename, cfg)
    plt.savefig(plot_path,
                bbox_inches='tight',
                orientation='landscape',
                additional_artists=[colorbar])
    logger.info("Wrote %s", plot_path)
    plt.close()

    return (plot_path, feedback_cube)


def _create_regression_file(tas_cube,
                            cube,
                            dataset_name,
                            cfg,
                            description=None):
    """Save regression plot as netcdf file for a given dataset."""
    var = cube.var_name
    reg = stats.linregress(tas_cube.data, cube.data)
    filename = f'{var}_regression_{dataset_name}'
    attrs = {
        'dataset': dataset_name,
        'regression_r_value': reg.rvalue,
        'regression_slope': reg.slope,
        'regression_interception': reg.intercept,
        'feedback_parameter': reg.slope,
    }
    attrs.update(cfg.get('output_attributes', {}))
    if description is not None:
        attrs['description'] = description
        filename += f"_{description.replace(' ', '_')}"
    if var in ('rtmt', 'rtnt'):
        attrs['ECS'] = -reg.intercept / (2.0 * reg.slope)
    tas_coord = iris.coords.AuxCoord(
        tas_cube.data,
        **extract_variables(cfg, as_iris=True)['tas'])
    cube = iris.cube.Cube(cube.data,
                          attributes=attrs,
                          aux_coords_and_dims=[(tas_coord, 0)],
                          **extract_variables(cfg, as_iris=True)[var])
    netcdf_path = get_diagnostic_filename(filename, cfg)
    io.iris_save(cube, netcdf_path)
    return netcdf_path


def _create_regression_plot(tas_cube,
                            cube,
                            dataset_name,
                            cfg,
                            description=None):
    """Create regression plot."""
    var = cube.var_name
    logger.debug("Plotting '%s' vs. 'tas' for '%s'", var, dataset_name)
    reg = stats.linregress(tas_cube.data, cube.data)

    # Regression line
    x_reg = np.linspace(-1.0, 9.0, 2)
    y_reg = reg.slope * x_reg + reg.intercept

    # Plot data
    title = (f'{FEEDBACK_PARAMETERS.get(var,var)} TOA radiance for '
             f'{dataset_name}')
    filename = f'{var}_regression_{dataset_name}'
    if description is not None:
        title += f' ({description})'
        filename += f"_{description.replace(' ', '_')}"
    plot_path = get_plot_filename(filename, cfg)
    text = r'r = {:.2f}, {} = {:.2f}'.format(reg.rvalue,
                                             NICE_SYMBOLS.get(var, var),
                                             reg.slope)
    if var in ('rtmt', 'rtnt'):
        text += ', F = {:.2f}, ECS = {:.2f}'.format(
            reg.intercept, -reg.intercept / (2.0 * reg.slope))
    plot.scatterplot(
        [tas_cube.data, x_reg],
        [cube.data, y_reg],
        plot_path,
        plot_kwargs=[{
            'linestyle': 'none',
            'markeredgecolor': 'b',
            'markerfacecolor': 'none',
            'marker': 's',
        }, {
            'color': 'k',
            'linestyle': '-',
        }],
        save_kwargs={
            'bbox_inches': 'tight',
            'orientation': 'landscape',
        },
        axes_functions={
            'axhline': {
                'args': [0.0],
                'kwargs': {
                    'color': 'black',
                    'linestyle': '--',
                },
            },
            'set_title': title,
            'set_xlabel':
            f"tas / {NICE_UNITS.get(tas_cube.units.origin, 'K')}",
            'set_ylabel':
            f"{var} / {NICE_UNITS.get(cube.units.origin, 'unknown')}",
            'set_xlim': [0.0, 8.0],
            'text': {
                'args': [0.05, 0.9, text],
                'kwargs': {
                    'transform': 'transAxes',
                },
            },
        },
    )
    return (plot_path, reg)


def _create_table(table, cfg, description=None):
    """Create summary table containing all climate feedback parameters."""
    logger.debug("Creating summary table")
    (cell_data, row_labels, col_labels, col_units) = _dict_to_array(table)

    # Create netcdf file
    cubes = _get_cube_list_for_table(cell_data, row_labels, col_labels,
                                     col_units)
    filename = 'summary_table'
    if description is not None:
        filename += f"_{description.replace(' ', '_')}"
    netcdf_path = get_diagnostic_filename(filename, cfg)
    for cube in cubes:
        cube.attributes.update(cfg.get('output_attributes', {}))
    io.iris_save(cubes, netcdf_path)

    # Create plot
    cell_text = np.vectorize('{:.2f}'.format)(cell_data)
    col_labels = [f"{NICE_SYMBOLS.get(l, l)} / "
                  f"{NICE_UNITS.get(col_units[i], 'unknown')}"
                  for (i, l) in enumerate(col_labels)]
    (_, axes) = plt.subplots()
    axes.axis('off')
    table = axes.table(
        cellText=cell_text,
        rowLabels=row_labels,
        colLabels=col_labels,
        loc='center',
        fontsize=8.0,
    )
    table.scale(1.7, 1.7)

    # Save plot
    plot_path = os.path.join(cfg['plot_dir'], filename + '.pdf')
    plt.savefig(plot_path, bbox_inches='tight', orientation='landscape')
    logger.info("Wrote %s", plot_path)
    plt.close()

    # Provenance
    caption = (
        'Forcing, Feedback and Equilibrium Climate Sensitivity (ECS) values. '
        'SW = short wave, LW = long wave, cs = clear sky, CRE = cloud '
        'radiative effect (similar to Andrews et al., Geophys. Res. Lett., '
        '39, 2012).')
    _write_provenance(
        netcdf_path,
        plot_path,
        caption,
        sorted([d['filename'] for d in cfg['input_data'].values()]),
        cfg,
    )


def _dict_to_array(dict_):
    """Convert (2D) dictionary to table."""
    row_labels = list(dict_.keys())

    # Columns (variables)
    all_cols = set()
    for row in dict_.values():
        all_cols |= set(row.keys())
    col_labels = list(all_cols)
    col_labels.sort()

    # Data
    cell_data = np.array(
        [[dict_.get(row, {}).get(col, np.nan) for col in col_labels]
         for row in row_labels])
    col_units = [UNITS.get(var, 'unknown') for var in col_labels]
    return (cell_data, row_labels, col_labels, col_units)


def _get_anomaly_data(input_data, year_idx=None):
    """Calculate anomaly data for all variables."""
    logger.info("Calculating anomaly data")
    project = input_data[0]['project']
    new_input_data = []
    for (var, var_data) in group_metadata(input_data, 'short_name').items():
        grouped_data = group_metadata(var_data, 'dataset')
        for (dataset_name, datasets) in grouped_data.items():
            logger.debug("Calculating '%s' anomaly for dataset '%s'", var,
                         dataset_name)
            data_4x = select_metadata(datasets, exp=EXP_4XCO2[project])
            data_pic = select_metadata(datasets, exp='piControl')

            # Check if all experiments are available
            if not data_4x:
                raise ValueError(
                    f"No '{EXP_4XCO2[project]}' data available for '{var}' of "
                    f"'{dataset_name}'")
            if not data_pic:
                raise ValueError(
                    f"No 'piControl' data available for '{var}' of "
                    f"'{dataset_name}'")

            # Calculate anomaly, extract correct years and save it
            cube = calculate_anomaly(data_4x, data_pic)
            _check_cube_dimensions(cube)
            cube = cube[year_idx]
            new_input_data.append({
                **data_4x[0],
                'ancestors': [data_4x[0]['filename'], data_pic[0]['filename']],
                'cube':
                cube,
            })
    msg = '' if not COORDS else f" with additional coordinates {COORDS['rad']}"
    logger.info("Found %iD 'tas' data and %iD radiation data%s",
                NDIMS.get('tas'), NDIMS.get('rad'), msg)
    return new_input_data


def _get_cube_list_for_table(cell_data, row_labels, col_labels, col_units):
    """Create :class:`iris.cube.CubeList` representing a table."""
    aux_coord = iris.coords.AuxCoord(row_labels, long_name='dataset')
    cubes = iris.cube.CubeList()
    for (idx, label) in enumerate(col_labels):
        if label in ('ECS', 'F', 'rtnt') and RTMT_DATASETS:
            rtmt_datasets = sorted(list(RTMT_DATASETS))
            attrs = {'net_toa_radiation': RTMT_TEXT.format(rtmt_datasets)}
        else:
            attrs = {}
        cube = iris.cube.Cube(
            np.ma.masked_invalid(cell_data[:, idx]),
            var_name=VAR_NAMES.get(label, label),
            long_name=LONG_NAMES.get(label, label),
            units=col_units[idx],
            aux_coords_and_dims=[(aux_coord, 0)],
            attributes=attrs,
        )
        cubes.append(cube)
    return cubes


def _get_data_time_last(cube):
    """Get data of :class:`iris.cube.Cube` with time axis as last dimension."""
    return np.moveaxis(cube.data, cube.coord_dims('time')[0], -1)


@partial(np.vectorize, excluded=['x_arr'], signature='(n),(n)->()')
def _get_intercept(x_arr, y_arr):
    """Get intercept of linear regression of two (masked) arrays."""
    if np.ma.is_masked(y_arr):
        x_arr = x_arr[~y_arr.mask]
        y_arr = y_arr[~y_arr.mask]
    if len(y_arr) < 2:
        return np.nan
    reg = stats.linregress(x_arr, y_arr)
    return reg.intercept


def _get_levels():
    """Get symmetric levels for contour plot.

    This function might be changed to consider cube in the future.

    """
    n_levels_per_sign = 50
    val_max = 8.0
    val_min = -8.0
    max_range = max([val_max, -val_min])
    range_ = np.linspace(0.0, max_range, n_levels_per_sign + 1)[1:]
    levels = list(-range_[::-1]) + [0.0] + list(range_)
    return levels


def _get_mmm_rad(rad_var, rad_datasets):
    """Get multi-model mean for radiation data."""
    logger.debug("Calculating multi-model mean for variable '%s'", rad_var)
    ancestors = []
    dataset_names = []
    mmm = []
    for dataset in rad_datasets:
        cube = dataset['cube']
        ancestors.extend(dataset['ancestors'])
        dataset_names.append(dataset['dataset'])
        mmm.append(cube.data)
    _check_array_shapes(mmm, rad_var)
    mmm = np.ma.array(mmm)
    mmm_cube = cube.copy(data=np.ma.mean(mmm, axis=0))
    attributes = {
        'ancestors': ancestors,
        'dataset': 'MultiModelMean',
        'datasets': '|'.join(dataset_names),
        'project': rad_datasets[0]['project'],
        'short_name': rad_var,
    }
    mmm_cube.attributes = attributes
    return {**attributes, 'cube': mmm_cube}


def _get_mmm_tas(rad_var, rad_datasets, tas_datasets):
    """Get multi-model mean for tas data."""
    logger.debug(
        "Calculating multi-model mean 'tas' for radiation variable '%s'",
        rad_var)
    ancestors = []
    dataset_names = []
    mmm = []
    for dataset_name in [d['dataset'] for d in rad_datasets]:
        tas_data = select_metadata(tas_datasets, dataset=dataset_name)
        if not tas_data:
            raise ValueError(
                f"No 'tas' data for dataset '{dataset_name}' available for "
                f"multi-model mean calculation")
        cube = tas_data[0]['cube']
        ancestors.extend(tas_data[0]['ancestors'])
        dataset_names.append(dataset_name)
        mmm.append(cube.data)
    _check_array_shapes(mmm, 'tas')
    mmm = np.ma.array(mmm)
    mmm_cube = cube.copy(data=np.ma.mean(mmm, axis=0))
    attributes = {
        'ancestors': ancestors,
        'dataset': 'MultiModelMean',
        'datasets': '|'.join(dataset_names),
        'project': rad_datasets[0]['project'],
        'short_name': _get_tas_var('MultiModelMean', rad_var),
    }
    mmm_cube.attributes = attributes
    return {**attributes, 'cube': mmm_cube}


def _get_multi_model_mean(input_data):
    """Get multi-model mean for all variables."""
    logger.info("Calculating multi-model means")
    mmm_data = []
    tas_data = select_metadata(input_data, short_name='tas')
    for (var, datasets) in group_metadata(input_data, 'short_name').items():
        if var == 'tas':
            continue
        mmm_rad = _get_mmm_rad(var, datasets)
        mmm_tas = _get_mmm_tas(var, datasets, tas_data)
        mmm_data.append(mmm_rad)
        mmm_data.append(mmm_tas)
    input_data.extend(mmm_data)
    return input_data


def _get_provenance_record(caption):
    """Create a provenance record."""
    record = {
        'caption': caption,
        'statistics': ['mean', 'diff'],
        'domains': ['global'],
        'authors': ['schlund_manuel'],
        'references': ['andrews12grl', 'gregory04grl'],
        'realms': ['atmos'],
        'themes': ['phys'],
    }
    return record


@partial(np.vectorize, excluded=['x_arr'], signature='(n),(n)->()')
def _get_slope(x_arr, y_arr):
    """Get slope of linear regression of two (masked) arrays."""
    if np.ma.is_masked(y_arr):
        x_arr = x_arr[~y_arr.mask]
        y_arr = y_arr[~y_arr.mask]
    if len(y_arr) < 2:
        return np.nan
    reg = stats.linregress(x_arr, y_arr)
    return reg.slope


def _get_tas_var(dataset_name, rad_var):
    """Get correct tas data for a certain radiation variable."""
    if dataset_name == 'MultiModelMean':
        return f'tas_{rad_var}'
    return 'tas'


def _vectorized_linregress(x_arr, y_arr):
    """Vectorized version if :func:`scipy.stats.linregress`."""
    slope = np.vectorize(lambda x, y: stats.linregress(x, y).slope,
                         signature='(n),(n)->()')
    intercept = np.vectorize(lambda x, y: stats.linregress(x, y).intercept,
                             signature='(n),(n)->()')
    return (slope(x_arr, y_arr), intercept(x_arr, y_arr))


def _write_scalar_data(data, ancestor_files, cfg, description=None):
    """Write scalar data for multiple datasets."""
    var_attrs = [
        {
            'short_name': 'ecs',
            'long_name': 'Equilibrium Climate Sensitivity (Gregory method)',
            'units': cf_units.Unit('K'),
        },
        {
            'short_name': 'lambda',
            'long_name': 'Climate Feedback Parameter',
            'units': cf_units.Unit('W m-2 K-1'),
        },
    ]
    global_attrs = {'project': list(cfg['input_data'].values())[0]['project']}
    if RTMT_DATASETS:
        rtmt_datasets = sorted(list(RTMT_DATASETS))
        global_attrs['net_toa_radiation'] = RTMT_TEXT.format(rtmt_datasets)
    for (idx, var_attr) in enumerate(var_attrs):
        caption = '{long_name} for multiple climate models'.format(**var_attr)
        if description is not None:
            filename = '{}_{}'.format(var_attr['short_name'],
                                      description.replace(' ', '_'))
            attributes = {'Description': description}
            caption += f' for {description}.'
        else:
            filename = var_attr['short_name']
            attributes = {}
            caption += '.'
        attributes.update(global_attrs)
        path = get_diagnostic_filename(filename, cfg)
        if not data[idx]:
            raise ValueError(f"Cannot write file {path}, no data for variable "
                             f"'{var_attr['short_name']}' given")

        # Scalar data
        if NDIMS['rad'] == 1:
            io.save_scalar_data({d: data[idx][d].data
                                 for d in data[idx]},
                                path,
                                var_attr,
                                attributes=attributes)

        # 1D data
        elif NDIMS['rad'] == 2:
            io.save_1d_data(data[idx],
                            path,
                            COORDS['rad'][0],
                            var_attr,
                            attributes=attributes)

        # Higher dimensions
        else:
            logger.info(
                "Writing netcdf summary file including ECS and feedback "
                "parameters for all datasets is not supported for %iD data "
                "yet", NDIMS['rad'])
            return

        # Provenance
        provenance_record = _get_provenance_record(caption)
        provenance_record['ancestors'] = ancestor_files
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(path, provenance_record)


def _write_provenance(netcdf_path, plot_path, caption, ancestors, cfg,
                      **kwargs):
    """Write provenance information for a single dataset cube."""
    provenance_record = _get_provenance_record(caption)
    provenance_record.update({
        'ancestors': ancestors,
        **kwargs,
    })
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)
        provenance_logger.log(plot_path, provenance_record)


def calculate_ecs(input_data, cfg, description=None):
    """Calculate ECS and net climate feedback parameters."""
    logger.info("Calculating ECS and net climate feedback parameter")
    msg = '' if description is None else f' for {description}'
    ancestors = []
    ecs = {}
    feedback_parameter = {}

    # Iterate over all datasets and save ECS and feedback parameters
    for dataset in select_metadata(input_data, short_name='tas'):
        dataset_name = dataset['dataset']
        logger.debug("Calculating ECS%s of dataset '%s'", msg, dataset_name)
        rtnt_data = select_metadata(input_data,
                                    short_name='rtnt',
                                    dataset=dataset_name)
        if not rtnt_data:
            logger.debug(
                "No 'rtmt' or 'rtnt' data for '%s' available, skipping ECS "
                "calculation for it", dataset_name)
            continue
        tas_cube = dataset['cube']
        rtnt_cube = rtnt_data[0]['cube']
        if rtnt_cube.ndim > 2:
            raise ValueError(
                f"Calculating ECS is only supported for cubes with less than "
                f"3 dimensions, got {rtnt_cube.ndim:d}D cube")
        ancestors.extend(dataset['ancestors'] + rtnt_data[0]['ancestors'])
        coords = [(coord, idx - 1)
                  for (idx,
                       coord) in enumerate(rtnt_cube.coords(dim_coords=True))
                  if coord.name() != 'time']

        # Calculate ECS (using linear regression)
        reg = _vectorized_linregress(_get_data_time_last(tas_cube),
                                     _get_data_time_last(rtnt_cube))
        ecs[dataset_name] = iris.cube.Cube(-reg[1] / (2 * reg[0]),
                                           dim_coords_and_dims=coords)
        feedback_parameter[dataset_name] = iris.cube.Cube(
            reg[0], dim_coords_and_dims=coords)
    ancestors = list(set(ancestors))
    if not ecs:
        logger.info(
            "No 'rtmt' or 'rtnt' data available, skipping ECS calculation")
        return

    # Write data
    _write_scalar_data([ecs, feedback_parameter], ancestors, cfg, description)


def check_input_data(cfg):
    """Check input data."""
    if not variables_available(cfg, ['tas']):
        raise ValueError("This diagnostic needs the variable 'tas'")
    input_data = cfg['input_data'].values()
    project_group = group_metadata(input_data, 'project')
    projects = list(project_group.keys())
    if len(projects) > 1:
        raise ValueError(
            f"This diagnostic supports only unique 'project' attributes, got "
            f"{projects}")
    project = projects[0]
    if project not in EXP_4XCO2:
        raise ValueError(f"Project '{project}' not supported yet")
    exp_group = group_metadata(input_data, 'exp')
    exps = set(exp_group.keys())
    if exps != {'piControl', EXP_4XCO2[project]}:
        raise ValueError(
            f"This diagnostic needs 'piControl' and '{EXP_4XCO2[project]}' "
            f"experiments, got {exps}")


def plot_feedback_parameters(input_data, cfg, description=None):
    """Plot feedback parameters vs. remaining dimension(s)."""
    # Iterate over radiation quantities (y axis)
    for (var, datasets) in group_metadata(input_data, 'short_name').items():
        if 'tas' in var:
            continue
        logger.info("Creating feedback parameter plots for variable '%s'", var)

        # Iterate over all available datasets
        for dataset in datasets:
            dataset_name = dataset['dataset']
            tas_data = select_metadata(input_data,
                                       short_name=_get_tas_var(
                                           dataset_name, var),
                                       dataset=dataset_name)
            if not tas_data:
                raise ValueError(
                    f"No 'tas' data for '{dataset_name}' available")
            cube = dataset['cube']
            tas_cube = tas_data[0]['cube']
            if cube.ndim not in (2, 3):
                raise ValueError(
                    f"Feedback plots are not supported for {cube.ndim:d}D "
                    f"input data, this requires 2D or 3D data")

            # Create plot
            (plot_path,
             feedback_cube) = _create_feedback_plot(tas_cube,
                                                    cube,
                                                    dataset_name,
                                                    cfg,
                                                    description=description)
            (netcdf_path,
             caption) = _create_feedback_file(feedback_cube,
                                              dataset_name,
                                              cfg,
                                              description=description)

            # Provenance
            if 'latitude' in COORDS['rad'] and 'longitude' in COORDS['rad']:
                plot_types = ['geo']
            elif 'latitude' in COORDS['rad']:
                plot_types = ['zonal']
            else:
                plot_types = ['other']
            _write_provenance(
                netcdf_path,
                plot_path,
                caption,
                dataset['ancestors'] + tas_data[0]['ancestors'],
                cfg,
                plot_types=plot_types,
            )


def plot_regressions(input_data, cfg, description=None):
    """Plot linear regressions used to calculate feedback parameters."""
    table = OrderedDict()

    # Iterate over radiation quantities (y axis)
    for (var, datasets) in group_metadata(input_data, 'short_name').items():
        if 'tas' in var:
            continue
        logger.info("Creating regression plots for variable '%s'", var)

        # Iterate over all available datasets
        for dataset in datasets:
            dataset_name = dataset['dataset']
            table.setdefault(dataset_name, {})
            tas_data = select_metadata(input_data,
                                       short_name=_get_tas_var(
                                           dataset_name, var),
                                       dataset=dataset_name)
            if not tas_data:
                raise ValueError(
                    f"No 'tas' data for '{dataset_name}' available")
            tas_cube = tas_data[0]['cube']
            if dataset['cube'].ndim > 1:
                raise ValueError(
                    "Regression plots are not supported for input data with "
                    "more than one dimension (which should be time)")

            # Save plot and netcdf file
            (plot_path, reg) = _create_regression_plot(tas_cube,
                                                       dataset['cube'],
                                                       dataset_name,
                                                       cfg,
                                                       description=description)
            netcdf_path = _create_regression_file(tas_cube,
                                                  dataset['cube'],
                                                  dataset_name,
                                                  cfg,
                                                  description=description)

            # Expand table
            table[dataset_name][var] = reg.slope
            if var == 'rtnt':
                table[dataset_name]['ECS'] = (-reg.intercept / 2.0 / reg.slope)
                table[dataset_name]['F'] = reg.intercept

            # Provenance
            caption = (
                'Scatterplot between {} TOA radiance and global mean surface '
                'temperature anomaly{} of the abrupt 4x CO2 experiment '
                'including linear regression for {} (following Andrews et '
                'al., Geophys. Res. Lett., 39, 2012).'.format(
                    FEEDBACK_PARAMETERS.get(var, var),
                    '' if description is None else f' for {description}',
                    dataset_name))
            _write_provenance(netcdf_path,
                              plot_path,
                              caption,
                              dataset['ancestors'] + tas_data[0]['ancestors'],
                              cfg,
                              plot_types=['scatter'])

    # Create summary table
    _create_table(table, cfg, description=description)


def preprocess_data(cfg, year_idx=None):
    """Calculate anomalies and multi-model mean."""
    input_data = deepcopy(list(cfg['input_data'].values()))
    input_data = sorted_metadata(input_data, ['short_name', 'exp', 'dataset'])

    # Use 'rtmt' instead of 'rtmt' if necessary
    for dataset in input_data:
        if dataset['short_name'] == 'rtmt':
            RTMT_DATASETS.add(dataset['dataset'])
            dataset['short_name'] = 'rtnt'
    if RTMT_DATASETS:
        logger.info("Using 'rtmt' instead of 'rtnt' for datasets '%s'",
                    RTMT_DATASETS)

    # Calculate anomalies for every dataset
    input_data = _get_anomaly_data(input_data, year_idx)

    # Calculate multi-model mean
    if cfg['calculate_mmm']:
        input_data = _get_multi_model_mean(input_data)

        # Remove other datasets if desired
        if cfg['only_consider_mmm']:
            logger.info("Removing all datasets except for 'MultiModelMean'")
            input_data = [d for d in input_data if d['dataset'] ==
                          'MultiModelMean']

    return input_data


def set_default_cfg(cfg):
    """Set default values for cfg."""
    cfg = deepcopy(cfg)
    cfg.setdefault('calculate_mmm', True)
    cfg.setdefault('only_consider_mmm', False)
    cfg.setdefault('seaborn_settings', {})
    if cfg['only_consider_mmm'] and not cfg['calculate_mmm']:
        logger.warning("Automatically setting 'calculate_mmm' to 'True' since "
                       "'only_consider_mmm' is set to 'True'")
        cfg['calculate_mmm'] = True
    return cfg


def main(cfg):
    """Run the diagnostic."""
    cfg = set_default_cfg(cfg)
    sns.set(cfg['seaborn_settings'])
    check_input_data(cfg)
    year_indices = {
        'all 150 years': slice(None),
        'first 20 years': slice(None, 20),
        'last 130 years': slice(20, None),
    }
    for (descr, year_idx) in year_indices.items():
        logger.info("Considering %s for all datasets", descr)
        input_data = preprocess_data(cfg, year_idx)

        # Calculate and save ECS
        if NDIMS['rad'] < 3:
            calculate_ecs(input_data, cfg, description=descr)
        else:
            logger.info("No ECS calculation for %iD data available",
                        NDIMS['rad'])

        # Plots
        if NDIMS['rad'] == 1:
            plot_regressions(input_data, cfg, description=descr)
        elif NDIMS['rad'] in (2, 3):
            plot_feedback_parameters(input_data, cfg, description=descr)
        else:
            logger.info("No plots for %iD data available", NDIMS['rad'])


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
