#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simple preprocessing of MLR model input.

Description
-----------
This diagnostic performs preprocessing operations for datasets used as MLR
model input in a desired way. It can also be used to process output of MLR
models for plotting.

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
aggregate_by: dict, optional
    Aggregate over given coordinates (dict values; given as  :obj:`list` of
    :obj:`str`) using a desired aggregator (dict key; given as  :obj:`str`).
    Allowed aggregators are ``'max'``, ``'mean'``, ``'median'``, ``'min'``,
    ``'sum'``, ``'std'``, ``'var'``, and ``'trend'``.
apply_common_mask: bool, optional (default: False)
    Apply common mask to all datasets. Requires identical shapes for all
    datasets.
area_weighted: bool, optional (default: True)
    Use weighted aggregation when collapsing over latitude and/or longitude
    using ``collapse``. Weights are estimated using grid cell boundaries. Only
    possible if the dataset contains ``latitude`` and ``longitude``
    coordinates.
argsort: dict, optional
    Calculate :func:`numpy.ma.argsort` along given coordinate to get ranking.
    The coordinate can be specified by the ``coord`` key. If ``descending`` is
    set to ``True``, use descending order instead of ascending.
collapse: dict, optional
    Collapse over given coordinates (dict values; given as :obj:`list` of
    :obj:`str`) using a desired aggregator (dict key; given as :obj:`str`).
    Allowed aggregators are ``'max'``, ``'mean'``, ``'median'``, ``'min'``,
    ``'sum'``, ``'std'``, ``'var'``, and ``'trend'``.
convert_units_to: str, optional
    Convert units of the input data. Can also be given as dataset option.
extract: dict, optional
    Extract certain values (dict values, given as :obj:`int`, :obj:`float` or
    iterable of them) for certain coordinates (dict keys, given as :obj:`str`).
extract_ignore_bounds: bool, optional (default: False)
    If ``True``, ignore coordinate bounds when using ``extract`` or
    ``extract_range``. If ``False``, consider coordinate bounds when using
    ``extract`` or ``extract_range``. For time coordinates, bounds are always
    ignored.
extract_range: dict, optional
    Like ``extract``, but instead of specific values extract ranges (dict
    values, given as iterable of exactly two :obj:`int` s or :obj:`float` s)
    for certain coordinates (dict keys, given as :obj:`str`).
ignore: list of dict, optional
    Ignore specific datasets by specifying multiple :obj:`dict` s of metadata.
landsea_fraction_weighted: str, optional
    When given, use land/sea fraction for weighted aggregation when collapsing
    over latitude and/or longitude using ``collapse``. Only possible if the
    dataset contains ``latitude`` and ``longitude`` coordinates. Must be one of
    ``'land'``, ``'sea'``.
mask: dict of dict
    Mask datasets. Keys have to be :mod:`numpy.ma` conversion operations (see
    `<https://docs.scipy.org/doc/numpy/reference/routines.ma.html>`_) and
    values all the keyword arguments of them.
n_jobs: int (default: 1)
    Maximum number of jobs spawned by this diagnostic script. Use ``-1`` to use
    all processors. More details are given `here
    <https://scikit-learn.org/stable/glossary.html#term-n-jobs>`_.
normalize_by_mean: bool, optional (default: False)
    Remove total mean of the dataset in the last step (resulting mean will be
    0.0). Calculates weighted mean if ``area_weighted``, ``time_weighted`` or
    ``landsea_fraction_weighted`` are set and the cube contains the
    corresponding coordinates. Does not apply to error datasets.
normalize_by_std: bool, optional (default: False)
    Scale total standard deviation of the dataset in the last step (resulting
    standard deviation will be 1.0).
output_attributes: dict, optional
    Write additional attributes to netcdf files, e.g. ``'tag'``.
pattern: str, optional
    Pattern matched against ancestor file names.
ref_calculation: str, optional
    Perform calculations involving reference dataset. Must be one of ``merge``
    (simply merge two datasets by adding the data of the reference dataset as
    :class:`iris.coords.AuxCoord` to the original dataset), ``add`` (add
    reference dataset), ``divide`` (divide by reference dataset), ``multiply``
    (multiply with reference dataset), ``subtract`` (subtract reference
    dataset) or ``trend`` (use reference dataset as x axis for calculation of
    linear trend along a specified axis, see ``ref_kwargs``).
ref_kwargs: dict, optional
    Keyword arguments for calculations involving reference datasets. Allowed
    keyword arguments are:

    * ``matched_by`` (:obj:`list` of :obj:`str`, default: ``[]``): Use a
      given set of attributes to match datasets with their corresponding
      reference datasets (specified by ``ref = True``).
    * ``collapse_over`` (:obj:`str`, default: ``'time'``): Coordinate which
      is collapsed. Only relevant when ``ref_calculation`` is set to ``trend``.
return_trend_stderr: bool, optional (default: True)
    Return standard error of slope in case of trend calculations (as
    ``var_type`` ``prediction_input_error``).
scalar_operations: dict, optional
    Operations involving scalars. Allowed keys are ``add``, ``divide``,
    ``multiply`` or ``subtract``. The corresponding values (:obj:`float` or
    :obj:`int`) are scalars that are used with the operations.
time_weighted: bool, optional (default: True)
    Use weighted aggregation when collapsing over time dimension using
    ``collapse``. Weights are estimated using grid cell boundaries.
unify_coords_to: dict, optional
    If given, replace coordinates of all datasets with that of a reference cube
    (if necessary and possible, broadcast beforehand). The reference dataset
    is determined by keyword arguments given to this option (keyword arguments
    must point to exactly one dataset).

"""

import datetime
import functools
import logging
import os
import warnings
from copy import deepcopy

import dask.array as da
import iris
import numpy as np
from cf_units import Unit
from joblib import Parallel, delayed
from scipy import stats

from esmvaltool.diag_scripts import mlr
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    io,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))

AGGREGATORS = {
    'max': iris.analysis.MAX,
    'mean': iris.analysis.MEAN,
    'median': iris.analysis.MEDIAN,
    'min': iris.analysis.MIN,
    'std': iris.analysis.STD_DEV,
    'sum': iris.analysis.SUM,
    'var': iris.analysis.VARIANCE,
}


def _add_categorized_time_coords(cube, coords, aggregator):
    """Add categorized time coordinates to cube."""
    for coord_name in coords:
        if cube.coords(coord_name):
            continue
        if hasattr(iris.coord_categorisation, f'add_{coord_name}'):
            getattr(iris.coord_categorisation, f'add_{coord_name}')(cube,
                                                                    'time')
            logger.debug("Added coordinate '%s' to cube", coord_name)
        else:
            raise ValueError(
                f"Cannot aggregate over coordinate(s) '{coords}' using "
                f"'{aggregator}': Categorized coordinate '{coord_name}' is "
                f"not a coordinate of cube {cube.summary(shorten=True)} and "
                f"cannot be added via iris.coord_categorisation")


def _apply_trend_aggregator(cfg, cube, data, coord_name):
    """Apply aggregator ``trend`` to cube."""
    return_stderr = _return_stderr(cfg, data)
    units = cube.units

    # Get corresponding dimensional coordinate
    coord_dims = cube.coord_dims(coord_name)
    if len(coord_dims) != 1:
        raise ValueError(
            f"Trend aggregation along coordinate '{coord_name}' requires 1D "
            f"coordinate, got {len(coord_dims):d}D coordinate")
    dim_coord = cube.coord(dim_coords=True, dimensions=coord_dims[0])

    # Calculate trends in parallel
    parallel = Parallel(n_jobs=cfg['n_jobs'])
    coord_values = np.unique(cube.coord(coord_name).points)
    cube_slices = [cube.extract(iris.Constraint(**{coord_name: val})) for
                   val in coord_values]
    all_cubes = parallel(
        [delayed(_calculate_slope_along_coord)(cube_slice, dim_coord.name(),
                                               return_stderr=return_stderr)
         for cube_slice in cube_slices]
    )

    # Merge output (Original units might get lost in pool)
    cubes = [tup[0] for tup in all_cubes]
    cube = iris.cube.CubeList(cubes).merge_cube()
    cube.units = units
    if return_stderr:
        cube_stderr = iris.cube.CubeList(
            [tup[1] for tup in all_cubes]).merge_cube()
        cube_stderr.units = units
    else:
        cube_stderr = None
    units = _get_coord_units(cube, coord_name)
    (cube, data) = _set_trend_metadata(cfg, cube, cube_stderr, data, units)
    data['trend'] = f'aggregated along coordinate {coord_name}'
    return (cube, data)


def _calculate_slope_along_coord(cube, coord_name, return_stderr=True):
    """Calculate slope of a cube along a given coordinate."""
    coord = cube.coord(coord_name)
    coord_dims = cube.coord_dims(coord_name)
    if len(coord_dims) != 1:
        raise ValueError(
            f"Trend calculation along coordinate '{coord_name}' requires "
            f"1D coordinate, got {len(coord_dims):d}D coordinate")

    # Get slope and error if desired
    x_data = coord.points
    y_data = np.moveaxis(cube.data, coord_dims[0], -1)
    calc_slope = np.vectorize(_get_slope, excluded=['x_arr'],
                              signature='(n),(n)->()')
    slope = calc_slope(x_data, y_data)
    if return_stderr:
        calc_slope_stderr = np.vectorize(_get_slope_stderr, excluded=['x_arr'],
                                         signature='(n),(n)->()')
        slope_stderr = calc_slope_stderr(x_data, y_data)
    else:
        slope_stderr = None

    # Apply dummy aggregator for correct cell method and set data
    aggregator = iris.analysis.Aggregator('trend', _remove_axis)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message='Collapsing a non-contiguous coordinate',
            category=UserWarning,
            module='iris',
        )
        cube = cube.collapsed(coord_name, aggregator)
    cube.data = np.ma.masked_invalid(slope)
    if slope_stderr is not None:
        cube_stderr = cube.copy()
        cube_stderr.data = np.ma.masked_invalid(slope_stderr)
    else:
        cube_stderr = None
    return (cube, cube_stderr)


def _check_cubes(cube, ref_cube, ref_option):
    """Check cube and reference cube."""
    if cube.shape != ref_cube.shape:
        raise ValueError(
            f"Expected identical shapes for data and reference data, got "
            f"{cube.shape} and {ref_cube.shape}")
    if ref_option == 'subtract' and cube.units != ref_cube.units:
        logger.warning(
            "Got different units for original dataset and the corresponding "
            "reference dataset ('%s' and '%s') for operation '%s'",
            cube.units, ref_cube.units, ref_option)


def _collapse_over(cfg, cube, data, coords, aggregator):
    """Collapse over cube."""
    coords = deepcopy(coords)
    iris_op = AGGREGATORS[aggregator]
    if aggregator not in ('mean', 'sum'):
        cube = cube.collapsed(coords, iris_op)
        return (cube, data)

    # Latitude and/or longitude (weighted if desired)
    horizontal_coords = _get_horizontal_coordinates(coords)
    if horizontal_coords:
        horizontal_weights = _get_horizontal_weights(cfg, cube)
        cube = cube.collapsed(horizontal_coords, iris_op,
                              weights=horizontal_weights)
        for coord in horizontal_coords:
            coords.remove(coord)
        if aggregator == 'sum' and cfg['area_weighted']:
            cube.units *= Unit('m2')
            data['units'] = str(cube.units)

    # Time (weighted if desired)
    if 'time' in coords:
        time_weights = _get_time_weights(cfg, cube)
        time_units = mlr.get_absolute_time_units(cube.coord('time').units)
        cube = cube.collapsed(['time'], iris_op, weights=time_weights)
        coords.remove('time')
        if aggregator == 'sum' and time_weights is not None:
            cube.units *= time_units
            data['units'] = str(cube.units)

    # Remaining operations
    if coords:
        cube = cube.collapsed(coords, iris_op)

    return (cube, data)


def _coord_constraint(cell, value, coord_name, ignore_bounds=False,
                      interpret_as_range=False):
    """Callable that can be used to form a :class:`iris.Constraint`."""
    if coord_name == 'time' or ignore_bounds:
        cell_object = cell.point
    else:
        cell_object = cell
    if interpret_as_range:
        return value[0] <= cell_object <= value[1]
    try:
        return cell_object in value
    except TypeError:
        return cell_object == value


def _fail_if_stderr(data, description):
    """Raise exception of data is a standard error."""
    if 'stderr' in data:
        raise ValueError(
            f"{description} is not supported with standard errors yet")


def _get_all_weights(cfg, cube):
    """Get all desired weights for a cube."""
    weights = mlr.get_all_weights(
        cube, area_weighted=cfg['area_weighted'],
        time_weighted=cfg['time_weighted'],
        landsea_fraction_weighted=cfg.get('landsea_fraction_weighted'))
    return weights


def _get_constrained_cube(cube, constraints):
    """Merge multiple :class:`iris.Constraint` s and apply them to cube."""
    constraint = constraints[0]
    for new_constraint in constraints[1:]:
        constraint &= new_constraint
    return cube.extract(constraint)


def _get_coord_units(cube, coord_name):
    """Get units of cube's coordinate."""
    coord = cube.coord(coord_name)
    if coord_name == 'time':
        units = mlr.get_absolute_time_units(coord.units)
    else:
        units = coord.units
    return units


def _get_error_datasets(input_data, **kwargs):
    """Extract error datasets from input data."""
    input_data = select_metadata(input_data, **kwargs)
    error_data = []
    for dataset in input_data:
        if dataset.get('stderr', False):
            error_data.append(dataset)
    return error_data


def _get_horizontal_coordinates(coords):
    """Extract horizontal coordinates from :obj:`list` of coordinates."""
    horizontal_coords = []
    if 'latitude' in coords:
        horizontal_coords.append('latitude')
    if 'longitude' in coords:
        horizontal_coords.append('longitude')
    return horizontal_coords


def _get_horizontal_weights(cfg, cube):
    """Get weights for horizontal dimensions."""
    weights = mlr.get_horizontal_weights(
        cube,
        area_weighted=cfg['area_weighted'],
        landsea_fraction_weighted=cfg.get('landsea_fraction_weighted'))
    return weights


def _get_landsea_fraction_weights(cfg, cube):
    """Calculate land/sea fraction weights."""
    landsea_fraction_weights = None
    if 'landsea_fraction_weighted' in cfg:
        landsea_fraction_weights = mlr.get_landsea_fraction_weights(
            cube, cfg['landsea_fraction_weighted'])
    return landsea_fraction_weights


def _get_ref_calc(cfg, dataset, ref_datasets, ref_option):
    """Perform calculations involving reference datasets for regular data."""
    ref_kwargs = cfg.get('ref_kwargs', {})
    ref_dataset = _get_ref_dataset(dataset, ref_datasets, **ref_kwargs)
    cube = dataset['cube']
    ref_cube = ref_dataset['cube']
    _check_cubes(cube, ref_cube, ref_option)
    dataset['original_cube'] = cube.copy()
    dataset['ref_cube'] = ref_cube.copy()
    if ref_option == 'merge':
        aux_coord = cube_to_aux_coord(ref_cube)
        cube.add_aux_coord(aux_coord, np.arange(cube.ndim))
        suffix = None
    elif ref_option == 'add':
        cube.data += ref_cube.data
        suffix = 'plus ref'
    elif ref_option == 'multiply':
        cube.data *= ref_cube.data
        cube.units *= ref_cube.units
        dataset['units'] = str(cube.units)
        suffix = 'multiplied by ref'
    elif ref_option == 'divide':
        cube.data /= ref_cube.data
        cube.units /= ref_cube.units
        dataset['units'] = str(cube.units)
        suffix = 'divided by ref'
    elif ref_option == 'subtract':
        cube.data -= ref_cube.data
        suffix = 'minus ref'
    elif ref_option == 'trend':
        (cube, cube_stderr) = _get_trend_relative_to_ref(
            cfg, dataset, ref_cube,
            collapse_over=ref_kwargs.get('collapse_over'))
        (cube, dataset) = _set_trend_metadata(cfg, cube, cube_stderr, dataset,
                                              ref_cube.units)
        suffix = 'relative to ref'
    else:
        raise ValueError(f"Got invalid ref option '{ref_option}'")
    if suffix is not None:
        suffix_no_space = suffix.replace(' ', '_')
        dataset['standard_name'] = None
        dataset['short_name'] += f'_{suffix_no_space}'
        dataset['long_name'] += f' ({suffix})'
        exp = ('' if ref_dataset.get('exp') is None else
               f", experiment {ref_dataset['exp']}")
        dataset['reference_data'] = (
            f"{ref_dataset['short_name']} of {ref_dataset['dataset']} "
            f"(project {ref_dataset['project']}{exp}) for years "
            f"{ref_dataset['start_year']} to {ref_dataset['end_year']}")
    dataset['cube'] = cube
    return dataset


def _get_ref_calc_stderr(cfg, dataset, ref_datasets, regular_datasets,
                         ref_option):
    """Perform calculations involving reference datasets for error data."""
    ref_kwargs = cfg.get('ref_kwargs', {})

    # Extract reference dataset (error)
    ref_dataset = _get_ref_dataset(dataset, ref_datasets, **ref_kwargs)

    # Extract regular dataset (corresponding mean to error)
    excluded_keys = ['var_type', 'short_name', 'standard_name', 'long_name',
                     'variable_group', 'diagnostic', 'filename', 'cube',
                     'recipe_dataset_index', 'stderr', 'alias', 'units']
    kwargs = {key: dataset[key] for key in dataset if key not in excluded_keys}
    reg_dataset = select_metadata(regular_datasets, **kwargs)
    if len(reg_dataset) != 1:
        raise ValueError(
            f"Expected exactly one regular dataset for error dataset "
            f"{dataset}, got {len(reg_dataset):d}")
    reg_dataset = reg_dataset[0]

    # Perform calculations
    cube = dataset['cube']
    ref_cube = ref_dataset['cube']
    reg_cube = reg_dataset['cube']
    _check_cubes(cube, ref_cube, ref_option)
    if ref_option == 'merge':
        aux_coord = cube_to_aux_coord(ref_cube)
        cube.add_aux_coord(aux_coord, np.arange(cube.ndim))
    if ref_option == 'divide':
        error = np.ma.abs(reg_cube.data) * np.ma.sqrt(
            (cube.data / reg_dataset['original_cube'].data)**2 +
            (ref_cube.data / reg_dataset['ref_cube'].data)**2)
        cube.data = error
    elif ref_option == 'subtract':
        cube.data = np.ma.sqrt(cube.data**2 + ref_cube.data**2)
    elif ref_option == 'trend':
        raise ValueError(
            "Calculations involving reference datasets with option 'trend' "
            "is not supported for error datasets yet; errors are calculated "
            "from the original dataset using the standard error of slopes")
    else:
        raise NotImplementedError(
            f"Calculations involving reference datasets with option "
            f"'{ref_option}' are not supported yet")
    cube.units = reg_cube.units
    dataset['standard_name'] = reg_dataset['standard_name']
    dataset['short_name'] = reg_dataset['short_name']
    dataset['long_name'] = reg_dataset['long_name']
    dataset['units'] = reg_dataset['units']
    dataset['reference_data'] = reg_dataset['reference_data']
    dataset['cube'] = cube
    return dataset


def _get_ref_dataset(dataset, ref_datasets, **ref_kwargs):
    """Extract reference dataset for a given dataset."""
    metadata = ref_kwargs.get('matched_by', [])
    kwargs = {m: dataset[m] for m in metadata if m in dataset}
    ref_dataset = select_metadata(ref_datasets, **kwargs)
    if len(ref_dataset) != 1:
        raise ValueError(
            f"Expected exactly one reference dataset (with attribute ref "
            f"== True) for dataset {dataset}, got {len(ref_dataset):d}. "
            f"Consider extending list of metadata for option 'matched_by' in "
            f"'ref_kwargs' (used {kwargs})")
    ref_dataset = ref_dataset[0]
    return ref_dataset


def _get_single_constraint(cube, coord_name, val,
                           ignore_bounds=False, interpret_as_range=False):
    """Get single :class:`iris.Constraint`."""
    if coord_name == 'time':
        time_units = cube.coord('time').units
        val = time_units.num2date(val)
    if interpret_as_range:
        try:
            len_range = len(val)
        except TypeError:
            raise TypeError(
                f"Expected iterable for values of 'extract_range' for "
                f"coordinate '{coord_name}', got '{val}'")
        if len_range != 2:
            raise ValueError(
                f"Expected exactly two elements for range of '{coord_name}' "
                f"in 'extract_range', got {len_range:d} ({val})")
        logger.debug("Extracting range %s for coordinate '%s'", val,
                     coord_name)
    coord_vals = functools.partial(_coord_constraint, value=val,
                                   coord_name=coord_name,
                                   ignore_bounds=ignore_bounds,
                                   interpret_as_range=interpret_as_range)
    return iris.Constraint(**{coord_name: coord_vals})


def _get_slope(x_arr, y_arr):
    """Get slope of linear regression of two (masked) arrays."""
    if np.ma.is_masked(y_arr):
        x_arr = x_arr[~y_arr.mask]
        y_arr = y_arr[~y_arr.mask]
    if len(y_arr) < 2:
        return np.nan
    reg = stats.linregress(x_arr, y_arr)
    return reg.slope


def _get_slope_stderr(x_arr, y_arr):
    """Get standard error of linear slope of two (masked) arrays."""
    if np.ma.is_masked(y_arr):
        x_arr = x_arr[~y_arr.mask]
        y_arr = y_arr[~y_arr.mask]
    if len(y_arr) < 2:
        return np.nan
    reg = stats.linregress(x_arr, y_arr)
    return reg.stderr


def _get_time_weights(cfg, cube):
    """Calculate time weights."""
    time_weights = None
    if cfg['time_weighted']:
        time_weights = mlr.get_time_weights(cube)
    return time_weights


def _get_trend_relative_to_ref(cfg, data, ref_cube, collapse_over=None):
    """Calculate linear trend relative to reference dataset."""
    if collapse_over is None:
        collapse_over = 'time'
    cube = data['cube']
    return_stderr = _return_stderr(cfg, data)

    # Get coordinate
    coord_dims = cube.coord_dims(collapse_over)
    if len(coord_dims) != 1:
        raise ValueError(
            f"Trend calculation involving reference dataset along coordinate "
            f"'{collapse_over}' requires 1D coordinate, got "
            f"{len(coord_dims):d}D coordinate")
    if ref_cube.coord_dims(collapse_over) != coord_dims:
        raise ValueError(
            f"Trend calculation involving reference dataset along coordinate "
            f"'{collapse_over}' requires that the coordinate covers identical "
            f"dimensions for the dataset and reference dataset, got "
            f"{coord_dims} and {ref_cube.coord_dims(collapse_over)}")

    # Get slope and error if desired
    x_data = np.moveaxis(ref_cube.data, coord_dims[0], -1)
    y_data = np.moveaxis(cube.data, coord_dims[0], -1)
    calc_slope = np.vectorize(_get_slope, signature='(n),(n)->()')
    slope = calc_slope(x_data, y_data)
    if return_stderr:
        calc_slope_stderr = np.vectorize(_get_slope_stderr,
                                         signature='(n),(n)->()')
        slope_stderr = calc_slope_stderr(x_data, y_data)
    else:
        slope_stderr = None

    # Apply dummy aggregator for correct cell method and set data
    aggregator = iris.analysis.Aggregator('trend using ref', _remove_axis)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message='Collapsing a non-contiguous coordinate',
            category=UserWarning,
            module='iris',
        )
        cube = cube.collapsed(collapse_over, aggregator)
    cube.data = np.ma.masked_invalid(slope)
    if slope_stderr is not None:
        cube_stderr = cube.copy()
        cube_stderr.data = np.ma.masked_invalid(slope_stderr)
    else:
        cube_stderr = None
    return (cube, cube_stderr)


def _remove_axis(data, axis=None):
    """Remove given axis of arrays by the first index of a given axis."""
    return np.take(data, 0, axis=axis)


def _return_stderr(cfg, data):
    """Check if standard error should be returned."""
    return (data.get('var_type') == 'prediction_input' and
            cfg['return_trend_stderr'])


def _promote_aux_coord(cube, data, coord_name):
    """Promote auxiliary coordinate to dimensional coordinate."""
    aux_coords = [coord.name() for coord in cube.coords(dim_coords=False)]
    if coord_name in aux_coords:
        try:
            iris.util.promote_aux_coord_to_dim_coord(cube, coord_name)
        except ValueError as exc:
            logger.debug(
                "Could not promote coordinate '%s' to dimensional "
                "coordinate: %s", coord_name, str(exc))
        else:
            if isinstance(data.get('stderr'), dict):
                stderr_cube = data['stderr']['cube']
                iris.util.promote_aux_coord_to_dim_coord(stderr_cube,
                                                         coord_name)


def _set_trend_metadata(cfg, cube, cube_stderr, data, units):
    """Set correct metadata for trend calculation."""
    cube.units /= units
    data['standard_name'] = None
    data['short_name'] += '_trend'
    data['long_name'] += ' (trend)'
    data['units'] += str(cube.units)
    if cube_stderr is not None:
        cube_stderr.units /= units
        stderr_data = deepcopy(data)
        stderr_data = cache_cube(cfg, cube_stderr, stderr_data)
        data['stderr'] = stderr_data
    return (cube, data)


def add_standard_errors(input_data):
    """Add calculated standard errors to list of data."""
    new_input_data = []
    for data in input_data:
        if isinstance(data.get('stderr'), dict):
            stderr_data = data.pop('stderr')
            stderr_data['stderr'] = True
            stderr_data['standard_name'] = None
            stderr_data['short_name'] += '_stderr'
            stderr_data['long_name'] += ' (Standard Error)'
            stderr_data['var_type'] += '_error'
            new_input_data.append(stderr_data)
            logger.info("Added standard error for %s", data['filename'])
        new_input_data.append(data)
    return new_input_data


def aggregate_by(cfg, cube, data):
    """Aggregate cube over specified coordinate."""
    for (aggregator, coords) in cfg.get('aggregate_by', {}).items():
        if not isinstance(coords, list):
            coords = [coords]
        if aggregator not in AGGREGATORS:
            raise ValueError(
                f"Expected one of {list(AGGREGATORS.keys())} as aggregator "
                f"for 'aggregate_by', got '{aggregator}'")
        iris_op = AGGREGATORS[aggregator]
        logger.debug("Aggregating over coordinate(s) %s by calculating %s",
                     coords, aggregator)
        _add_categorized_time_coords(cube, coords, aggregator)
        cube = cube.aggregated_by(coords, iris_op)
        if len(coords) == 1:
            _promote_aux_coord(cube, data, coords[0])
    return (cube, data)


def aggregate_by_trend(cfg, cube, data):
    """Aggregate cube over specified coordinate using ``trend``."""
    if 'trend' not in cfg.get('aggregate_by', {}):
        return (cube, data)
    coords = cfg['aggregate_by']['trend']
    if not isinstance(coords, list):
        coords = [coords]
    logger.debug("Aggregating over coordinate(s) %s by calculating 'trend'",
                 coords)
    if len(coords) != 1:
        raise ValueError(
            f"Aggregation using 'trend' is currently only supported with a "
            f"single coordinate, got {coords}")
    _add_categorized_time_coords(cube, coords, 'trend')
    coord_name = coords[0]
    (cube, data) = _apply_trend_aggregator(cfg, cube, data, coord_name)
    _promote_aux_coord(cube, data, coord_name)
    return (cube, data)


def apply_common_mask(cfg, input_data):
    """Apply common mask to all datasets."""
    if not cfg.get('apply_common_mask'):
        return input_data
    logger.info("Applying common mask to all cubes")
    shapes = {data['cube'].shape for data in input_data}
    if len(shapes) > 1:
        raise ValueError(
            f"Expected cubes with identical shapes when 'apply_common_mask' "
            f"is set to 'True', got shapes {shapes}")
    common_mask = da.full(list(shapes)[0], False)
    for data in input_data:
        common_mask |= da.ma.getmaskarray(data['cube'].core_data())
    for data in input_data:
        data['cube'].data = da.ma.masked_array(data['cube'].core_data(),
                                               mask=common_mask)
    return input_data


def argsort(cfg, cube, data):
    """Calculate :func:`numpy.ma.argsort` along given axis (= Ranking)."""
    if not cfg.get('argsort'):
        return (cube, data)
    _fail_if_stderr(data, "'argsort'")
    coord = cfg['argsort'].get('coord')
    if not coord:
        raise ValueError(
            "When 'argsort' is given, a valid 'coord' needs to specified as "
            "key")
    logger.debug("Calculating argsort along coordinate '%s' to get ranking",
                 coord)
    axis = cube.coord_dims(coord)[0]
    original_mask = np.ma.getmaskarray(cube.data)
    if cfg['argsort'].get('descending'):
        ranking = np.ma.argsort(-cube.data, axis=axis, fill_value=-np.inf)
        cube.attributes['order'] = 'descending'
    else:
        ranking = np.ma.argsort(cube.data, axis=axis, fill_value=np.inf)
        cube.attributes['order'] = 'ascending'
    cube.data = np.ma.array(ranking, mask=original_mask, dtype=cube.dtype)
    cube.units = Unit('no unit')
    data['standard_name'] = None
    data['short_name'] += '_ranking'
    data['long_name'] += ' (ranking)'
    data['units'] = str(cube.units)
    return (cube, data)


def cache_cube(cfg, cube, data):
    """Cache cube in :obj:`dict`."""
    path = data['filename']
    basename = os.path.splitext(os.path.basename(path))[0]
    if cube.var_name is not None:
        basename = basename.replace(cube.var_name, data['short_name'])
        cube.var_name = data['short_name']
    if 'var_type' in data:
        for var_type in mlr.VAR_TYPES:
            if basename.endswith(f'_{var_type}'):
                basename = basename.replace(f'_{var_type}', '')
        basename += f"_{data['var_type']}"
    new_path = get_diagnostic_filename(basename, cfg)
    data['filename'] = new_path
    data['cube'] = cube
    new_attrs = cfg.get('output_attributes', {})
    data.update(new_attrs)
    return data


def collapse(cfg, cube, data):
    """Collapse data over specified coordinates."""
    for (aggregator, coords) in cfg.get('collapse', {}).items():
        if not isinstance(coords, list):
            coords = [coords]
        if aggregator not in AGGREGATORS:
            raise ValueError(
                f"Expected one of {list(AGGREGATORS.keys())} as aggregator "
                f"for 'collapse', got '{aggregator}'")
        logger.debug("Collapsing coordinate(s) %s by calculating %s", coords,
                     aggregator)
        if coords == ['all']:
            coords = [coord.name() for coord in cube.coords(dim_coords=True)]
        (cube, data) = _collapse_over(cfg, cube, data, coords, aggregator)
    return (cube, data)


def collapse_with_trend(cfg, cube, data):
    """Collapse data over specified coordinates using ``trend``."""
    if 'trend' not in cfg.get('collapse', {}):
        return (cube, data)
    coords = cfg['collapse']['trend']
    if not isinstance(coords, list):
        coords = [coords]
    logger.debug("Collapsing coordinate(s) %s by calculating 'trend'", coords)
    if coords == ['all']:
        coords = [coord.name() for coord in cube.coords(dim_coords=True)]
    if len(coords) != 1:
        raise ValueError(
            f"Collapsing using 'trend' is currently only supported with a "
            f"single coordinate, got {coords}")
    coord_name = coords[0]
    if not cube.coords(coord_name):
        raise iris.exceptions.CoordinateNotFoundError(
            f"Cannot calculate trend along '{coord_name}', cube "
            f"{cube.summary(shorten=True)} does not contain a coordinate "
            f"with that name")
    return_stderr = _return_stderr(cfg, data)
    (cube,
     cube_stderr) = _calculate_slope_along_coord(
         cube, coord_name, return_stderr=return_stderr)
    units = _get_coord_units(cube, coord_name)
    (cube, data) = _set_trend_metadata(cfg, cube, cube_stderr, data, units)
    data['trend'] = f'along coordinate {coord_name}'
    return (cube, data)


def convert_units_to(cfg, cube, data):
    """Convert units if desired."""
    cfg_settings = cfg.get('convert_units_to')
    data_settings = data.get('convert_units_to')
    if cfg_settings or data_settings:
        units_to = cfg_settings
        if data_settings:
            units_to = data_settings
        logger.debug("Converting units from '%s' to '%s'", cube.units,
                     units_to)
        try:
            cube.convert_units(units_to)
        except ValueError:
            raise ValueError(
                f"Cannot convert units of cube {cube.summary(shorten=True)} "
                f"from '{cube.units}' to '{units_to}'")
        data['units'] = str(cube.units)
    return (cube, data)


def cube_to_aux_coord(cube):
    """Convert :class:`iris.cube.Cube` to :class:`iris.coords.AuxCoord`."""
    aux_coord = iris.coords.AuxCoord(cube.data,
                                     var_name=cube.var_name,
                                     standard_name=cube.standard_name,
                                     long_name=cube.long_name,
                                     units=cube.units)
    return aux_coord


def extract(cfg, cube):
    """Extract specific coordinate values."""
    if not cfg.get('extract'):
        return cube
    constraints = []
    for (coord_name, val) in cfg['extract'].items():
        constraint = _get_single_constraint(
            cube, coord_name, val, ignore_bounds=cfg['extract_ignore_bounds'])
        constraints.append(constraint)
    new_cube = _get_constrained_cube(cube, constraints)
    if new_cube is None:
        raise ValueError(
            f"Extracting {cfg['extract']} from cube "
            f"{cube.summary(shorten=True)} yielded empty cube")
    return new_cube


def extract_range(cfg, cube):
    """Extract range of coordinate values."""
    if not cfg.get('extract_range'):
        return cube
    constraints = []
    for (coord_name, coord_range) in cfg['extract_range'].items():
        constraint = _get_single_constraint(
            cube, coord_name, coord_range,
            ignore_bounds=cfg['extract_ignore_bounds'],
            interpret_as_range=True)
        constraints.append(constraint)
    new_cube = _get_constrained_cube(cube, constraints)
    if new_cube is None:
        raise ValueError(
            f"Extracting range {cfg['extract_range']} from cube "
            f"{cube.summary(shorten=True)} yielded empty cube")
    return new_cube


def get_ref_cube(input_data, **kwargs):
    """Extract reference dataset."""
    logger.info("Using keyword arguments %s to extract reference datasets for "
                "unifying coordinates", kwargs)
    datasets = select_metadata(input_data, **kwargs)
    if len(datasets) != 1:
        raise ValueError(
            f"Expected exactly one reference dataset for unifying coords "
            f"matching {kwargs}, got {len(datasets):d}")
    ref_cube = iris.load_cube(datasets[0]['filename'])
    return ref_cube


def load_cubes(input_data):
    """Load cubes into :obj:`dict`."""
    for data in input_data:
        path = data['filename']
        logger.info("Loading %s", path)
        cube = iris.load_cube(path)
        data['cube'] = cube
        data['original_filename'] = path
    return input_data


def mask(cfg, cube):
    """Perform masking operations."""
    n_masked_values_old = np.count_nonzero(np.ma.getmaskarray(cube.data))
    for (masking_op, kwargs) in cfg.get('mask', {}).items():
        if not hasattr(np.ma, masking_op):
            raise AttributeError(
                f"Invalid masking operation, '{masking_op}' is not a function "
                f"of module numpy.ma")
        logger.debug("Applying mask operation '%s' using arguments %s",
                     masking_op, kwargs)
        masked_data = getattr(np.ma, masking_op)(cube.data, **kwargs)
        cube = cube.copy(masked_data)
    n_masked_values_new = np.count_nonzero(np.ma.getmaskarray(cube.data))
    n_total = cube.data.size
    diff = n_masked_values_new - n_masked_values_old
    if diff:
        logger.info(
            "Additionally masked %i values by operations %s (before: %i "
            "non-masked values, after: %i non-masked values", diff,
            cfg['mask'], n_total - n_masked_values_old,
            n_total - n_masked_values_new)
    return cube


def normalize_by_mean(cfg, cube, data):
    """Normalize final dataset by mean."""
    if cfg.get('normalize_by_mean') and '_error' not in data['var_type']:
        units = cube.units
        logger.debug("Normalizing mean")
        weights = _get_all_weights(cfg, cube)
        mean = np.ma.average(cube.data, weights=weights)
        cube.data -= mean
        data['long_name'] += ' (mean normalized)'
        data['normalize_by_mean'] = (
            f"Mean normalized to 0.0 {units} by subtraction, original mean "
            f"was {mean} {units}")
        data['original_mean'] = mean
    return (cube, data)


def normalize_by_std(cfg, cube, data):
    """Normalize final dataset by standard deviation."""
    if not cfg.get('normalize_by_std'):
        return (cube, data)
    units = cube.units
    logger.debug("Normalizing by standard_deviation")
    std = np.ma.std(cube.data)
    cube.data /= std
    cube.units = '1'
    data['long_name'] += ' (std normalized)'
    data['units'] = str(cube.units)
    data['normalize_by_std'] = (
        f"Standard deviation scaled to 1.0 by division, original std was "
        f"{std} {units}")
    data['original_units'] = str(units)
    data['original_std'] = std
    return (cube, data)


def ref_calculation(cfg, input_data):
    """Perform all calculation involving reference datasets."""
    if not cfg.get('ref_calculation'):
        return input_data
    ref_option = cfg['ref_calculation']
    ref_options = ['merge', 'add', 'divide', 'multiply', 'subtract', 'trend']
    if ref_option not in ref_options:
        raise ValueError(
            f"Expected one of {ref_options} for 'ref_calculation', got "
            f"'{ref_option}'")
    ref_kwargs = cfg.get('ref_kwargs', {})
    metadata = ref_kwargs.get('matched_by', [])
    logger.info("Performing calculation '%s' involving reference datasets",
                ref_option)
    logger.info("Retrieving reference dataset attributes %s to match datasets",
                metadata)
    ref_datasets = select_metadata(input_data, ref=True)
    regular_datasets_errors = _get_error_datasets(input_data, ref=False)
    regular_datasets = []
    for dataset in select_metadata(input_data, ref=False):
        if dataset not in regular_datasets_errors:
            regular_datasets.append(dataset)
    new_data = []
    logger.info(
        "Performing calculations involving reference datasets for %i regular "
        "dataset(s)", len(regular_datasets))
    for dataset in regular_datasets:
        dataset = _get_ref_calc(cfg, dataset, ref_datasets, ref_option)
        new_data.append(dataset)
    logger.info(
        "Performing calculations involving reference datasets for %i error "
        "dataset(s)", len(regular_datasets_errors))
    for dataset in regular_datasets_errors:
        dataset = _get_ref_calc_stderr(cfg, dataset, ref_datasets, new_data,
                                       ref_option)
        new_data.append(dataset)
    return new_data


def scalar_operations(cfg, cube):
    """Perform scalar operations."""
    allowed_operations = ('add', 'divide', 'multiply', 'divide')
    for (operation, constant) in cfg.get('scalar_operations', {}).items():
        if operation == 'add':
            cube.data += constant
            logger.debug("Added %f to data", constant)
        elif operation == 'divide':
            cube.data /= constant
            logger.debug("Divided %f from data", constant)
        elif operation == 'multiply':
            cube.data *= constant
            logger.debug("Multiplied %f to data", constant)
        elif operation == 'subtract':
            cube.data -= constant
            logger.debug("Subtracted %f from data", constant)
        else:
            raise ValueError(
                f"Expected one of {allowed_operations} for operation in "
                f"'scalar_operations', got '{operation}'")
    return cube


def unify_coords_to(cube, ref_cube):
    """Unify coordinates."""
    if ref_cube is None:
        return cube

    # Broadcast if necessary/possible
    if cube.shape != ref_cube.shape:
        logger.info(
            "Broadcasting %s to shape of reference cube %s",
            cube.summary(shorten=True), ref_cube.shape)
        old_cube = cube.copy()
        broadcasted_data = np.broadcast_to(
            np.ma.array(old_cube.data).filled(np.nan), ref_cube.shape)
        cube = iris.cube.Cube(np.ma.masked_invalid(broadcasted_data))
        cube.metadata = old_cube.metadata

    # Set new coordinates
    for coord in cube.coords():
        cube.remove_coord(coord)
    for dim_coord in ref_cube.coords(dim_coords=True):
        coord_dims = ref_cube.coord_dims(dim_coord)
        cube.add_dim_coord(dim_coord, coord_dims)
    for aux_coord in ref_cube.coords(dim_coords=False):
        coord_dims = ref_cube.coord_dims(aux_coord)
        cube.add_aux_coord(aux_coord, coord_dims)
    return cube


def write_cube(cfg, cube, data):
    """Write cube (check for MLR attributes and existing files first)."""
    if not mlr.datasets_have_mlr_attributes([data], log_level='error'):
        raise ValueError(
            f"Cannot write cube {cube.summary(shorten=True)} using metadata "
            f"{data}")

    # Get new path
    new_path = data['filename']
    if os.path.exists(new_path):
        now = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S%f")
        data['filename'] = new_path.replace('.nc', f'_{now}.nc')

    # Provenance
    ancestors = [data.pop('original_filename')]
    opts = [opt for opt in cfg if opt in globals()]
    caption = (f"{cube.long_name} for {mlr.get_alias(data)} preprocessed with "
               f"operations {opts}.")
    record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': caption,
        'references': ['schlund20jgr'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(data['filename'], record)

    # Write file
    io.metadata_to_netcdf(cube, data)


def main(cfg):
    """Run the diagnostic."""
    cfg = deepcopy(cfg)
    warnings.filterwarnings(
        'ignore',
        message='Using DEFAULT_SPHERICAL_EARTH_RADIUS',
        category=UserWarning,
        module='iris',
    )
    input_data = mlr.get_input_data(cfg,
                                    pattern=cfg.get('pattern'),
                                    ignore=cfg.get('ignore'),
                                    check_mlr_attributes=False)

    # Default options
    cfg.setdefault('area_weighted', True)
    cfg.setdefault('extract_ignore_bounds', False)
    cfg.setdefault('n_jobs', 1)
    cfg.setdefault('return_trend_stderr', True)
    cfg.setdefault('time_weighted', True)
    logger.info("Using at most %i processes", cfg['n_jobs'])

    # Get reference dataset for unifying coordinates if necessary
    if 'unify_coords_to' in cfg:
        ref_cube = get_ref_cube(input_data, **cfg['unify_coords_to'])
    else:
        ref_cube = None

    # Load cubes and apply common mask
    input_data = load_cubes(input_data)
    input_data = apply_common_mask(cfg, input_data)

    # Operations that add additional datasets (standard errors)
    for data in input_data:
        data.setdefault('ref', False)
        if data['ref'] == 'True':
            data['ref'] = True
        if data['ref'] == 'False':
            data['ref'] = False
        cube = data['cube']
        cube = unify_coords_to(cube, ref_cube)
        cube = mask(cfg, cube)
        cube = scalar_operations(cfg, cube)
        cube = extract_range(cfg, cube)
        cube = extract(cfg, cube)
        (cube, data) = aggregate_by_trend(cfg, cube, data)
        (cube, data) = collapse_with_trend(cfg, cube, data)
        data = cache_cube(cfg, cube, data)
    input_data = add_standard_errors(input_data)
    cfg.get('collapse', {}).pop('trend', None)
    cfg.get('aggregate_by', {}).pop('trend', None)

    # Remaining operations
    for data in input_data:
        cube = data['cube']
        (cube, data) = aggregate_by(cfg, cube, data)
        (cube, data) = collapse(cfg, cube, data)
        (cube, data) = argsort(cfg, cube, data)
        data = cache_cube(cfg, cube, data)

    # Calculations involving reference datasets
    input_data = ref_calculation(cfg, input_data)
    input_data = add_standard_errors(input_data)

    # Convert units
    for data in input_data:
        cube = data['cube']
        (cube, data) = convert_units_to(cfg, cube, data)
        data = cache_cube(cfg, cube, data)

    # Save cubes
    for data in input_data:
        cube = data.pop('cube')
        data.pop('original_cube', None)
        data.pop('ref_cube', None)
        data.pop('stderr', None)

        # Normalize and write cubes
        (cube, data) = normalize_by_mean(cfg, cube, data)
        (cube, data) = normalize_by_std(cfg, cube, data)
        write_cube(cfg, cube, data)


# Run main function when this script is called
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
