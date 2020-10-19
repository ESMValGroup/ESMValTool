#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simple postprocessing of MLR model output.

Description
-----------
This diagnostic performs postprocessing operations for MLR model output (mean
and error).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Notes
-----
Prior to postprocessing, this diagnostic groups input datasets according to
``tag`` and ``prediction_name``. For each group, accepts datasets with three
different ``var_type`` s:

* ``prediction_output``: **Exactly one** necessary, refers to
  the mean prediction and serves as reference dataset (regarding shape).
* ``prediction_output_error``: Arbitrary number of error datasets. If not
  given, error calculation is skipped. May be squared errors (marked by the
  attribute ``squared``) or not. In addition, a single covariance dataset can
  be specified (``short_name`` ending with ``_cov``).
* ``prediction_input``: Dataset used to estimate covariance structure of
  the mean prediction (i.e. matrix of Pearson correlation coefficients) for
  error estimation. At most one dataset allowed. Ignored when no
  ``prediction_output_error`` is given. This is only possible when (1) the
  shape of the ``prediction_input`` dataset is identical to the shape of the
  ``prediction_output_error`` datasets, (2) the number of dimensions of the
  ``prediction_input`` dataset is higher than the number of dimensions of the
  ``prediction_output_error`` datasets and they have identical trailing
  (rightmost) dimensions or (3) the number of dimensions of the
  ``prediction_input`` dataset is higher than the number of dimensions of
  ``prediction_output_error`` datasets and all dimensions of the
  ``prediction_output_error`` datasets are mapped to a corresponding dimension
  of the ``prediction_input`` using the ``cov_estimate_dim_map`` option (e.g.
  when ``prediction_input`` has shape ``(10, 5, 100, 20)`` and
  ``prediction_output_error`` has shape ``(5, 20)``, you can use
  ``cov_estimate_dim_map: [1, 3]`` to map the dimensions of
  ``prediction_output_error`` to dimension 1 and 3 of ``prediction_input``).

All data with other ``var_type`` s is ignored (``feature``, ``label``, etc.).

Real error calculation (using covariance dataset given as
``prediction_output_error``) and estimation (using ``prediction_input`` dataset
to estimate covariance structure) is only possible if the mean prediction cube
is collapsed completely during postprocessing, i.e. all coordinates are listed
for either ``mean`` or ``sum``.

Configuration options in recipe
-------------------------------
add_var_from_cov: bool, optional (default: True)
    Calculate variances from covariance matrix (diagonal elements) and add
    those to (squared) error datasets. Set to ``False`` if variance is already
    given separately in prediction output.
area_weighted: bool, optional (default: True)
    Calculate weighted averages/sums when collapsing over latitude and/or
    longitude coordinates using grid cell areas (calculated using grid cell
    boundaries). Only possible if the datasets contains ``latitude`` and
    ``longitude`` coordinates.
convert_units_to: str, optional
    Convert units of the input data.
cov_estimate_dim_map: list of int, optional
    Map dimensions of ``prediction_output_error`` datasets to corresponding
    dimensions of ``prediction_input`` used for estimating covariance. Only
    relevant if both dataset types are given. See notes above for more
    information.
ignore: list of dict, optional
    Ignore specific datasets by specifying multiple :obj:`dict` s of metadata.
landsea_fraction_weighted: str, optional
    When given, calculate weighted averages/sums when collapsing over latitude
    and/or longitude coordinates using land/sea fraction (calculated using
    Natural Earth masks). Only possible if the datasets contains ``latitude``
    and ``longitude`` coordinates. Must be one of ``'land'``, ``'sea'``.
mean: list of str, optional
    Perform mean over the given coordinates.
pattern: str, optional
    Pattern matched against ancestor file names.
sum: list of str, optional
    Perform sum over the given coordinates.
time_weighted: bool, optional (default: True)
    Calculate weighted averages/sums for time (using grid cell boundaries).

"""

import logging
import os
import warnings
from copy import deepcopy
from pprint import pformat

import iris
import numpy as np
from cf_units import Unit

from esmvaltool.diag_scripts import mlr
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    io,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))

OPS = {
    'mean': iris.analysis.MEAN,
    'sum': iris.analysis.SUM,
}


def _calculate_lower_error_bound(cfg, squared_error_cube, basepath):
    """Calculate lower error bound."""
    logger.debug("Calculating lower error bound")
    lower_bound = _collapse_regular_cube(cfg, squared_error_cube, power=2)
    lower_bound.data = np.ma.sqrt(lower_bound.data)
    mlr.square_root_metadata(lower_bound)
    _convert_units(cfg, lower_bound)
    lower_bound.attributes['error_type'] = 'lower_bound'
    new_path = basepath.replace('.nc', '_lower_bound.nc')
    io.iris_save(lower_bound, new_path)
    logger.info("Lower bound of error: %s %s", lower_bound.data,
                lower_bound.units)
    ancestors = _get_ancestors(squared_error_cube)
    _write_provenance(cfg, 'Lower bound of', lower_bound, new_path, ancestors)


def _calculate_real_error(cfg, ref_cube, cov_cube, basepath):
    """Calculate real error using covariance."""
    logger.debug("Calculating real error using covariance")
    real_error = _collapse_covariance_cube(cfg, cov_cube, ref_cube)
    real_error.data = np.ma.sqrt(real_error.data)
    real_error.var_name = cov_cube.var_name.replace('_cov', '_error')
    real_error.long_name = cov_cube.long_name.replace('(covariance)',
                                                      '(error)')
    real_error.units = real_error.units.root(2)
    _convert_units(cfg, real_error)
    real_error.attributes['error_type'] = 'real_error'
    new_path = basepath.replace('.nc', '_real.nc')
    io.iris_save(real_error, new_path)
    logger.info("Real error (using covariance): %s %s", real_error.data,
                real_error.units)
    ancestors = _get_ancestors(cov_cube)
    _write_provenance(cfg, 'Real', real_error, new_path, ancestors)


def _calculate_upper_error_bound(cfg, squared_error_cube, basepath):
    """Calculate upper error bound."""
    logger.debug("Calculating upper error bound")
    upper_bound = squared_error_cube.copy()
    upper_bound.data = np.ma.sqrt(upper_bound.data)
    mlr.square_root_metadata(upper_bound)
    upper_bound = _collapse_regular_cube(cfg, upper_bound)
    _convert_units(cfg, upper_bound)
    upper_bound.attributes['error_type'] = 'upper_bound'
    new_path = basepath.replace('.nc', '_upper_bound.nc')
    io.iris_save(upper_bound, new_path)
    logger.info("Upper bound of error: %s %s", upper_bound.data,
                upper_bound.units)
    ancestors = _get_ancestors(squared_error_cube)
    _write_provenance(cfg, 'Upper bound of', upper_bound, new_path, ancestors)


def _convert_units(cfg, cube):
    """Convert units if desired."""
    cfg_settings = cfg.get('convert_units_to')
    if cfg_settings:
        units_to = cfg_settings
        logger.debug("Converting units from '%s' to '%s'", cube.units,
                     units_to)
        try:
            cube.convert_units(units_to)
        except ValueError:
            raise ValueError(
                f"Cannot convert units of cube {cube.summary(shorten=True)} "
                f"from '{cube.units}' to '{units_to}'")


def _collapse_covariance_cube(cfg, cov_cube, ref_cube):
    """Collapse covariance cube with using desired operations."""
    (weights, units, coords) = _get_all_weights(cfg, ref_cube)
    if len(coords) < ref_cube.ndim:
        raise ValueError(
            f"Calculating real error using covariance dataset "
            f"('prediction_output_error') is only possible if all "
            f"{ref_cube.ndim:d} dimensions of the cube are collapsed, got "
            f"only {len(coords):d} ({coords})")
    weights = weights.ravel()
    weights = weights[~np.ma.getmaskarray(ref_cube.data).ravel()]
    weights = np.outer(weights, weights)
    cov_cube = cov_cube.collapsed(cov_cube.coords(dim_coords=True),
                                  iris.analysis.SUM,
                                  weights=weights)
    cov_cube.units *= units**2
    return cov_cube


def _collapse_regular_cube(cfg, cube, power=1):
    """Collapse cube with using desired operations."""
    (weights, units, coords) = _get_all_weights(cfg, cube, power=power)
    cube = cube.collapsed(coords, iris.analysis.SUM, weights=weights)
    cube.units *= units
    return cube


def _corrcoef(array, rowvar=True, weights=None):
    """Fast version of :func:`numpy.ma.corrcoef`."""
    if not rowvar:
        array = array.T
        if weights is not None:
            weights = weights.T
    mean = np.ma.average(array, axis=1, weights=weights).reshape(-1, 1)
    if weights is None:
        sqrt_weights = 1.0
    else:
        sqrt_weights = np.ma.sqrt(weights)
    demean = (array - mean) * sqrt_weights
    res = np.ma.dot(demean, demean.T)
    row_norms = np.ma.sqrt(np.ma.sum(demean**2, axis=1))
    res /= np.ma.outer(row_norms, row_norms)
    return res


def _estim_cov_differing_shape(cfg, squared_error_cube, cov_est_cube, weights):
    """Collapse estimated covariance.

    Estimate error by estimating covariance from dataset with at least one
    dimension more than errors themselves.

    """
    logger.info(
        "Estimating true error from covariance derived from "
        "'prediction_input' dataset with shape %s for errors "
        "('prediction_output_error') with shape %s", cov_est_cube.shape,
        squared_error_cube.shape)

    # Load data
    error = np.ma.sqrt(squared_error_cube.data)
    error = np.ma.filled(error, 0.0)
    cov_est = np.ma.filled(cov_est_cube.data, np.nan)

    # Reshape if necessary
    if 'cov_estimate_dim_map' in cfg:
        dim_map = tuple(cfg['cov_estimate_dim_map'])
        cov_est = _reshape_covariance(cov_est, error, dim_map)
    cov_est = np.ma.masked_invalid(cov_est)
    if not _identical_trailing_dimensions(cov_est, error):
        raise ValueError(
            f"Expected identical trailing (rightmost) dimensions of "
            f"'prediction_input' data used to estimate covariance structure "
            f"and 'prediction_output_error' datasets, got {cov_est.shape} and "
            f"{error.shape}")

    # Estimate covariance
    error = error.ravel()
    cov_est = cov_est.reshape(-1, *error.shape)
    pearson_coeffs = _corrcoef(cov_est, rowvar=False)
    covariance = pearson_coeffs * np.ma.outer(error, error)

    # Collapse covariance
    weights = weights.ravel()
    weights = np.outer(weights, weights)
    error = np.ma.sqrt(np.ma.sum(covariance * weights))
    return error


def _estim_cov_identical_shape(squared_error_cube, cov_est_cube, weights):
    """Collapse estimated covariance.

    Estimate error by approximating covariance from dataset with identical
    shape as errors themselves (for a better estimate, the dataset used to
    estimate covariance needs at least one dimension more than the errors).

    """
    logger.info(
        "Estimating true error from covariance derived from "
        "'prediction_input' dataset with same shape as errors "
        "('prediction_output_error')")
    error = np.ma.sqrt(squared_error_cube.data)
    error = np.ma.filled(error, 0.0)
    cov_est = np.ma.array(cov_est_cube.data)
    if cov_est.ndim > 2:
        error = error.reshape(error.shape[0], -1)
        cov_est = cov_est.reshape(cov_est.shape[0], -1)
        weights = weights.reshape(weights.shape[0], -1)

    # Pearson coefficients (= normalized covariance) over both dimensions
    pearson_dim0 = _corrcoef(cov_est, weights=weights)
    pearson_dim1 = _corrcoef(cov_est, rowvar=False, weights=weights)

    # Covariances
    cov_dim0 = (np.einsum('...i,...j->...ij', error, error) *
                np.einsum('...i,...j->...ij', weights, weights) * pearson_dim1)
    cov_dim1 = (np.einsum('i...,j...->...ij', error, error) *
                np.einsum('i...,j...->...ij', weights, weights) * pearson_dim0)

    # Errors over dimensions
    error_dim0 = np.ma.sqrt(np.ma.sum(cov_dim0, axis=(1, 2)))
    error_dim1 = np.ma.sqrt(np.ma.sum(cov_dim1, axis=(1, 2)))

    # Collaps further (all weights are already included in first step)
    cov_order_0 = pearson_dim0 * np.ma.outer(error_dim0, error_dim0)
    cov_order_1 = pearson_dim1 * np.ma.outer(error_dim1, error_dim1)
    error_order_0 = np.ma.sqrt(np.ma.sum(cov_order_0))
    error_order_1 = np.ma.sqrt(np.ma.sum(cov_order_1))
    logger.debug(
        "Found real errors %e and %e after collapsing with different "
        "orderings, using maximum", error_order_0, error_order_1)

    # Ordering of collapsing matters, maximum is used
    return max([error_order_0, error_order_1])


def _estimate_real_error(cfg, squared_error_cube, cov_est_dataset, basepath):
    """Estimate real error using estimated covariance."""
    logger.debug(
        "Estimating real error using estimated covariance from "
        "'prediction_input' dataset %s", cov_est_dataset['filename'])
    cov_est_cube = iris.load_cube(cov_est_dataset['filename'])

    # Check dimensions
    if cov_est_cube.ndim < 2:
        raise ValueError(
            f"Expected at least 2D 'prediction_input' dataset for covariance "
            f"structure estimation, got {cov_est_cube.ndim:d}D dataset")
    if cov_est_cube.ndim < squared_error_cube.ndim:
        raise ValueError(
            f"Expected number of dimensions of 'prediction_input' dataset "
            f"used for covariance structure estimation to be greater than or "
            f"equal the number of dimensions of the errors datasets, got "
            f"{cov_est_cube.ndim:d} and {squared_error_cube.ndim}")

    # Check if all dimensions are collapsed
    (weights, units, coords) = _get_all_weights(cfg, squared_error_cube)
    if len(coords) < squared_error_cube.ndim:
        raise ValueError(
            f"Estimating real error using 'prediction_input' dataset for "
            f"covariance structure estimation is only possible if all "
            f"{squared_error_cube.ndim:d} dimensions of the error cube are "
            f"collapsed, got only {len(coords):d} ({coords})")

    # Estimate error
    if cov_est_cube.shape == squared_error_cube.shape:
        error = _estim_cov_identical_shape(squared_error_cube, cov_est_cube,
                                           weights)
    else:
        error = _estim_cov_differing_shape(cfg, squared_error_cube,
                                           cov_est_cube, weights)

    # Create cube (collapse using dummy operation)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message="Collapsing spatial coordinate 'latitude' without "
                    "weighting",
            category=UserWarning,
            module='iris',
        )
        real_error = squared_error_cube.collapsed(coords, iris.analysis.MEAN)
    real_error.data = error
    mlr.square_root_metadata(real_error)
    real_error.units *= units
    _convert_units(cfg, real_error)
    real_error.attributes['error_type'] = 'estimated_real_error'

    # Save cube
    new_path = basepath.replace('.nc', '_estimated.nc')
    io.iris_save(real_error, new_path)
    logger.info("Estimated real error (using estimated covariance): %s %s",
                real_error.data, real_error.units)

    # Provenance
    ancestors = [
        *_get_ancestors(squared_error_cube),
        cov_est_dataset['filename'],
    ]
    _write_provenance(cfg, 'Estimated', real_error, new_path, ancestors)


def _get_all_weights(cfg, cube, power=1):
    """Get all necessary weights (including norm for mean calculation)."""
    cfg = deepcopy(cfg)
    all_coords = []
    weights = np.ones(cube.shape)
    units = Unit('1')

    # Iterate over operations
    for operation in ('sum', 'mean'):
        normalize = (operation == 'mean')
        coords = cfg.get(operation, [])
        all_coords.extend(coords)
        if coords == 'all':
            coords = [c.name() for c in cube.coords(dim_coords=True)]
        horizontal_coords = _get_horizontal_coordinates(coords)

        # Horizontal coordinates
        if horizontal_coords:
            (horizontal_weights, area_units) = _get_horizontal_weights(
                cfg, cube, power=power)
            if operation == 'sum':
                units *= area_units
            if horizontal_weights is not None:
                weights *= horizontal_weights
            weights /= _get_normalization_factor(
                horizontal_weights, horizontal_coords, cube,
                normalize=normalize)**power
            for coord in horizontal_coords:
                coords.remove(coord)

        # Time coordinate
        if 'time' in coords:
            (time_weights, time_units) = _get_time_weights(cfg, cube,
                                                           power=power)
            if operation == 'sum':
                units *= time_units
            if time_weights is not None:
                weights *= time_weights
            weights /= _get_normalization_factor(
                time_weights, ['time'], cube, normalize=normalize)**power
            coords.remove('time')

        # Remaining coordinates
        weights /= _get_normalization_factor(
            None, coords, cube, normalize=normalize)**power

    # Apply mask of cube to weights
    weights = np.ma.array(weights, mask=np.ma.getmaskarray(cube.data))
    logger.debug("Found coordinates %s to collapse over", all_coords)
    logger.debug("Found units '%s' for weights", units)
    return (weights, units, all_coords)


def _get_ancestors(cube):
    """Extract ancestors from ``filename`` attribute of cube."""
    ancestors = cube.attributes['filename'].split('|')
    return ancestors


def _get_covariance_dataset(error_datasets, ref_cube):
    """Extract covariance dataset."""
    explanation = ("i.e. dataset with short_name == '*_cov' among "
                   "'prediction_output_error' datasets")
    cov_datasets = []
    other_datasets = []

    # Get covariance dataset(s)
    for dataset in error_datasets:
        if '_cov' in dataset['short_name']:
            cov_datasets.append(dataset)
        else:
            other_datasets.append(dataset)
    if not cov_datasets:
        logger.warning(
            "No covariance dataset (%s) found, calculation of real error not "
            "possible", explanation)
        return (None, other_datasets)
    if len(cov_datasets) > 1:
        filenames = [d['filename'] for d in cov_datasets]
        raise ValueError(
            f"Expected at most one covariance dataset ({explanation}), got "
            f"{len(cov_datasets):d}:\n{pformat(filenames)}")

    # Check shape
    cov_cube = iris.load_cube(cov_datasets[0]['filename'])
    cov_cube.attributes['filename'] = cov_datasets[0]['filename']
    ref_size = np.ma.array(ref_cube.data).compressed().shape[0]
    if cov_cube.shape != (ref_size, ref_size):
        raise ValueError(
            f"Expected shape of covariance dataset to be "
            f"{(ref_size, ref_size)}, got {cov_cube.shape} (after removal of "
            f"all missing values)")
    return (cov_cube, other_datasets)


def _get_horizontal_weights(cfg, cube, power=1):
    """Calculate weights for horizontal coordinates."""
    weights = mlr.get_horizontal_weights(
        cube,
        area_weighted=cfg['area_weighted'],
        landsea_fraction_weighted=cfg.get('landsea_fraction_weighted'),
    )
    if weights is not None:
        weights = weights**power
    if cfg['area_weighted']:
        units = Unit('m2')**power
    else:
        units = Unit('1')**power
    return (weights, units)


def _get_horizontal_coordinates(coords):
    """Extract horizontal coordinates from :obj:`list` of coordinates."""
    horizontal_coords = []
    if 'latitude' in coords:
        horizontal_coords.append('latitude')
    if 'longitude' in coords:
        horizontal_coords.append('longitude')
    return horizontal_coords


def _get_normalization_factor(weights, coords, cube, normalize=False):
    """Get normalization constant for calculation of means."""
    if not normalize:
        return 1.0
    if not coords:
        return 1.0
    if weights is None:
        weights = np.ones(cube.shape)
    weights = np.ma.array(weights.copy(), mask=np.ma.getmaskarray(cube.data))
    coord_dims = []
    for coord in coords:
        if cube.coord_dims(coord):
            coord_dims.extend(cube.coord_dims(coord))
    coord_dims = tuple(coord_dims)
    if not coord_dims:
        norm = np.ma.ravel(weights)
    else:
        norm = np.ma.ravel(np.ma.sum(weights, axis=coord_dims))
    norm = norm[~np.ma.getmaskarray(norm)]
    return norm[0]


def _get_time_weights(cfg, cube, power=1):
    """Calculate time weights."""
    time_weights = None
    time_units = mlr.get_absolute_time_units(cube.coord('time').units)
    if cfg.get['time_weighted']:
        time_weights = mlr.get_time_weights(cube)
        time_weights = time_weights**power
    return (time_weights, time_units**power)


def _identical_trailing_dimensions(larger_array, smaller_array):
    """Check if trailing dimensions of two arrays are identical."""
    if larger_array.ndim < smaller_array.ndim:
        raise ValueError(
            f"Expected array with higher number of dimensions as first "
            f"argument, got {larger_array.ndim:d}D array as first argument, "
            f"{smaller_array.ndim:d}D array as second")
    return larger_array.shape[-smaller_array.ndim:] == smaller_array.shape


def _write_provenance(cfg, title, error_cube, netcdf_path, ancestors):
    """Write provenance record."""
    caption = f'{title} {error_cube.long_name}'
    attrs = error_cube.attributes
    if 'mlr_model_name' in attrs:
        caption += f" of MLR model {attrs['mlr_model_name']}"
    if 'prediction_name' in attrs:
        caption += f" for prediction {attrs['prediction_name']}"
    caption += '.'
    record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': caption,
        'references': ['schlund20jgr'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, record)


def check_cfg(cfg):
    """Check options of configuration and catch errors."""
    for operation in ('sum', 'mean'):
        if operation in cfg:
            cfg[operation] = list(set(cfg[operation]))
    for coord in cfg.get('sum', []):
        if coord in cfg.get('mean', []):
            raise ValueError(f"Coordinate '{coord}' given in 'sum' and 'mean'")


def postprocess_errors(cfg, ref_cube, error_datasets, cov_estim_datasets):
    """Postprocess errors."""
    logger.info(
        "Postprocessing errors using mean prediction cube %s as reference",
        ref_cube.summary(shorten=True))

    # Extract covariance
    (cov_cube,
     error_datasets) = _get_covariance_dataset(error_datasets, ref_cube)

    # Extract squared errors
    squared_error_cube = mlr.get_squared_error_cube(ref_cube, error_datasets)

    # Extract variance from covariance if desired
    if cfg.get('add_var_from_cov', True) and cov_cube is not None:
        var = np.ma.empty(ref_cube.shape, dtype=ref_cube.dtype)
        mask = np.ma.getmaskarray(ref_cube.data)
        var[mask] = np.ma.masked
        var[~mask] = np.diagonal(cov_cube.data.copy())
        squared_error_cube.data += var
        logger.debug(
            "Added variance calculated from covariance to squared error "
            "datasets")
        if not error_datasets:
            error_datasets = True

    # Extract basename for error cubes
    basepath = mlr.get_new_path(cfg, ref_cube.attributes['filename'])
    basepath = basepath.replace('.nc', '_error.nc')

    # Lower and upper error bounds
    if error_datasets:
        _calculate_lower_error_bound(cfg, squared_error_cube, basepath)
        _calculate_upper_error_bound(cfg, squared_error_cube, basepath)

        # Estimated real error using estimated covariance
        if cov_estim_datasets:
            _estimate_real_error(cfg, squared_error_cube,
                                 cov_estim_datasets[0], basepath)

    # Real error
    if cov_cube is not None:
        _calculate_real_error(cfg, ref_cube, cov_cube, basepath)


def postprocess_mean(cfg, cube, data):
    """Postprocess mean prediction cube."""
    logger.info("Postprocessing mean prediction cube %s",
                cube.summary(shorten=True))
    cube = _collapse_regular_cube(cfg, cube)
    _convert_units(cfg, cube)
    new_path = mlr.get_new_path(cfg, data['filename'])
    io.iris_save(cube, new_path)
    logger.info("Mean prediction: %s %s", cube.data, cube.units)
    _write_provenance(cfg, "Postprocessed", cube, new_path, [data['filename']])


def _reshape_covariance(cov_est, error, dim_map):
    """Reshape covariance estimation input to match errors."""
    if len(dim_map) != error.ndim:
        raise ValueError(
            f"Dimension mapping for covariance estimation "
            f"'cov_estimate_dim_map' needs to cover all dimensions of "
            f"'prediction_output_error' cubes with shape {error.shape}, "
            f"got {dim_map}")
    if len(set(dim_map)) != len(dim_map):
        raise ValueError(
            f"Duplicate dimensions in 'cov_estimate_dim_map' are not "
            f"allowed, got {dim_map}")
    logger.debug(
        "Reshaping 'prediction_input' with shape %s to contain dimensions "
        "of 'prediction_output_error' %s as last dimensions using mapping "
        "%s", cov_est.shape, error.shape, dim_map)

    # Get target shape
    indices = list(range(cov_est.ndim))
    for dim in dim_map:
        if dim not in indices:
            raise ValueError(
                f"Dimensional index {dim:d} in 'cov_estimate_dim_map' is out "
                f"of range for {cov_est.ndim:d}D 'prediction_input' dataset "
                f"used for covariance estimation")
        indices.remove(dim)
        indices.append(dim)
    new_shape = tuple(np.array(cov_est.shape)[indices])

    # Broadcast to new shape
    dim_map_for_broadcasting = []
    for idx in range(cov_est.ndim):
        dim_map_for_broadcasting.append(indices.index(idx))
    cov_est = iris.util.broadcast_to_shape(cov_est, new_shape,
                                           dim_map_for_broadcasting)
    logger.info(
        "Reshaped 'prediction_input' for covariance estimation to %s",
        cov_est.shape)
    return cov_est


def split_datasets(datasets, tag, pred_name):
    """Split datasets into mean and error."""
    grouped_data = group_metadata(datasets, 'var_type')

    # Mean/reference dataset
    mean = grouped_data.get('prediction_output', [])
    if len(mean) != 1:
        filenames = [d['filename'] for d in mean]
        raise ValueError(
            f"Expected exactly one 'prediction_output' dataset for tag "
            f"'{tag}' of prediction '{pred_name}', got {len(mean):d}:\n"
            f"{pformat(filenames)}")
    logger.info(
        "Found mean prediction dataset ('prediction_output') for tag '%s' of "
        "prediction '%s': %s (used as reference)", tag, pred_name,
        mean[0]['filename'])

    # Errors
    error = grouped_data.get('prediction_output_error', [])
    if not error:
        logger.warning(
            "No 'prediction_output_error' datasets for tag '%s' of prediction "
            "'%s' found, error calculation not possible (not searching for "
            "'prediction_input' datasets for covariance estimation, either)",
            tag, pred_name)
        cov_estimation = []
    else:
        logger.info(
            "Found error datasets ('prediction_output_error') for tag '%s' of "
            "prediction '%s':", tag, pred_name)
        logger.info(pformat([d['filename'] for d in error]))

        # Estimation for covariance
        cov_estimation = grouped_data.get('prediction_input', [])
        if not cov_estimation:
            logger.warning(
                "No 'prediction_input' dataset for tag '%s' of prediction "
                "'%s' found, real error estimation using estimated covariance "
                "structure not possible", tag, pred_name)
        elif len(cov_estimation) > 1:
            filenames = [d['filename'] for d in cov_estimation]
            raise ValueError(
                f"Expected at most one 'prediction_input' dataset for tag "
                f"'{tag}' of prediction '{pred_name}', got "
                f"{len(cov_estimation):d}:\n{pformat(filenames)}")
        else:
            logger.info(
                "Found 'prediction_input' dataset for covariance structure "
                "estimation for tag '%s' of prediction '%s': %s", tag,
                pred_name, cov_estimation[0]['filename'])

    return (mean[0], error, cov_estimation)


def main(cfg):
    """Run the diagnostic."""
    warnings.filterwarnings(
        'ignore',
        message='Using DEFAULT_SPHERICAL_EARTH_RADIUS',
        category=UserWarning,
        module='iris',
    )
    input_data = mlr.get_input_data(cfg,
                                    pattern=cfg.get('pattern'),
                                    ignore=cfg.get('ignore'))

    # Check cfg
    check_cfg(cfg)
    cfg.setdefault('area_weighted', True)
    cfg.setdefault('time_weighted', True)

    # Process data
    for (tag, tag_datasets) in group_metadata(input_data, 'tag').items():
        logger.info("Processing tag '%s'", tag)
        grouped_data = group_metadata(tag_datasets, 'prediction_name')
        for (pred_name, datasets) in grouped_data.items():
            logger.info("Processing prediction '%s'", pred_name)
            (dataset, error_datasets,
             cov_estim_datastets) = split_datasets(datasets, tag, pred_name)

            # Extract cubes
            logger.debug(
                "Loaded mean prediction cube from '%s' (used as reference)",
                dataset['filename'])
            cube = iris.load_cube(dataset['filename'])
            cube.attributes['filename'] = dataset['filename']
            if cube.ndim < 1:
                raise ValueError(
                    f"Postprocessing scalar dataset '{dataset['filename']}' "
                    f"not supported yet")

            # Process mean prediction
            postprocess_mean(cfg, cube, dataset)

            # Process errors
            postprocess_errors(cfg, cube, error_datasets, cov_estim_datastets)


# Run main function when this script is called
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
