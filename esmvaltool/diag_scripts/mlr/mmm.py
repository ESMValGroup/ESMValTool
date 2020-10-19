#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Use simple multi-model mean for predictions.

Description
-----------
This diagnostic calculates the (unweighted) mean over all given datasets for a
given target variable.

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
convert_units_to: str, optional
    Convert units of the input data. Can also be given as dataset option.
dtype: str (default: 'float64')
    Internal data type which is used for all calculations, see
    `<https://docs.scipy.org/doc/numpy/user/basics.types.html>`_ for a list of
    allowed values.
ignore: list of dict, optional
    Ignore specific datasets by specifying multiple :obj:`dict` s of metadata.
mlr_model_name: str, optional (default: 'MMM')
    Human-readable name of the MLR model instance (e.g used for labels).
mmm_error_type: str, optional
    If given, additionally saves estimated squared MMM model error. If the
    option is set to ``'loo'``, the (constant) error is estimated as RMSEP
    using leave-one-out cross-validation. No other options are supported at the
    moment.
pattern: str, optional
    Pattern matched against ancestor file names.
prediction_name: str, optional
    Default ``prediction_name`` of output cubes if no 'prediction_reference'
    dataset is given.
weighted_samples: dict
    If specified, use weighted mean square error to estimate prediction error.
    The given keyword arguments are directly passed to
    :func:`esmvaltool.diag_scripts.mlr.get_all_weights` to calculate the sample
    weights. By default, area weights and time weights are used.

"""

import logging
import os
from copy import deepcopy
from pprint import pformat

import iris
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import LeaveOneOut

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts import mlr
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    io,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))


def _add_dataset_attributes(cube, datasets, cfg):
    """Add dataset-related attributes to cube."""
    dataset_names = sorted(list({d['dataset'] for d in datasets}))
    projects = sorted(list({d['project'] for d in datasets}))
    start_years = list({d['start_year'] for d in datasets})
    end_years = list({d['end_year'] for d in datasets})
    cube.attributes['dataset'] = '|'.join(dataset_names)
    cube.attributes['description'] = 'MMM prediction'
    cube.attributes['end_year'] = min(end_years)
    cube.attributes['mlr_model_name'] = cfg['mlr_model_name']
    cube.attributes['mlr_model_type'] = 'mmm'
    cube.attributes['project'] = '|'.join(projects)
    cube.attributes['start_year'] = min(start_years)
    cube.attributes['var_type'] = 'prediction_output'


def _load_cube(cfg, dataset):
    """Load single :class:`iris.cube.Cube`."""
    path = dataset['filename']
    cube = iris.load_cube(path)
    cube.data = cube.core_data().astype(cfg['dtype'], casting='same_kind')
    convert_units(cfg, cube, dataset)
    return (cube, path)


def add_general_attributes(cube, **kwargs):
    """Add general attributes to cube."""
    for (key, val) in kwargs.items():
        if val is not None:
            cube.attributes[key] = val


def convert_units(cfg, cube, data):
    """Convert units if desired."""
    cfg_settings = cfg.get('convert_units_to')
    data_settings = data.get('convert_units_to')
    if cfg_settings or data_settings:
        units_to = cfg_settings
        if data_settings:
            units_to = data_settings
        logger.info("Converting units from '%s' to '%s'", cube.units, units_to)
        cube.convert_units(units_to)


def get_loo_error_cube(cfg, label_datasets):
    """Estimate prediction error using cross-validation."""
    loo = LeaveOneOut()
    logger.info("Estimating prediction error using cross-validator %s",
                str(loo.__class__))
    label_datasets = np.array(label_datasets)
    errors = []
    for (train_idx, test_idx) in loo.split(label_datasets):
        ref_cube = get_mmm_cube(cfg, label_datasets[test_idx])
        mmm_cube = get_mmm_cube(cfg, label_datasets[train_idx])

        # Apply mask
        mask = np.ma.getmaskarray(ref_cube.data).ravel()
        mask |= np.ma.getmaskarray(mmm_cube.data).ravel()

        y_true = ref_cube.data.ravel()[~mask]
        y_pred = mmm_cube.data.ravel()[~mask]
        weights = mlr.get_all_weights(ref_cube, **cfg['weighted_samples'])
        weights = weights.ravel()[~mask]

        # Calculate mean squared error
        error = mean_squared_error(y_true, y_pred, sample_weight=weights)
        errors.append(error)

    # Get error cube
    error_cube = get_mmm_cube(cfg, label_datasets)
    error_array = np.empty(error_cube.shape).ravel()
    mask = np.ma.getmaskarray(error_cube.data).ravel()
    error_array[mask] = np.nan
    error_array[~mask] = np.mean(errors)
    error_array = np.ma.masked_invalid(error_array)
    error_cube.data = error_array.reshape(error_cube.shape)

    # Cube metadata
    error_cube.attributes['error_type'] = 'loo'
    error_cube.attributes['squared'] = 1
    error_cube.attributes['var_type'] = 'prediction_output_error'
    error_cube.var_name += '_squared_mmm_error_estim'
    error_cube.long_name += ' (squared MMM error estimation using CV)'
    error_cube.units = mlr.units_power(error_cube.units, 2)
    return error_cube


def get_grouped_data(cfg, input_data=None):
    """Get input files."""
    if input_data is None:
        logger.debug("Loading input data from 'cfg' argument")
        input_data = mlr.get_input_data(cfg,
                                        pattern=cfg.get('pattern'),
                                        ignore=cfg.get('ignore'))
    else:
        logger.debug("Loading input data from 'input_data' argument")
        if not mlr.datasets_have_mlr_attributes(input_data, log_level='error'):
            raise ValueError("At least one input dataset does not have valid "
                             "MLR attributes")
    if not input_data:
        raise ValueError("No input data found")
    paths = [d['filename'] for d in input_data]
    logger.debug("Found files")
    logger.debug(pformat(paths))

    # Extract necessary data
    label_data = select_metadata(input_data, var_type='label')
    if not label_data:
        raise ValueError("No data with var_type 'label' found")
    prediction_reference_data = select_metadata(
        input_data, var_type='prediction_reference')
    extracted_data = label_data + prediction_reference_data
    logger.debug("Found 'label' data")
    logger.debug(pformat([d['filename'] for d in label_data]))
    logger.debug("Found 'prediction_reference' data")
    logger.debug(pformat([d['filename'] for d in prediction_reference_data]))

    # Return grouped data
    return group_metadata(extracted_data, 'tag')


def get_mmm_cube(cfg, label_datasets):
    """Get multi-model mean data."""
    cubes = iris.cube.CubeList()
    paths = []
    (ref_cube, _) = _load_cube(cfg, label_datasets[0])
    for dataset in label_datasets:
        (cube, path) = _load_cube(cfg, dataset)
        ih.prepare_cube_for_merging(cube, path)
        cubes.append(cube)
        paths.append(path)
    mmm_cube = cubes.merge_cube()
    if len(paths) > 1:
        mmm_cube = mmm_cube.collapsed(['cube_label'], iris.analysis.MEAN)
    for aux_coord in ref_cube.coords(dim_coords=False):
        mmm_cube.add_aux_coord(aux_coord, ref_cube.coord_dims(aux_coord))
    mmm_cube.remove_coord('cube_label')
    _add_dataset_attributes(mmm_cube, label_datasets, cfg)
    return mmm_cube


def get_reference_dataset(datasets, tag):
    """Get ``prediction_reference`` dataset."""
    ref_datasets = select_metadata(datasets, var_type='prediction_reference')
    if not ref_datasets:
        logger.warning(
            "Calculating residuals for '%s' not possible, no "
            "'prediction_reference' dataset given", tag)
        return (None, None)
    if len(ref_datasets) > 1:
        filenames = [d['filename'] for d in ref_datasets]
        raise ValueError(
            f"Expected at most one 'prediction_reference' dataset for "
            f"'{tag}', got {len(ref_datasets):d}:\n{pformat(filenames)}")
    return (ref_datasets[0], ref_datasets[0].get('prediction_name'))


def get_residual_cube(mmm_cube, ref_cube):
    """Calculate residuals."""
    if mmm_cube.shape != ref_cube.shape:
        raise ValueError(
            f"Expected identical shapes for 'label' and "
            f"'prediction_reference' datasets, got {mmm_cube.shape} and "
            f"{ref_cube.shape}, respectively")
    res_cube = ref_cube.copy()
    res_cube.data -= mmm_cube.data
    res_cube.attributes = mmm_cube.attributes
    res_cube.attributes['residuals'] = 'true minus predicted values'
    res_cube.attributes['var_type'] = 'prediction_residual'
    res_cube.var_name += '_residual'
    res_cube.long_name += ' (residual)'
    return res_cube


def save_error(cfg, label_datasets, mmm_path, **cube_attrs):
    """Save estimated error of MMM."""
    if len(label_datasets) < 2:
        logger.warning(
            "Estimating MMM prediction error not possible, at least 2 'label' "
            "datasets are needed, only %i is given", len(label_datasets))
        return
    error_type = cfg['mmm_error_type']
    allowed_error_types = ['loo']
    logger.info("Calculating error using error type '%s'", error_type)
    if error_type == 'loo':
        err_cube = get_loo_error_cube(cfg, label_datasets)
    else:
        raise NotImplementedError(
            f"mmm_error_type '{error_type}' is currently not supported, "
            f"supported types are {allowed_error_types}")
    add_general_attributes(err_cube, **cube_attrs)
    err_path = mmm_path.replace('_prediction', '_squared_prediction_error')
    io.iris_save(err_cube, err_path)
    write_provenance(cfg, err_path,
                     [d['filename'] for d in label_datasets],
                     f"{err_cube.long_name} of MMM model "
                     f"{cfg['mlr_model_name']} using error type {error_type}.")


def save_residuals(cfg, mmm_cube, ref_dataset, label_datasets, **cube_attrs):
    """Save residuals."""
    logger.info("Calculating residuals")
    (ref_cube, _) = _load_cube(cfg, ref_dataset)
    res_cube = get_residual_cube(mmm_cube, ref_cube)
    add_general_attributes(res_cube, **cube_attrs)
    mmm_path = mmm_cube.attributes['filename']
    res_path = mmm_path.replace('_prediction', '_prediction_residual')
    io.iris_save(res_cube, res_path)
    ancestors = ([d['filename'] for d in label_datasets] +
                 [ref_dataset['filename']])
    caption = (f"Residuals of predicted {res_cube.long_name} of MMM model "
               f"{cfg['mlr_model_name']}")
    if 'prediction_name' in cube_attrs:
        caption += f" for prediction {cube_attrs['prediction_name']}"
    caption += '.'
    write_provenance(cfg, res_path, ancestors, caption)


def write_provenance(cfg, netcdf_path, ancestors, caption):
    """Write provenance information."""
    record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': caption,
        'references': ['schlund20jgr'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, record)


def main(cfg, input_data=None, description=None):
    """Run the diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('dtype', 'float64')
    cfg.setdefault('mlr_model_name', 'MMM')
    cfg.setdefault('weighted_samples',
                   {'area_weighted': True, 'time_weighted': True})

    # Get data
    grouped_data = get_grouped_data(cfg, input_data=input_data)
    description = '' if description is None else f'_for_{description}'

    # Loop over all tags
    for (tag, datasets) in grouped_data.items():
        logger.info("Processing label '%s'", tag)

        # Get label datasets and reference dataset if possible
        label_datasets = select_metadata(datasets, var_type='label')
        (ref_dataset, pred_name) = get_reference_dataset(datasets, tag)
        if pred_name is None:
            pred_name = cfg.get('prediction_name')

        # Calculate multi-model mean
        logger.info("Calculating multi-model mean")
        mmm_cube = get_mmm_cube(cfg, label_datasets)
        add_general_attributes(mmm_cube, tag=tag, prediction_name=pred_name)
        mmm_path = get_diagnostic_filename(
            f"mmm_{tag}_prediction{description}", cfg)
        io.iris_save(mmm_cube, mmm_path)
        write_provenance(cfg, mmm_path,
                         [d['filename'] for d in label_datasets],
                         f"Predicted {mmm_cube.long_name} of MMM model "
                         f"{cfg['mlr_model_name']}.")

        # Estimate prediction error using cross-validation
        if 'mmm_error_type' in cfg:
            save_error(cfg, label_datasets, mmm_path, tag=tag,
                       prediction_name=pred_name)

        # Calculate residuals
        if ref_dataset is not None:
            save_residuals(cfg, mmm_cube, ref_dataset, label_datasets, tag=tag,
                           prediction_name=pred_name)


# Run main function when this script is called
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
