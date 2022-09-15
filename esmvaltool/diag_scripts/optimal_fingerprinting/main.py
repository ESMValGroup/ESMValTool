#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Optimal fingerprinting.

Description
-----------
Separate forced and unforced signal.

yi(t) = fi(t) + vi(t)

Assumption:

fi(t) = M(t) * bi + ci

Author
------
Manuel Schlund (DLR, Germany)

Configuration options in recipe
-------------------------------
dim_obs: str, optional (default: 'year')
    Dimension across which the different "observations" of the variable are
    taken (i.e., the dimension used to setup the linear regressions). This is
    usually some kind of time coordinate. May not be confused with
    observational data (i.e., observations of the Earth system).
mmm_facets: dict, optional
    Facets to determine the multi-model mean dataset. These are passed as
    keyword arguments to
    :func:`esmvaltool.diag_scripts.shared.select_metadata`. Exactly one
    multi-model mean dataset need to be supplied. By default,
    ``{'multi_model_statistics': 'MultiModelMean'}`` is used.

"""
import logging
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
import iris.coord_categorisation
import numpy as np
from iris.exceptions import CoordinateNotFoundError
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


def _add_categorized_time_coords(cube, coord):
    """Add categorized time coordinates to cube."""
    if cube.coords(coord):
        return
    if hasattr(iris.coord_categorisation, f'add_{coord}'):
        getattr(iris.coord_categorisation, f'add_{coord}')(cube, 'time')
        logger.debug("Added coordinate '%s' to cube", coord)
        return
    raise ValueError(
        f"Failed to add categorized time coordinate '{coord}' to cube "
        f"{cube.summary(shorten=True)}"
    )


def _check_dims(cfg, cube, cube_mmm):
    """Check if dimensions of `cube` and `cube_mmm` match."""
    dim_obs = cfg['dim_obs']

    # Check if observation dimension is present
    if not cube.coords(dim_obs):
        try:
            _add_categorized_time_coords(cube, dim_obs)
        except Exception as exc:
            raise CoordinateNotFoundError(
                f"Cube {cube.summary(shorten=True)} does not have coordinate "
                f"'{dim_obs}', make sure the option 'dim_obs' is set "
                f"correctly in the recipe"
            ) from exc
    if not cube_mmm.coords(dim_obs):
        try:
            _add_categorized_time_coords(cube_mmm, dim_obs)
        except Exception as exc:
            raise CoordinateNotFoundError(
                f"Multi-model mean cube {cube_mmm.summary(shorten=True)} does "
                f"not have coordinate '{dim_obs}', make sure the option "
                f"'dim_obs' is set correctly in the recipe"
            ) from exc

    # Check that observational coordinates are 1D
    dim_obs_coord = cube.coord(dim_obs)
    if dim_obs_coord.ndim != 1:
        raise ValueError(
            f"Expected 1D '{dim_obs}' coordinate for cube "
            f"{cube.summary(shorten=True)}, got {dim_obs_coord.ndim:d}D "
            f"coordinate"
        )
    dim_obs_coord_mmm = cube_mmm.coord(dim_obs)
    if dim_obs_coord_mmm.ndim != 1:
        raise ValueError(
            f"Expected 1D '{dim_obs}' coordinate for multi-model mean cube "
            f"{cube_mmm.summary(shorten=True)}, got "
            f"{dim_obs_coord_mmm.ndim:d}D coordinate"
        )

    # Equalize coordinates
    dim_obs_coord.var_name = dim_obs
    dim_obs_coord.long_name = dim_obs
    dim_obs_coord.standard_name = None
    dim_obs_coord.attributes = {}
    dim_obs_coord_mmm.var_name = dim_obs
    dim_obs_coord_mmm.long_name = dim_obs
    dim_obs_coord_mmm.standard_name = None
    dim_obs_coord_mmm.attributes = {}

    # It is allowed for the cube to only contain a subset of the data of the
    # MMM cube, but all points in cube need to be present in the MMM cube.
    cube_mmm_subset = cube_mmm.subset(dim_obs_coord)
    if cube_mmm_subset.coord(dim_obs).shape[0] != dim_obs_coord.shape[0]:
        raise ValueError(
            f"The points of dimension '{dim_obs}' of cube "
            f"{cube.summary(shorten=True)} must be a subset of the points of "
            f"multi-model mean cube {cube_mmm.summary(shorten=True)}, "
            f"got\n{dim_obs_coord.points}\nfor cube, and\n"
            f"{dim_obs_coord_mmm.points}\nfor the multi-model mean cube"
        )

    # Check that remaining coordinates are identical
    if cube.shape != cube_mmm_subset.shape:
        raise ValueError(
            f"Expected identical shapes for cube and multi-model mean cube "
            f"after equalizing '{dim_obs}' dimension, got {cube.shape} and "
            f"{cube_mmm_subset.shape}, respectively"
        )
    coord_dims = cube.coord_dims(dim_obs)
    for coord in cube.coords(dim_coords=True):
        if cube.coord_dims(coord) == coord_dims:
            continue
        if not cube_mmm.coords(coord):
            raise CoordinateNotFoundError(
                f"Expected dimensional coordinate '{coord.name()}' of cube "
                f"{cube.summary(shorten=True)} also in multi-model mean cube "
                f"{cube_mmm.summary(shorten=True)}"
            )
        coord_mmm = cube_mmm.coord(coord)
        if coord != coord_mmm:
            raise ValueError(
                f"Expected identical '{coord.name()}' coordinates for cube "
                f"{cube.summary(shorten=True)} and multi-model mean cube "
                f"{cube_mmm.summary(shorten=True)}, got\n{coord}\nand\n"
                f"{coord_mmm},\nrespectively"
            )

    return (cube, cube_mmm_subset)


def _get_default_cfg(cfg):
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)
    cfg.setdefault('dim_obs', 'year')
    cfg.setdefault('mmm_facets', {'multi_model_statistics': 'MultiModelMean'})

    logger.info("Using '%s' as observational dimension", cfg['dim_obs'])

    return cfg


def _get_mmm_dataset(cfg, input_data):
    """Extract multi-model mean dataset."""
    logger.info("Using facets to extract multi-model mean dataset:\n%s",
                pformat(cfg['mmm_facets']))
    datasets = select_metadata(input_data, **cfg['mmm_facets'])
    if len(datasets) != 1:
        raise ValueError(
            f"Expected exactly 1 multi-model mean dataset, found "
            f"{len(datasets):d}. Consider extending the option 'mmm_facets', "
            f"currently used:\n{pformat(cfg['mmm_facets'])}"
        )
    logger.info("Using multi-model mean dataset:\n%s", pformat(datasets[0]))
    return datasets[0]


def _get_lin_reg(x_arr, y_arr):
    """Get linear regression of two (masked) arrays."""
    if np.ma.is_masked(y_arr):
        x_arr = x_arr[~y_arr.mask]
        y_arr = y_arr[~y_arr.mask]
    if len(y_arr) < 2:
        return np.nan
    return stats.linregress(x_arr, y_arr)


def _get_calc_slope():
    """Get slope of linear regression of two (masked) arrays."""

    def calc_slope(x_arr, y_arr):
        """Calculate slope of linear regression."""
        return _get_lin_reg(x_arr, y_arr).slope

    calc_slope = np.vectorize(calc_slope, excluded=['x_arr'],
                              signature='(n),(n)->()')
    return calc_slope


def _get_calc_intercept():
    """Get intercept of linear regression of two (masked) arrays."""

    def calc_intercept(x_arr, y_arr):
        """Calculate intercept of linear regression."""
        return _get_lin_reg(x_arr, y_arr).intercept

    calc_intercept = np.vectorize(calc_intercept, excluded=['x_arr'],
                                  signature='(n),(n)->()')
    return calc_intercept


def _split_response(cfg, cube, cube_mmm):
    """Split response into forced and unforced component."""
    (cube, cube_mmm_subset) = _check_dims(cfg, cube, cube_mmm)
    coord_dim = cube.coord_dims(cfg['dim_obs'])[0]

    # Move time to rightmost position to use np.vectorize effectively
    cube_data = np.moveaxis(cube.data, coord_dim, -1)
    cube_mmm_subset_data = np.moveaxis(cube_mmm_subset.data, coord_dim, -1)

    # Perform linear regression (the actual optimal fingerprinting) assuming
    #
    # yi(t) = fi(t) + vi(t) (signal decomposition)
    #
    # with fi(t) = M(t) * bi + ci (multi-model mean is approx. of forced resp.)
    #
    # --> yi(t) = M(t) * bi + vi(t) + ci
    #
    # for each position (e.g., grid cells, or global mean). Since vi(t) is
    # error-like (e.g., mean of zero), perform linear regression between yi(t)
    # and M(t) to calculate bi and ci.
    #
    # Here,
    #
    # yi(t): Response of dataset i at time t
    # fi(t): Forced component of response of dataset i at time t
    # vi(t): Unforced component of dataset i at time t
    # M(t) : Multi-model response at time t
    # bi   : Single scaling factor of dataset i
    # ci   : Constant for dataset i
    slope = _get_calc_slope()(cube_mmm_subset_data, cube_data)
    intercept = _get_calc_intercept()(cube_mmm_subset_data, cube_data)

    # Reshape arrays so that they are broadcastable to the original shape
    slope = np.expand_dims(slope, axis=coord_dim)
    intercept = np.expand_dims(intercept, axis=coord_dim)

    # Calculate forced component
    cube_forced_comp = cube.copy()
    cube_forced_comp.data = cube_mmm_subset.data * slope + intercept
    cube_forced_comp.var_name += "_forced_comp"
    cube_forced_comp.long_name += " (Forced Component)"
    cube_forced_comp.attributes['optimal_fingerprinting'] = 'forced'

    # Calculate unforced component
    cube_unforced_comp = cube.copy()
    cube_unforced_comp.data = cube.data - cube_forced_comp.data
    cube_unforced_comp.var_name += "_unforced_comp"
    cube_unforced_comp.long_name += " (Unforced Component)"
    cube_unforced_comp.attributes['optimal_fingerprinting'] = 'unforced'

    return (cube_forced_comp, cube_unforced_comp)


def main(cfg):
    """Run diagnostic."""
    cfg = _get_default_cfg(cfg)
    input_data = list(cfg['input_data'].values())

    # Extract multi-model mean dataset that is used as proxy for forced
    # response
    dataset_mmm = _get_mmm_dataset(cfg, input_data)
    filename_mmm = dataset_mmm['filename']
    cube_mmm = iris.load_cube(filename_mmm)
    logger.info("Loaded MMM data from %s", filename_mmm)

    # Remove forced response for each dataset to separate the signal into the
    # forced and unforced component
    for dataset in input_data:
        filename = dataset['filename']
        cube = iris.load_cube(filename)
        logger.info("Processing %s", filename)
        (cube_forced_comp, cube_unforced_comp) = _split_response(
            cfg,
            cube,
            cube_mmm,
        )

        # Provencance tracking
        general_caption = (
            f"component of {dataset['long_name']} of dataset "
            f"{dataset['dataset']} (experiment {dataset['exp']}) from "
            f"{dataset['start_year']} to {dataset['end_year']}"
        )
        record = {
            'ancestors': [filename, filename_mmm],
            'authors': ['schlund_manuel'],
        }

        # Save forced component
        record['caption'] = "Forced " + general_caption
        path_forced_comp = get_diagnostic_filename(
            Path(filename).stem + '_forced_comp', cfg)
        iris.save(cube_forced_comp, path_forced_comp)
        logger.info("Wrote %s", path_forced_comp)
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(path_forced_comp, record)

        # Save unforced component
        record['caption'] = "Unforced " + general_caption
        path_unforced_comp = get_diagnostic_filename(
            Path(filename).stem + '_unforced_comp', cfg)
        iris.save(cube_unforced_comp, path_unforced_comp)
        logger.info("Wrote %s", path_unforced_comp)
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(path_unforced_comp, record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
