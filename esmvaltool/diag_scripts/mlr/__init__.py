"""Convenience functions for MLR diagnostics."""

import logging
import os
import re
import warnings
from copy import deepcopy
from pprint import pformat

import iris
import numpy as np
import shapely.vectorized as shp_vect
from cartopy.io import shapereader
from cf_units import Unit
from iris.fileformats.netcdf import UnknownCellMethodWarning

import esmvalcore.preprocessor
from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    io,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))

NECESSARY_KEYS = io.NECESSARY_KEYS + [
    'tag',
    'var_type',
]
VAR_TYPES = [
    'feature',
    'label',
    'label_to_rescale',
    'prediction_input',
    'prediction_input_error',
    'prediction_output',
    'prediction_output_error',
    'prediction_output_misc',
    'prediction_reference',
    'prediction_residual',
]
WARNINGS_TO_IGNORE = [
    {
        'message': ".* contains unknown cell method 'trend'",
        'category': UnknownCellMethodWarning,
        'module': 'iris',
    },
    {
        'message': "Using DEFAULT_SPHERICAL_EARTH_RADIUS",
        'category': UserWarning,
        'module': 'iris',
    },
]


def _check_coords(cube, coords, weights_type):
    """Check coordinates prior to weights calculations."""
    cube_str = cube.summary(shorten=True)
    for coord_name in coords:
        try:
            coord = cube.coord(coord_name)
        except iris.exceptions.CoordinateNotFoundError:
            logger.error(
                "Calculation of %s for cube %s failed, coordinate "
                "'%s' not found", weights_type, cube_str, coord_name)
            raise
        if not coord.has_bounds():
            logger.debug(
                "Guessing bounds of coordinate '%s' of cube %s for "
                "calculation of %s", coord_name, cube_str, weights_type)
            coord.guess_bounds()


def _get_datasets(input_data, **kwargs):
    """Get datasets according to ``**kwargs``."""
    datasets = []
    for dataset in input_data:
        dataset_copy = deepcopy(dataset)
        for key in kwargs:
            if key not in dataset_copy:
                dataset_copy[key] = None
        if select_metadata([dataset_copy], **kwargs):
            datasets.append(dataset)
    return datasets


def _get_ne_land_mask_cube(n_lats=1000, n_lons=2000):
    """Get Natural Earth land mask."""
    ne_dir = os.path.join(
        os.path.dirname(os.path.realpath(esmvalcore.preprocessor.__file__)),
        'ne_masks',
    )
    ne_file = os.path.join(ne_dir, 'ne_10m_land.shp')
    reader = shapereader.Reader(ne_file)
    geometries = list(reader.geometries())

    # Setup grid
    lat_coord = iris.coords.DimCoord(
        np.linspace(-90.0, 90.0, n_lats), var_name='lat',
        standard_name='latitude', long_name='latitude', units='degrees')
    lon_coord = iris.coords.DimCoord(
        np.linspace(-180.0, 180.0, n_lons), var_name='lon',
        standard_name='longitude', long_name='longitude', units='degrees')
    (lats, lons) = np.meshgrid(lat_coord.points, lon_coord.points)

    # Setup mask (1: land, 0: sea)
    mask = np.full(lats.shape, False, dtype=bool)
    for geometry in geometries:
        mask |= shp_vect.contains(geometry, lons, lats)
    land_mask = np.swapaxes(np.where(mask, 1, 0), 0, 1)

    # Setup cube
    cube = iris.cube.Cube(land_mask,
                          var_name='land_mask',
                          long_name='Land mask (1: land, 0: sea)',
                          units='no_unit',
                          dim_coords_and_dims=[(lat_coord, 0), (lon_coord, 1)])
    return cube


def _has_valid_coords(cube, coords):
    """Check if cube has valid coords for calculating weights."""
    for coord_name in coords:
        try:
            cube.coord(coord_name)
        except iris.exceptions.CoordinateNotFoundError:
            return False
    return True


def check_predict_kwargs(predict_kwargs):
    """Check keyword argument for ``predict()`` functions.

    Parameters
    ----------
    predict_kwargs : keyword arguments, optional
        Keyword arguments for a ``predict()`` function.

    Raises
    ------
    RuntimeError
        ``return_var`` and ``return_cov`` are both set to ``True`` in the
        keyword arguments.

    """
    return_var = predict_kwargs.get('return_var', False)
    return_cov = predict_kwargs.get('return_cov', False)
    if return_var and return_cov:
        raise RuntimeError(
            "Cannot return variance (return_cov=True) and full covariance "
            "matrix (return_cov=True) simultaneously")


def create_alias(dataset, attributes, delimiter='-'):
    """Create alias key of a dataset using a list of attributes.

    Parameters
    ----------
    dataset : dict
        Metadata dictionary representing a single dataset.
    attributes : list of str
        List of attributes used to create the alias.
    delimiter : str, optional (default: '-')
        Delimiter used to separate different attributes in the alias.

    Returns
    -------
    str
        Dataset alias.

    Raises
    ------
    AttributeError
        ``dataset`` does not contain one of the ``attributes``.

    """
    alias = []
    if not attributes:
        raise ValueError(
            "Expected at least one element for attributes, got empty list")
    for attribute in attributes:
        if attribute not in dataset:
            raise AttributeError(
                f"Dataset {dataset} does not contain attribute '{attribute}' "
                f"for alias creation")
        alias.append(dataset[attribute])
    return delimiter.join(alias)


def datasets_have_mlr_attributes(datasets, log_level='debug', mode='full'):
    """Check (MLR) attributes of ``datasets``.

    Parameters
    ----------
    datasets : list of dict
        Datasets to check.
    log_level : str, optional (default: 'debug')
        Verbosity level of the logger.
    mode : str, optional (default: 'full')
        Checking mode. Must be one of ``'only_missing'`` (only check if
        attributes are missing), ``'only_var_type'`` (check only `var_type`) or
        ``'full'`` (check both).

    Returns
    -------
    bool
        ``True`` if all required attributes are available, ``False`` if not.

    Raises
    ------
    ValueError
        Invalid value for argument ``mode`` is given.

    """
    output = True
    accepted_modes = ('full', 'only_missing', 'only_var_type')
    if mode not in accepted_modes:
        raise ValueError(
            f"'mode' must be one of {accepted_modes}, got '{mode}'")
    for dataset in datasets:
        if mode != 'only_var_type':
            for key in NECESSARY_KEYS:
                if key not in dataset:
                    getattr(logger, log_level)(
                        "Dataset '%s' does not have necessary (MLR) attribute "
                        "'%s'", dataset, key)
                    output = False
        if mode != 'only_missing' and dataset.get('var_type') not in VAR_TYPES:
            getattr(logger, log_level)(
                "Dataset '%s' has invalid var_type '%s', must be one of %s",
                dataset, dataset.get('var_type'), VAR_TYPES)
            output = False
    return output


def get_1d_cube(x_data, y_data, x_kwargs=None, y_kwargs=None):
    """Convert 2 arrays to :class:`iris.cube.Cube` (with single coordinate).

    Parameters
    ----------
    x_data : numpy.ndarray
        Data for coordinate.
    y_data : numpy.ndarray
        Data for cube.
    x_kwargs : dict
        Keyword arguments passed to :class:`iris.coords.AuxCoord`.
    y_kwargs : dict
        Keyword arguments passed to :class:`iris.cube.Cube`.

    Returns
    -------
    iris.cube.Cube
        1D cube with single auxiliary coordinate.

    Raises
    ------
    ValueError
        Arrays are not 1D and do not have matching shapes.

    """
    if x_kwargs is None:
        x_kwargs = {}
    if y_kwargs is None:
        y_kwargs = {}
    x_data = np.ma.array(x_data)
    y_data = np.ma.array(y_data)
    if x_data.ndim != 1:
        raise ValueError(
            f"Expected 1D array for 'x_data', got {x_data.ndim:d}D array")
    if y_data.ndim != 1:
        raise ValueError(
            f"Expected 1D array for 'y_data', got {y_data.ndim:d}D array")
    if x_data.shape != y_data.shape:
        raise ValueError(
            f"Expected identical shapes for 'x_data' and 'y_data', got "
            f"{x_data.shape} and {y_data.shape}, respectively")
    aux_coord = iris.coords.AuxCoord(x_data, **x_kwargs)
    cube = iris.cube.Cube(y_data, aux_coords_and_dims=[(aux_coord, 0)],
                          **y_kwargs)
    return cube


def get_absolute_time_units(units):
    """Convert time reference units to absolute ones.

    This function converts reference time units (like ``'days since YYYY'``) to
    absolute ones (like ``'days'``).

    Parameters
    ----------
    units : cf_units.Unit
        Time units to convert.

    Returns
    -------
    cf_units.Unit
        Absolute time units.

    Raises
    ------
    ValueError
        If conversion failed (e.g. input units are not time units).

    """
    if units.is_time_reference():
        units = Unit(units.symbol.split()[0])
    if not units.is_time():
        raise ValueError(
            f"Cannot convert units '{units}' to reasonable time units")
    return units


def get_alias(dataset):
    """Get alias for dataset.

    Parameters
    ----------
    dataset : dict
        Dataset metadata.

    Returns
    -------
    str
        Alias.

    """
    alias = f"{dataset['project']} dataset {dataset['dataset']}"
    additional_info = []
    for key in ('mip', 'exp', 'ensemble'):
        if key in dataset:
            additional_info.append(dataset[key])
    if additional_info:
        alias += f" ({', '.join(additional_info)})"
    if 'start_year' in dataset and 'end_year' in dataset:
        alias += f" from {dataset['start_year']:d} to {dataset['end_year']:d}"
    return alias


def get_all_weights(cube, area_weighted=True, time_weighted=True,
                    landsea_fraction_weighted=None, normalize=False):
    """Get all possible weights of cube.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube.
    area_weighted : bool, optional (default: True)
        Use area weights.
    time_weighted : bool, optional (default: True)
        Use time weights.
    landsea_fraction_weighted : str, optional
        If given, use land/sea fraction weights. Must be one of ``'land'``,
        ``'sea'``.
    normalize : bool, optional (default: False)
        Normalize weights with total area and total time range.

    Returns
    -------
    numpy.ndarray
        Area weights.

    """
    logger.debug("Calculating all weights of cube %s",
                 cube.summary(shorten=True))
    weights = np.ones(cube.shape)

    # Horizontal weights
    if _has_valid_coords(cube, ['latitude', 'longitude']):
        horizontal_weights = get_horizontal_weights(
            cube, area_weighted=area_weighted,
            landsea_fraction_weighted=landsea_fraction_weighted,
            normalize=normalize)
        if horizontal_weights is not None:
            weights *= horizontal_weights

    # Time weights
    if _has_valid_coords(cube, ['time']) and time_weighted:
        time_weights = get_time_weights(cube, normalize=normalize)
        weights *= time_weights

    return weights


def get_area_weights(cube, normalize=False):
    """Get area weights of cube.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube.
    normalize : bool, optional (default: False)
        Normalize weights with total area.

    Returns
    -------
    numpy.ndarray
        Area weights.

    Raises
    ------
    iris.exceptions.CoordinateNotFoundError
        Cube does not contain the coordinates ``latitude`` and ``longitude``.

    """
    logger.debug("Calculating area weights")
    _check_coords(cube, ['latitude', 'longitude'], 'area weights')
    area_weights = iris.analysis.cartography.area_weights(cube,
                                                          normalize=normalize)
    return area_weights


def get_horizontal_weights(cube, area_weighted=True,
                           landsea_fraction_weighted=None, normalize=False):
    """Get horizontal weights of cube.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube.
    area_weighted : bool, optional (default: True)
        Use area weights.
    landsea_fraction_weighted : str, optional
        If given, use land/sea fraction weights. Must be one of ``'land'``,
        ``'sea'``.
    normalize : bool, optional (default: False)
        Normalize weights with sum of weights over latitude and longitude (i.e.
        if only ``area_weighted`` is given, this is equal to the total area).

    Returns
    -------
    numpy.ndarray
        Area weights.

    Raises
    ------
    iris.exceptions.CoordinateMultiDimError
        Dimension of latitude or longitude coordinate is greater than 1.
    iris.exceptions.CoordinateNotFoundError
        Cube does not contain the coordinates ``latitude`` and ``longitude``.
    ValueError
        ``landsea_fraction_weighted`` is not one of ``'land'``, ``'sea'``.

    """
    logger.debug("Calculating horizontal weights")
    if not (area_weighted or landsea_fraction_weighted):
        return None

    # Get weights
    weights = np.ones(cube.shape)
    if area_weighted:
        weights *= get_area_weights(cube, normalize=False)
    if landsea_fraction_weighted is not None:
        weights *= get_landsea_fraction_weights(
            cube, landsea_fraction_weighted, normalize=False)

    # No normalization
    if not normalize:
        return weights

    # Get horizontal dimensions
    horizontal_dims = []
    if cube.coord_dims('latitude'):
        horizontal_dims.append(cube.coord_dims('latitude')[0])
    if cube.coord_dims('longitude'):
        horizontal_dims.append(cube.coord_dims('longitude')[0])

    # Normalization
    horizontal_dims = tuple(horizontal_dims)
    if not horizontal_dims:
        norm = np.ravel(weights)[0]
    else:
        norm = np.ravel(np.sum(weights, axis=horizontal_dims))[0]
    return weights / norm


def get_input_data(cfg, pattern=None, check_mlr_attributes=True, ignore=None):
    """Get input data and check MLR attributes if desired.

    Use ``input_data`` and ancestors to get all relevant input files.

    Parameters
    ----------
    cfg : dict
        Recipe configuration.
    pattern : str, optional
        Pattern matched against ancestor file names.
    check_mlr_attributes : bool, optional (default: True)
        If ``True``, only returns datasets with valid MLR attributes. If
        ``False``, returns all found datasets.
    ignore : list of dict, optional
        Ignore specific datasets by specifying multiple :obj:`dict`s of
        metadata. By setting an attribute to ``None``, ignore all datasets
        which do not have that attribute.

    Returns
    -------
    list of dict
        List of input datasets.

    Raises
    ------
    ValueError
        No input data found or at least one dataset has invalid attributes.

    """
    logger.debug("Extracting input files")
    input_data = list(cfg['input_data'].values())
    input_data.extend(io.netcdf_to_metadata(cfg, pattern=pattern))
    input_data = deepcopy(input_data)
    if ignore is not None:
        valid_data = []
        ignored_datasets = []
        logger.info("Ignoring files with %s", ignore)
        for kwargs in ignore:
            ignored_datasets.extend(_get_datasets(input_data, **kwargs))
        for dataset in input_data:
            if dataset not in ignored_datasets:
                valid_data.append(dataset)
    else:
        valid_data = input_data
    if not valid_data:
        raise ValueError("No input data found")
    if check_mlr_attributes:
        if not datasets_have_mlr_attributes(valid_data, log_level='error'):
            raise ValueError("At least one input dataset does not have valid "
                             "MLR attributes")
    logger.debug("Found files:")
    logger.debug(pformat([d['filename'] for d in valid_data]))
    return valid_data


def get_landsea_fraction_weights(cube, area_type, normalize=False):
    """Get land/sea fraction weights of cube using Natural Earth files.

    Note
    ----
    The implementation of this feature is not optimal. For large cubes,
    calculating the land/sea fraction weights might be very slow.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube.
    area_type : str
        Area type. Must be one of ``'land'`` (land fraction weighting) or
        ``'sea'`` (sea fraction weighting).
    normalize : bool, optional (default: False)
        Normalize weights with total land/sea fraction.

    Raises
    ------
    iris.exceptions.CoordinateMultiDimError
        Dimension of latitude or longitude coordinate is greater than 1.
    iris.exceptions.CoordinateNotFoundError
        Cube does not contain the coordinates ``latitude`` and ``longitude``.
    ValueError
        ``area_type`` is not one of ``'land'``, ``'sea'``.

    """
    allowed_types = ('land', 'sea')
    if area_type not in allowed_types:
        raise ValueError(
            f"Expected one of {allowed_types} for 'area_type' of land/sea "
            f"fraction weighting, got '{area_type}'")
    logger.debug("Calculating %s fraction weights", area_type)
    _check_coords(cube, ['latitude', 'longitude'],
                  f'{area_type} fraction weights')
    lat_coord = cube.coord('latitude')
    lon_coord = cube.coord('longitude')
    for coord in (lat_coord, lon_coord):
        if coord.ndim > 1:
            raise iris.exceptions.CoordinateMultiDimError(
                f"Calculating {area_type} fraction weights for "
                f"multidimensional coordinate '{coord.name}' is not supported")

    # Calculate land fractions on coordinate grid of cube
    ne_land_mask_cube = _get_ne_land_mask_cube()
    land_fraction = np.empty((lat_coord.shape[0], lon_coord.shape[0]),
                             dtype=np.float)
    for lat_idx in range(lat_coord.shape[0]):
        for lon_idx in range(lon_coord.shape[0]):
            lat_bounds = lat_coord.bounds[lat_idx]
            lon_bounds = lon_coord.bounds[lon_idx]
            submask = ne_land_mask_cube.intersection(latitude=lat_bounds,
                                                     longitude=lon_bounds)
            land_fraction[lat_idx, lon_idx] = (submask.data.sum() /
                                               submask.data.size)
    if area_type == 'sea':
        fraction_weights = 1.0 - land_fraction
    else:
        fraction_weights = land_fraction
    if normalize:
        fraction_weights /= np.ma.sum(fraction_weights)

    # Broadcast to original shape
    coord_dims = []
    if cube.coord_dims(lon_coord):
        coord_dims.append(cube.coord_dims(lon_coord)[0])
    else:
        fraction_weights = np.squeeze(fraction_weights, axis=1)
    if cube.coord_dims(lat_coord):
        coord_dims.insert(0, cube.coord_dims(lat_coord)[0])
    else:
        fraction_weights = np.squeeze(fraction_weights, axis=0)
    fraction_weights = iris.util.broadcast_to_shape(fraction_weights,
                                                    cube.shape,
                                                    tuple(coord_dims))

    return fraction_weights


def get_new_path(cfg, old_path):
    """Convert old path to new diagnostic path.

    Parameters
    ----------
    cfg : dict
        Recipe configuration.
    old_path : str
        Old path.

    Returns
    -------
    str
        New diagnostic path.

    """
    basename = os.path.splitext(os.path.basename(old_path))[0]
    new_path = get_diagnostic_filename(basename, cfg)
    return new_path


def get_squared_error_cube(ref_cube, error_datasets):
    """Get array of squared errors.

    Parameters
    ----------
    ref_cube : iris.cube.Cube
        Reference cube (determines mask, coordinates and attributes of output).
    error_datasets : list of dict
        List of metadata dictionaries where each dictionary represents a single
        dataset.

    Returns
    -------
    iris.cube.Cube
        Cube containing squared errors.

    Raises
    ------
    ValueError
        Shape of a dataset does not match shape of reference cube.

    """
    squared_error_cube = ref_cube.copy()

    # Fill cube with zeros
    squared_error_cube.data = np.ma.array(
        np.full(squared_error_cube.shape, 0.0),
        mask=np.ma.getmaskarray(squared_error_cube.data),
    )

    # Adapt cube metadata
    if 'error' in squared_error_cube.attributes.get('var_type', ''):
        if not squared_error_cube.attributes.get('squared'):
            squared_error_cube.var_name += '_squared'
            squared_error_cube.long_name += ' (squared)'
            squared_error_cube.units = units_power(squared_error_cube.units, 2)
    else:
        if squared_error_cube.attributes.get('squared'):
            squared_error_cube.var_name += '_error'
            squared_error_cube.long_name += ' (error)'
        else:
            squared_error_cube.var_name += '_squared_error'
            squared_error_cube.long_name += ' (squared error)'
            squared_error_cube.units = units_power(squared_error_cube.units, 2)
    squared_error_cube.attributes['squared'] = 1
    squared_error_cube.attributes['var_type'] = 'prediction_output_error'

    # Aggregate errors
    filenames = []
    for dataset in error_datasets:
        path = dataset['filename']
        cube = iris.load_cube(path)
        filenames.append(path)

        # Check shape
        if cube.shape != ref_cube.shape:
            raise ValueError(
                f"Expected shape {ref_cube.shape} for error cubes, got "
                f"{cube.shape} for dataset '{path}'")

        # Add squared error
        new_data = cube.data
        if not cube.attributes.get('squared'):
            new_data **= 2
        squared_error_cube.data += new_data
        logger.debug("Added '%s' to squared error datasets", path)
    squared_error_cube.attributes['filename'] = '|'.join(filenames)
    return squared_error_cube


def get_time_weights(cube, normalize=False):
    """Get time weights of cube.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube.
    normalize : bool, optional (default: False)
        Normalize weights with total time range.

    Returns
    -------
    numpy.ndarray
        Time weights.

    Raises
    ------
    iris.exceptions.CoordinateNotFoundError
        Cube does not contain the coordinate ``time``.

    """
    logger.debug("Calculating time weights")
    _check_coords(cube, ['time'], 'time weights')
    coord = cube.coord('time')
    time_weights = coord.bounds[:, 1] - coord.bounds[:, 0]
    time_weights = time_weights.squeeze()
    if normalize:
        time_weights /= np.ma.sum(time_weights)
    if time_weights.shape == ():
        time_weights = np.broadcast_to(time_weights, cube.shape)
    else:
        time_weights = iris.util.broadcast_to_shape(time_weights, cube.shape,
                                                    cube.coord_dims('time'))
    return time_weights


def ignore_warnings():
    """Ignore warnings given by ``WARNINGS_TO_IGNORE``."""
    for warning_kwargs in WARNINGS_TO_IGNORE:
        warning_kwargs.setdefault('action', 'ignore')
        warnings.filterwarnings(**warning_kwargs)


def square_root_metadata(cube):
    """Take the square root of the cube metadata.

    Parameters
    ----------
    cube : iris.cube.Cube
        Cube (will be modified in-place).

    """
    if 'squared_' in cube.var_name:
        cube.var_name = cube.var_name.replace('squared_', '')
    elif '_squared' in cube.var_name:
        cube.var_name = cube.var_name.replace('_squared', '')
    else:
        cube.var_name = 'root_' + cube.var_name
    if 'squared ' in cube.long_name:
        cube.long_name = cube.long_name.replace('squared ', '')
    elif 'Squared ' in cube.long_name:
        cube.long_name = cube.long_name.replace('Squared ', '')
    elif ' squared' in cube.long_name:
        cube.long_name = cube.long_name.replace(' squared', '')
    elif ' Squared' in cube.long_name:
        cube.long_name = cube.long_name.replace(' Squared', '')
    elif ' (squared)' in cube.long_name:
        cube.long_name = cube.long_name.replace(' (squared)', '')
    elif ' (Squared)' in cube.long_name:
        cube.long_name = cube.long_name.replace(' (Squared)', '')
    else:
        cube.long_name = 'Root ' + cube.long_name
    cube.units = cube.units.root(2)
    if cube.attributes.get('squared'):
        cube.attributes.pop('squared')


def units_power(units, power):
    """Raise a :class:`cf_units.Unit` to given power preserving symbols.

    Raise :class:`cf_units.Unit` to given power without expanding it first. For
    example, using ``units_power(Unit('J'), 2)`` gives ``Unit('J2')``. In
    contrast, simply using ``Unit('J')**2`` would yield ``'kg2 m4 s-4'``.

    Parameters
    ----------
    units : cf_units.Unit
        Input units.
    power : int
        Desired exponent.

    Returns
    -------
    cf_units.Unit
        Input units raised to given power.

    Raises
    ------
    TypeError
        Argument ``power`` is not :obj:`int`-like.
    ValueError
        Invalid unit given.

    """
    if round(power) != power:
        raise TypeError(
            f"Expected integer-like power for units exponentiation, got "
            f"{power}")
    power = int(power)
    if any([units.is_no_unit(), units.is_unknown()]):
        raise ValueError(
            f"Cannot raise units '{units.name}' to power {power:d}")
    if units.origin is None:
        logger.warning(
            "Symbol-preserving exponentiation of units '%s' is not "
            "supported, origin is not given", units)
        return units**power
    if units.origin.isdigit():
        return units**power
    if units.origin.split()[0][0].isdigit():
        logger.warning(
            "Symbol-preserving exponentiation of units '%s' is not "
            "supported yet because of leading numbers", units)
        return units**power
    new_units_list = []
    for split in units.origin.split():
        for elem in split.split('.'):
            if elem[-1].isdigit():
                exp = [int(d) for d in re.findall(r'-?\d+', elem)][0]
                val = ''.join(list(re.findall(r'[A-Za-z]', elem)))
                new_units_list.append(f'{val}{exp * power}')
            else:
                new_units_list.append(f'{elem}{power}')
    new_units = ' '.join(new_units_list)
    return Unit(new_units)
