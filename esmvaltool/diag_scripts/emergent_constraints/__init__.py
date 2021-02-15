"""Convenience functions for emergent constraints diagnostics."""
import logging
import os
from copy import deepcopy
from pprint import pformat

import iris
import iris.pandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy import integrate
from scipy.stats import linregress

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    io,
)

logger = logging.getLogger(__name__)

ALLOWED_VAR_TYPES = [
    'feature',
    'label',
    'prediction_input',
    'prediction_input_error',
]
COLOR_COMBINED_GROUPS = 'gray'
LEGEND_KWARGS = {
    'loc': 'center left',
    'bbox_to_anchor': [1.05, 0.5],
    'borderaxespad': 0.0,
}
PANDAS_PRINT_OPTIONS = ['display.max_rows', None, 'display.max_colwidth', -1]


def _check_x_y_arrays(x_array, y_array):
    """Ensure that the X and Y arrays have correct shapes."""
    x_array = np.ma.array(x_array)
    y_array = np.ma.array(y_array)

    # Check shapes
    if x_array.ndim != 1:
        raise ValueError(
            f"Expected 1D array for X training data, got {x_array.ndim:d}D "
            f"array")
    if y_array.ndim != 1:
        raise ValueError(
            f"Expected 1D array for Y training data, got {y_array.ndim:d}D "
            f"array")
    if x_array.shape != y_array.shape:
        raise ValueError(
            f"Expected identical shapes for X and Y training data, got "
            f"{x_array.shape} and {y_array.shape}, respectively")

    # Remove masked values
    mask = np.ma.getmaskarray(x_array)
    mask |= np.ma.getmaskarray(y_array)
    x_array = np.array(x_array[~mask])
    y_array = np.array(y_array[~mask])

    return (x_array, y_array)


def _add_column(data_frame, series, column_name):
    """Add column to :class:`pandas.DataFrame` (expands index if necessary)."""
    for row in series.index.difference(data_frame.index):
        data_frame = data_frame.append(pd.Series(name=row))
    if column_name in data_frame.columns:
        for row in series.index:
            if np.isnan(data_frame.loc[row, column_name]):
                data_frame.loc[row, column_name] = series.loc[row]
            else:
                if not np.isclose(data_frame.loc[row, column_name],
                                  series.loc[row]):
                    raise ValueError(
                        f"Got duplicate data for tag '{column_name}' of "
                        f"'{row}': {series.loc[row]:e} and "
                        f"{data_frame.loc[row, column_name]:e}")
    else:
        data_frame[column_name] = series
    return data_frame


def _crop_data_frame(data_frame, ref_data_frame, data_name, ref_data_name):
    """Crop columns of a data_frame so that it matches a given reference."""
    diff_not_in_data_frame = list(
        ref_data_frame.columns.difference(data_frame.columns))
    if diff_not_in_data_frame:
        raise ValueError(
            f"No '{data_name}' given for tags {diff_not_in_data_frame}")
    diff_not_in_ref = list(
        data_frame.columns.difference(ref_data_frame.columns))
    if diff_not_in_ref:
        logger.warning(
            "Ignoring '%s' of tags %s: no corresponding '%s' data available",
            data_name, diff_not_in_ref, ref_data_name)
        data_frame = data_frame[ref_data_frame.columns]
    return data_frame


def _check_data_frames(features, label, pred_input, pred_input_err):
    """Check indexes and columns of the input data."""
    if not list(features.columns):
        raise ValueError("Expected at least one feature")
    if len(label.columns) != 1:
        raise ValueError(
            f"Expected exactly 1 'label' variable, got {len(label.columns):d}")

    # Compare features and label
    if list(features.index) != list(label.index):
        raise ValueError(
            f"Expected identical datasets (climate models; independent "
            f"observations) for 'feature' and 'label', got "
            f"{features.index.values} and {label.index.values}")
    if len(features.index) < 2:
        raise ValueError("Expected at least two training points for features")

    # Compare features and prediction input data
    pred_input = _crop_data_frame(pred_input, features, 'prediction_input',
                                  'feature')
    pred_input_err = _crop_data_frame(pred_input_err, features,
                                      'prediction_input_error', 'feature')

    # Compare prediction_input and prediction_input_error
    if not list(pred_input.index):
        raise ValueError("Expected at least one prediction input point")
    if list(pred_input.index) != list(pred_input_err.index):
        raise ValueError(
            f"Expected identical training points for 'prediction_input' and "
            f"'prediction_input_error', got {pred_input.index.values} "
            f"and {pred_input_err.index.values}")

    return (features, label, pred_input, pred_input_err)


def _combine_dicts(old_dict, new_dict):
    """Combine two :obj:`dict` by adding values for identical keys to lists."""
    old_dict = deepcopy(old_dict)
    new_dict = deepcopy(new_dict)
    for (key, val) in new_dict.items():
        if key not in old_dict:
            old_dict[key] = val
            continue
        if isinstance(old_dict[key], list):
            if not isinstance(val, list):
                val = [val]
            old_dict[key] = list(set([*old_dict[key], *val]))
        else:
            if not isinstance(val, list):
                if not np.array_equal(val, old_dict[key]):
                    old_dict[key] = [old_dict[key], val]
            else:
                old_dict[key] = list(set([old_dict[key], *val]))
    return old_dict


def _get_additional_data(additional_data, recipe):
    """Get :class:`iris.cube.CubeList` from additional data."""
    if additional_data is None:
        return iris.cube.CubeList()
    cubes = _metadata_list_to_cube_list(additional_data, 'additional_data')
    for cube in cubes:
        cube.attributes['filename'] = recipe
    return cubes


def _get_attributes(cubes):
    """Extract attributes for features and labels."""
    attributes = {}

    # Extract attributes
    for cube in cubes:
        cube_attrs = cube.attributes
        tag = cube_attrs['tag']
        attributes.setdefault(tag, {})
        if cube_attrs['var_type'] in ('feature', 'label'):
            attributes[tag] = _combine_dicts(attributes[tag],
                                             _metadata_to_dict(cube.metadata))
        elif cube_attrs['var_type'] in ('prediction_input',
                                        'prediction_input_error'):
            attributes[tag] = _combine_dicts(
                attributes[tag], {'filename': cube_attrs['filename']})
        else:
            raise ValueError(
                f"File '{cube_attrs['filename']}' has invalid var_type "
                f"'{cube_attrs['var_type']}'")

    # Set default attributes and remove lengthy 'provenance' entry
    for tag in attributes:
        attributes[tag].pop('provenance', None)
        attributes[tag].setdefault('plot_title', f"Emergent constraint {tag}")
        if 'units' in attributes[tag]:
            axis_label = f"{tag} [{attributes[tag]['units']}]"
            attributes[tag].setdefault('plot_xlabel', axis_label)
            attributes[tag].setdefault('plot_ylabel', axis_label)
        else:
            attributes[tag].setdefault('plot_xlabel', tag)
            attributes[tag].setdefault('plot_ylabel', tag)
        attributes[tag].setdefault('plot_xlim', None)
        attributes[tag].setdefault('plot_ylim', None)

    return attributes


def _get_cube_list(input_files,
                   recipe,
                   additional_data=None,
                   external_file=None,
                   merge_identical_pred_input=True):
    """Get :class:`iris.cube.CubeList` of input files."""
    cubes = iris.cube.CubeList()

    # Input files
    for filename in input_files:
        logger.info("Loading '%s'", filename)
        cube = _load_cube_with_dataset_coord(filename)
        cube.attributes['filename'] = filename
        (feature_cube,
         prediction_cube) = _split_cube(cube, merge_identical_pred_input)
        if feature_cube is not None:
            cubes.append(feature_cube)
        if prediction_cube is not None:
            cubes.append(prediction_cube)

    # Additional data
    cubes.extend(_get_additional_data(additional_data, recipe))

    # External file
    cubes.extend(_get_external_cube_list(external_file))

    # Check metadata of cubes
    for cube in cubes:
        check_metadata(cube.attributes)

    return cubes


def _get_external_cube_list(external_file):
    """Get external :class:`iris.cube.CubeList`."""
    if external_file is None:
        return iris.cube.CubeList()
    with open(external_file, 'r') as infile:
        metadata_list = yaml.safe_load(infile)
    cubes = _metadata_list_to_cube_list(metadata_list, external_file)
    for cube in cubes:
        cube.attributes['filename'] = external_file
    return cubes


def _get_external_file(filepath, auxiliary_data_dir):
    """Get full path to external file (if available)."""
    if not filepath:
        return None
    filepath = os.path.expanduser(os.path.expandvars(filepath))
    if not os.path.isabs(filepath):
        filepath = os.path.join(auxiliary_data_dir, filepath)
    if not os.path.isfile(filepath):
        raise FileNotFoundError(
            f"Desired external file '{filepath}' does not exist")
    logger.info("Found external file '%s'", filepath)
    return filepath


def _get_data_frame(var_type, cubes, label_all_data, group_by=None):
    """Extract :class:`pandas.DataFrame` for a given ``var_type``."""
    data_frame = pd.DataFrame()
    for cube in cubes:
        cube_attrs = cube.attributes
        if var_type != cube_attrs['var_type']:
            continue
        if var_type in ('feature', 'label'):
            if group_by is not None and group_by not in cube_attrs:
                raise AttributeError(
                    f"Group attribute '{group_by}' not available in input "
                    f"file '{cube_attrs['filename']}'")
            group = cube_attrs.get(group_by, label_all_data)
            index = pd.MultiIndex.from_product(
                [[group], cube.coord('dataset').points],
                names=[group_by, 'dataset'])
        else:
            index = cube.coord('dataset').points
        series = pd.Series(data=cube.data, index=index)
        data_frame = _add_column(data_frame, series, cube_attrs['tag'])
    return data_frame


def _metadata_to_dict(metadata):
    """Convert :class:`iris.cube.CubeMetadata` to :obj:`dict`."""
    new_dict = {}
    for (key, val) in metadata._asdict().items():
        if isinstance(val, dict):
            new_dict.update(val)
        else:
            new_dict[key] = val
    return new_dict


def _split_cube(cube, merge_identical_pred_input=True):
    """Split cube in features and prediction_input."""
    if not cube.attributes.get('reference_dataset'):
        return (cube, None)

    # Get feature and prediction_input datasets
    features_datasets = []
    predicton_datasets = []
    references = cube.attributes['reference_dataset']
    for dataset in cube.coord('dataset').points:
        if dataset in references:
            predicton_datasets.append(dataset)
            logger.info(
                "Using dataset '%s' as prediction_input for variable '%s' "
                "with index %d", dataset, cube.var_name,
                len(predicton_datasets) - 1)
        else:
            features_datasets.append(dataset)

    # Extract cubes
    feature_cube = cube.extract(iris.Constraint(dataset=features_datasets))
    prediction_cube = cube.extract(iris.Constraint(dataset=predicton_datasets))
    feature_cube.attributes['var_type'] = 'feature'
    prediction_cube.attributes['var_type'] = 'prediction_input'

    # Merge identical prediction_input if desired
    if merge_identical_pred_input:
        (_, unique_idx) = np.unique(prediction_cube.data, return_index=True)
        diff = len(prediction_cube.coord('dataset').points) - len(unique_idx)
        if diff > 0:
            prediction_cube = prediction_cube[unique_idx]
            logger.info(
                "Removed %d identical prediction_input points for variable "
                "'%s'", diff, prediction_cube.var_name)

    # Set new index for prediction input
    prediction_cube.coord('dataset').points = np.arange(
        len(prediction_cube.coord('dataset').points))
    return (feature_cube, prediction_cube)


def _cube_to_dataset_coord(cube):
    """Convert :class:`iris.cube.Cube` to :class:`iris.coords.AuxCoord`."""
    if cube.ndim == 1:
        datasets = cube.data
    elif cube.ndim == 2:
        cube.data = cube.data.astype(str, casting='same_kind')
        datasets = [''.join(d.compressed()) for d in cube.data]
    else:
        raise ValueError(
            f"Only 1D and 2D cubes supported, got {cube.ndim:d}D cube")
    return iris.coords.AuxCoord(datasets,
                                var_name='dataset',
                                long_name='dataset')


def _get_first_cube_with_coord(cubes, accepted_coord_names):
    """Load single cube of :class:`iris.cube.CubeList` with specific coords."""
    returned_cube = None
    returned_coord = None
    for cube in cubes:
        for coord_name in accepted_coord_names:
            try:
                coord = cube.coord(coord_name)
                returned_cube = cube
                returned_coord = coord
                break
            except iris.exceptions.CoordinateNotFoundError:
                pass
        if returned_cube is not None:
            break
    else:
        raise ValueError(
            f"No cube of {cubes} contains 'dataset' coordinate (i.e. one of "
            f"{accepted_coord_names})")
    return (returned_cube, returned_coord)


def _load_cube_with_dataset_coord(filename):
    """Load cube with single ``dataset``-like coordinate.

    Files created by NCL cannot be read using a simple
    :func:`iris.load_cube`.

    """
    cubes = iris.load(filename)
    accepted_coord_names = ('dataset', 'model')

    # Handle single cube
    if len(cubes) == 1:
        (cube, coord) = _get_first_cube_with_coord(cubes, accepted_coord_names)
        if cube.ndim != 1:
            raise ValueError(
                f"Only 1D cubes supported, got {cube.ndim:d}D cube in file "
                f"'{filename}'")
        coord.var_name = 'dataset'
        coord.standard_name = None
        coord.long_name = 'dataset'
        return cube

    # At most two cubes are supported
    if len(cubes) > 2:
        raise ValueError(
            f"Loading NCL file '{filename}' failed, at most 2 cubes are "
            f"supported, got {len(cubes):d}")

    # Get 'model' or 'dataset' cube
    dataset_cube = None
    for cube in cubes:
        if cube.var_name in accepted_coord_names:
            logger.debug("Found coordinate cube '%s'", cube.var_name)
            dataset_cube = cube
        else:
            data_cube = cube
    if dataset_cube is None:
        raise ValueError(
            f"No 'dataset' coordinate (one of {accepted_coord_names}) in "
            f"file '{filename}' available")

    # Create new coordinate
    if data_cube.ndim != 1:
        raise ValueError(
            f"Only 1D cubes supported, got {data_cube.ndim:d}D cube in file "
            f"'{filename}'")
    if data_cube.shape[0] != dataset_cube.shape[0]:
        raise ValueError(
            f"Got differing sizes for first dimension of data cube "
            f"({data_cube.shape[0]:d}) and dataset cube "
            f"({dataset_cube.shape[0]:d}) in file '{filename}'")
    aux_coord = _cube_to_dataset_coord(dataset_cube)
    data_cube.add_aux_coord(aux_coord, 0)
    return data_cube


def _create_scatterplot(x_data,
                        y_data,
                        numbers_as_markers=True,
                        plot_regression_line_mean=False,
                        axes=None,
                        **kwargs):
    """Create single scatterplot including regression line."""
    if axes is None:
        (_, axes) = plt.subplots()

    # Scatterplots
    scatter_kwargs = dict(kwargs)
    scatter_kwargs.pop('label', None)
    for (idx, _) in enumerate(x_data):
        if numbers_as_markers:
            axes.text(x_data[idx],
                      y_data[idx],
                      x_data.index.get_level_values(-1)[idx],
                      size=7,
                      **scatter_kwargs)
        else:
            axes.scatter(x_data[idx], y_data[idx], **scatter_kwargs)

    # Regression line
    line_kwargs = {**kwargs, 'linestyle': '-'}
    fill_between_kwargs = {**kwargs, 'alpha': 0.3}
    fill_between_kwargs.pop('label', None)
    if plot_regression_line_mean:
        mean_kwargs = {**kwargs, 'marker': 'o'}
        mean_kwargs.pop('label', None)
        mean_kwargs.pop('linestyle', None)
    else:
        mean_kwargs = None
    axes = _create_regplot(x_data,
                           y_data,
                           axes=axes,
                           line_kwargs=line_kwargs,
                           fill_between_kwargs=fill_between_kwargs,
                           mean_kwargs=mean_kwargs)
    return axes


def _create_pred_input_plot(x_pred,
                            x_pred_error,
                            axes,
                            vline_kwargs=None,
                            vspan_kwargs=None):
    """Create plot for prediction input data (vertical lines)."""
    if vline_kwargs is None:
        vline_kwargs = {'color': 'k', 'linestyle': ':', 'label': 'Observation'}
    if vspan_kwargs is None:
        vspan_kwargs = {'color': 'k', 'alpha': 0.1}
    axes.axvline(x_pred, **vline_kwargs)
    axes.axvspan(x_pred - x_pred_error, x_pred + x_pred_error, **vspan_kwargs)
    return axes


def _create_pred_output_plot(x_data,
                             y_data,
                             x_pred,
                             x_pred_error,
                             axes,
                             hline_kwargs=None):
    """Create plot for prediction input data (vertical lines)."""
    if hline_kwargs is None:
        hline_kwargs = {'color': 'k', 'linestyle': ':'}
    (_, y_mean, _) = get_constraint(x_data, y_data, x_pred, x_pred_error)
    axes.axhline(y_mean, **hline_kwargs)
    return axes


def _create_regplot(x_data,
                    y_data,
                    axes=None,
                    line_kwargs=None,
                    fill_between_kwargs=None,
                    mean_kwargs=None):
    """Create single regression line plot."""
    if axes is None:
        (_, axes) = plt.subplots()
    if line_kwargs is None:
        line_kwargs = {'linestyle': '-', 'label': 'Linear regression'}
    if fill_between_kwargs is None:
        fill_between_kwargs = {'alpha': 0.3}

    # Create regression line
    reg = regression_line(x_data, y_data)

    # Add R2 and p-value to label if possible
    text = rf"$R^2={reg['rvalue']**2:.2f}, p={reg['pvalue']:.4f}$"
    if 'label' in line_kwargs:
        line_kwargs['label'] += rf' ({text})'
    else:
        if reg['rvalue'] < 0.0:
            axes.text(0.62, 0.93, text, transform=axes.transAxes)
        else:
            axes.text(0.02, 0.93, text, transform=axes.transAxes)

    # Plots regression
    axes.plot(reg['x'], reg['y'], **line_kwargs)
    axes.fill_between(reg['x'], reg['y_minus_err'], reg['y_plus_err'],
                      **fill_between_kwargs)

    # Plot means if desired
    if mean_kwargs is not None:
        x_mean = np.mean(reg['x'])
        y_mean = np.mean(reg['y'])
        axes.scatter(x_mean, y_mean, **mean_kwargs)
    return axes


def _get_pandas_cube(pandas_object):
    """Convert :mod:`pandas` object to cube and fix coordinates."""
    cube = iris.pandas.as_cube(pandas_object)
    for coord_name in ('index', 'columns'):
        try:
            names = getattr(pandas_object, coord_name).names
        except AttributeError:
            continue
        coord = cube.coord(coord_name)
        if not np.issubdtype(coord.dtype, np.number):
            coord.points = coord.points.astype(str)
            if coord.bounds is not None:
                coord.bounds = coord.bounds.astype(str)
        names = [n for n in names if n is not None]
        if not names:
            continue
        new_coord_name = '-'.join(names)
        coord.var_name = new_coord_name
        coord.long_name = new_coord_name
    return cube


def _metadata_list_to_cube_list(metadata_list, source):
    """Convert :obj:`list` of :obj:`dict` to :class:`iris.cube.CubeList`."""
    cubes = iris.cube.CubeList()
    for metadata in metadata_list:
        for attr in ('data', 'dataset'):
            if attr not in metadata:
                raise AttributeError(
                    f"Entry {metadata} from source '{source}' does not "
                    f"contain necessary attribute '{attr}'")
        aux_coord = iris.coords.AuxCoord(metadata.pop('dataset'),
                                         var_name='dataset',
                                         long_name='dataset')
        data_of_cube = metadata.pop('data')
        if data_of_cube is None:
            data_of_cube = np.nan
        cube = iris.cube.Cube(data_of_cube,
                              aux_coords_and_dims=[(aux_coord, ())])
        for key in ('var_name', 'standard_name', 'long_name', 'units'):
            if key in metadata:
                setattr(cube, key, metadata.pop(key))
        cube.attributes = metadata
        cubes.append(cube)
    return cubes


def _gaussian_pdf(x_val, x_mean, x_std):
    """Return Gaussian probability density."""
    norm = np.sqrt(2.0 * np.pi * x_std**2)
    return np.exp(-(x_val - x_mean)**2 / 2.0 / x_std**2) / norm


def _get_target_pdf(x_data,
                    y_data,
                    obs_mean,
                    obs_std,
                    n_points=1000,
                    necessary_p_value=None):
    """Get PDF of target variable including linear regression information."""
    (x_data, y_data) = _check_x_y_arrays(x_data, y_data)
    spe = standard_prediction_error(x_data, y_data)
    reg = linregress(x_data, y_data)

    # Get evenly spaced range of y
    y_range = 1.5 * (np.max(y_data) - np.min(y_data))
    y_lin = np.linspace(
        np.min(y_data) - y_range,
        np.max(y_data) + y_range, n_points)

    # Use unconstrained value of desired and necessary
    if necessary_p_value is not None:
        if reg.pvalue > necessary_p_value:
            y_pdf = _gaussian_pdf(y_lin, np.mean(y_data), np.std(y_data))
            return (y_lin, y_pdf, reg)

    # Helper functions for calculation of constrained target variable
    def obs_pdf(x_new):
        """Return PDF of observations P(x)."""
        return _gaussian_pdf(x_new, obs_mean, obs_std)

    def cond_pdf(x_new, y_new):
        """Return conditional PDF P(y|x)."""
        y_pred = reg.slope * x_new + reg.intercept
        y_std = spe(x_new)
        return _gaussian_pdf(y_new, y_pred, y_std)

    def comb_pdf(x_new, y_new):
        """Return combined PDF P(y,x)."""
        return obs_pdf(x_new) * cond_pdf(x_new, y_new)

    # PDF of target variable P(y)
    x_range = 3 * obs_std
    y_pdf = [
        integrate.quad(comb_pdf,
                       obs_mean - x_range,
                       obs_mean + x_range,
                       args=(y, ))[0] for y in y_lin
    ]
    return (y_lin, np.array(y_pdf), reg)


def check_metadata(metadata, allowed_var_types=None):
    """Check metadata.

    Parameters
    ----------
    metadata : dict
        Metadata to check.
    allowed_var_types : list of str, optional
        Allowed var_types, defaults to ``ALLOWED_VAR_TYPES``.

    Raises
    ------
    KeyError
        Metadata does not contain necessary keys ``'var_type'`` and ``'tag'``.
    ValueError
        Got invalid value for key ``'var_type'``.

    """
    if allowed_var_types is None:
        allowed_var_types = ALLOWED_VAR_TYPES
    filename = metadata.get('filename', metadata)
    for key in ('var_type', 'tag'):
        if key not in metadata:
            raise KeyError(
                f"Necessary key '{key}' not given in metadata of file "
                f"'{filename}'")
    if metadata['var_type'] not in allowed_var_types:
        raise ValueError(
            f"Expected one of {allowed_var_types} for 'var_type' of file "
            f"'{filename}', got '{metadata['var_type']}'")


def get_input_files(cfg, patterns=None, ignore_patterns=None):
    """Get input files.

    Parameters
    ----------
    cfg : dict
        Recipe configuration.
    patterns : list of str, optional
        Use only ancestor files that match these patterns as input files.
    ignore_patterns : list of str, optional
        Ignore input files that match these patterns.

    Returns
    -------
    list of str
        Input files.

    """
    input_files = []

    # Include only files that match patterns
    if patterns is None:
        patterns = []
    if not patterns:
        patterns.append('*.nc')
    for pattern in patterns:
        logger.debug("Looking for files matching the pattern '%s'", pattern)
        input_files.extend(io.get_all_ancestor_files(cfg, pattern=pattern))

    # Ignore files
    if not ignore_patterns:
        return input_files
    ignore_files = []
    for pattern in ignore_patterns:
        logger.debug("Ignoring for files matching the pattern '%s'", pattern)
        ignore_files.extend(io.get_all_ancestor_files(cfg, pattern=pattern))
    valid_files = []
    for filename in input_files:
        if filename not in ignore_files:
            valid_files.append(filename)
    return valid_files


def get_xy_data_without_nans(data_frame, feature, label):
    """Get (X, Y) data for ``(feature, label)`` combination without nans.

    Parameters
    ----------
    data_frame : pandas.DataFrame
        Training data.
    feature : str
        Name of the feature data.
    label : str
        Name of the label data.

    Returns
    -------
    tuple
        Tuple containing a :class:`pandas.DataFrame` for the X axis (feature)
        and a :class:`pandas.DataFrame` for the Y axis (label) without
        missing values.

    """
    idx_slice = pd.IndexSlice[:, [feature, label]]
    data_frame_xy = data_frame.loc[:, idx_slice]
    data_frame_xy.columns = data_frame_xy.columns.droplevel()
    data_frame_xy = data_frame_xy.dropna()
    x_data = data_frame_xy[feature]
    y_data = data_frame_xy[label]
    return (x_data, y_data)


def get_input_data(cfg):
    """Extract input data.

    Return training data, prediction input data and corresponding attributes.

    Parameters
    ----------
    cfg : dict
        Recipe configuration.

    Returns
    -------
    tuple
        A tuple containing the training data (:class:`pandas.DataFrame`), the
        prediction input data (:class:`pandas.DataFrame`) and the
        corresponding attributes (:obj:`dict`).

    """
    input_files = get_input_files(cfg,
                                  patterns=cfg.get('patterns'),
                                  ignore_patterns=cfg.get('ignore_patterns'))
    logger.debug("Found files:\n%s", pformat(input_files))

    # Get cubes
    external_file = _get_external_file(cfg.get('read_external_file'),
                                       cfg['auxiliary_data_dir'])
    cubes = _get_cube_list(
        input_files,
        recipe=cfg['recipe'],
        additional_data=cfg.get('additional_data'),
        external_file=external_file,
        merge_identical_pred_input=cfg.get('merge_identical_pred_input', True),
    )

    # Extract attributes for features and labels
    attributes = _get_attributes(cubes)

    # Extract DataFrames
    label_all_data = cfg.get('all_data_label', 'all')
    group_by = cfg.get('group_by')
    if group_by:
        logger.info("Grouping features and labels by '%s'", group_by)
    else:
        logger.info("Using label '%s' to label data in plots", label_all_data)
    features = _get_data_frame('feature', cubes, label_all_data, group_by)
    label = _get_data_frame('label', cubes, label_all_data, group_by)
    pred_input = _get_data_frame('prediction_input', cubes, label_all_data,
                                 group_by)
    pred_input_err = _get_data_frame('prediction_input_error', cubes,
                                     label_all_data, group_by)

    # Unify indices of features and label
    for row in features.index.difference(label.index):
        label = label.append(pd.Series(name=row))
    for row in label.index.difference(features.index):
        features = features.append(pd.Series(name=row))

    # Sort data frames
    for data_frame in (features, label, pred_input, pred_input_err):
        data_frame.sort_index(axis=0, inplace=True)
        data_frame.sort_index(axis=1, inplace=True)

    # Check data
    (features, label, pred_input,
     pred_input_err) = _check_data_frames(features, label, pred_input,
                                          pred_input_err)
    training_data = pd.concat([features, label], axis=1, keys=['x', 'y'])
    training_data['idx'] = np.arange(len(training_data.index)) + 1
    training_data.set_index('idx', append=True, inplace=True)
    training_data.index.names = [group_by, 'dataset', 'idx']
    prediction_data = pd.concat([pred_input, pred_input_err],
                                axis=1,
                                keys=['mean', 'error'])
    if training_data.dropna().shape[0] < 2:
        logger.error("Invalid training data:\n%s", training_data)
        raise ValueError(
            f"Expected at least 2 independent observations (=climate models) "
            f"where all training data (features and target label) is "
            f"available, got {training_data.dropna().shape[0]:d}")

    # Logger output
    with pd.option_context(*PANDAS_PRINT_OPTIONS):
        logger.info("Found training data:\n%s", training_data)
        logger.info("Found prediction data:\n%s", prediction_data)
    return (training_data, prediction_data, attributes)


def combine_groups(groups):
    """Combine :obj:`list` of groups to a single :obj:`str`.

    Parameters
    ----------
    groups : list of str
        List of group names.

    Returns
    -------
    str
        Combined :obj:`str`.

    """
    new_str = ', '.join(groups)
    return new_str


def pandas_object_to_cube(pandas_object,
                          index_droplevel=None,
                          columns_droplevel=None,
                          **kwargs):
    """Convert pandas object to :class:`iris.cube.Cube`.

    Parameters
    ----------
    pandas_object : pandas.DataFrame or pandas.Series
        Data to convert.
    index_droplevel : int or list of int, optional
        Drop levels of index if not ``None``.
    columns_droplevel : int or list of int, optional
        Drop levels of columns if not ``None``. Can only be used if
        ``pandas_object`` is a :class:`pandas.DataFrame`.
    **kwargs : Keyword arguments
        Keyword arguments used for the cube metadata, e.g. ``standard_name``,
        ``var_name``, etc.

    Returns
    -------
    iris.cube.Cube
        Data cube.

    Raises
    ------
    TypeError
        ``columns_droplevel`` is used when ``pandas_object`` is not a
        :class:`pandas.DataFrame`.

    """
    pandas_object = pandas_object.copy()
    if index_droplevel is not None:
        pandas_object.index = pandas_object.index.droplevel(index_droplevel)
    if columns_droplevel is not None:
        try:
            pandas_object.columns = pandas_object.columns.droplevel(
                columns_droplevel)
        except AttributeError:
            raise TypeError(
                f"'columns_droplevel' only supported for pandas.DataFrame "
                f"object, got {type(pandas_object)}")
    cube = _get_pandas_cube(pandas_object)
    for (key, val) in kwargs.items():
        setattr(cube, key, val)
    return cube


def set_plot_appearance(axes, attributes, **kwargs):
    """Set appearance of a plot.

    Parameters
    ----------
    axes : matplotlib.axes.Axes
        Matplotlib Axes object which contains the plot.
    attributes : dict
        Plot attributes.
    **kwargs : Keyword arguments
        Keyword arguments of the form ``plot_option=tag`` where ``plot_option``
        is something like ``plot_title``, ``plot_xlabel``, ``plot_xlim``, etc.
        and ``tag`` a key for the plot attributes :obj:`dict` that describes
        which attributes should be considered for that ``plot_option``.

    """
    for (plot_option, tag) in kwargs.items():
        plot_func = plot_option.replace('plot_', 'set_')
        value = attributes[tag][plot_option]
        getattr(axes, plot_func)(value)


def get_caption(attributes, feature, label, group=None):
    """Construct caption from plotting attributes for (feature, label) pair.

    Parameters
    ----------
    attributes : dict
        Plot attributes.
    feature : str
        Feature.
    label : str
        Label.
    group : str, optional
        Group.

    Returns
    -------
    str
        Caption.

    Raises
    ------
    KeyError
        ``attributes`` does not include necessary keys.

    """
    group_str = '' if group is None else f' ({group})'
    if feature not in attributes:
        raise KeyError(
            f"Attributes do not include necessary key for feature '{feature}'")
    if label not in attributes:
        raise KeyError(
            f"Attributes do not include necessary key for label '{label}'")
    feature_attrs = attributes[feature]
    label_attrs = attributes[label]
    if 'plot_title' not in feature_attrs:
        raise KeyError(
            f"Attributes for feature '{feature}' does not include necessary "
            f"key 'plot_title'")
    if 'plot_xlabel' not in feature_attrs:
        raise KeyError(
            f"Attributes for feature '{feature}' does not include necessary "
            f"key 'plot_xlabel'")
    if 'plot_ylabel' not in label_attrs:
        raise KeyError(
            f"Attributes for label '{label}' does not include necessary "
            f"key 'plot_ylabel'")
    caption = (f"{attributes[feature]['plot_title']}: "
               f"{attributes[label]['plot_ylabel']} vs. "
               f"{attributes[feature]['plot_xlabel']}{group_str}.")
    return caption


def get_provenance_record(attributes, tags, **kwargs):
    """Get provenance record.

    Parameters
    ----------
    attributes : dict
        Plot attributes. All provenance keys need to start with
        ``'provenance_'``.
    tags : list of str
        Tags used to retrieve data from the ``attributes`` :obj:`dict`, i.e.
        features and/or label.
    **kwargs : Keyword arguments
        Additional ``key:value`` pairs directly passed to the provenance record
        :obj:`dict`. All values may include the format strings ``{feature}``
        and ``{label}``.

    Returns
    -------
    dict
        Provenance record.

    """
    record = {}
    for tag in tags:
        for (key, value) in attributes[tag].items():
            if key.startswith('provenance_'):
                key = key.replace('provenance_', '')
                record.setdefault(key, [])
                if isinstance(value, str):
                    record[key].append(value)
                else:
                    record[key].extend(value)
            record.setdefault('ancestors', [])
            if key == 'filename':
                if isinstance(value, str):
                    record['ancestors'].append(value)
                else:
                    record['ancestors'].extend(value)
    for (key, value) in record.items():
        if isinstance(value, list):
            record[key] = list(set(value))
    record.update(kwargs)
    return record


def get_colors(cfg, groups=None):
    """Get color palette.

    Parameters
    ----------
    cfg : dict
        Recipe configuration.
    groups : list, optional
        Use to check whether color for combining groups has to be added.

    Returns
    -------
    list
        List of colors that can be used for :mod:`matplotlib`.

    """
    palette = cfg.get('seaborn_settings', {}).get('palette')
    colors = sns.color_palette(palette=palette)
    if groups is None:
        return colors
    if len(groups) > 1 and cfg.get('combine_groups', False):
        return [COLOR_COMBINED_GROUPS] + colors
    return colors


def get_groups(training_data, add_combined_group=False):
    """Extract groups from training data.

    Parameters
    ----------
    training_data : pandas.DataFrame
        Training data (features, label).
    add_combined_group : bool, optional (default: False)
        Add combined group of all other groups at the beginning of the
        returned :obj:`list`.

    Returns
    -------
    list of str
        Groups.

    """
    groups = list(set(training_data.index.get_level_values(0)))
    groups.sort()
    if add_combined_group and len(groups) > 1:
        groups.insert(0, combine_groups(groups))
    return groups


def plot_individual_scatterplots(training_data, pred_input_data, attributes,
                                 basename, cfg):
    """Plot individual scatterplots for the different groups.

    Plot scatterplots for all pairs of ``(feature, label)`` data (Separate plot
    for each group).

    Parameters
    ----------
    training_data : pandas.DataFrame
        Training data (features, label).
    pred_input_data : pandas.DataFrame
        Prediction input data (mean and error).
    attributes : dict
        Plot attributes for the different features and the label data.
    basename : str
        Basename for the name of the file.
    cfg : dict
        Recipe configuration.

    """
    logger.info("Plotting individual scatterplots")
    label = training_data.y.columns[0]
    groups = get_groups(training_data,
                        add_combined_group=cfg.get('combine_groups', False))

    # Iterate over features
    for feature in training_data.x.columns:
        (x_data, y_data) = get_xy_data_without_nans(training_data, feature,
                                                    label)

        # Individual plots
        colors = get_colors(cfg, groups=groups)
        for (idx, group) in enumerate(groups):
            try:
                x_sub_data = x_data.loc[group]
                y_sub_data = y_data.loc[group]
                index_droplevel = 1
            except KeyError:
                x_sub_data = x_data
                y_sub_data = y_data
                index_droplevel = [0, 2]
            axes = _create_scatterplot(
                x_sub_data,
                y_sub_data,
                numbers_as_markers=cfg.get('numbers_as_markers', False),
                plot_regression_line_mean=cfg.get('plot_regression_line_mean',
                                                  False),
                color=colors[idx],
                label=group)
            axes = _create_pred_input_plot(
                pred_input_data['mean'][feature].values,
                pred_input_data['error'][feature].values, axes)
            axes = _create_pred_output_plot(
                x_sub_data,
                y_sub_data,
                pred_input_data['mean'][feature].values,
                pred_input_data['error'][feature].values,
                axes,
                hline_kwargs={
                    'color': colors[idx],
                    'linestyle': ':'
                },
            )
            set_plot_appearance(axes,
                                attributes,
                                plot_title=feature,
                                plot_xlabel=feature,
                                plot_ylabel=label,
                                plot_xlim=feature,
                                plot_ylim=label)
            plt.legend(**LEGEND_KWARGS)
            filename = (f"scatterplot_{basename}_{feature}_"
                        f"{group.replace(', ', '-')}")
            plot_path = get_plot_filename(filename, cfg)
            plt.savefig(plot_path,
                        **cfg.get('savefig_kwargs', {}))
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Provenance
            cubes = iris.cube.CubeList([
                pandas_object_to_cube(
                    x_sub_data,
                    index_droplevel=index_droplevel,
                    var_name=feature,
                    long_name=attributes[feature]['plot_xlabel'],
                    units=attributes[feature]['units']),
                pandas_object_to_cube(
                    y_sub_data,
                    index_droplevel=index_droplevel,
                    var_name=label,
                    long_name=attributes[label]['plot_ylabel'],
                    units=attributes[label]['units']),
            ])
            netcdf_path = get_diagnostic_filename(filename, cfg)
            io.iris_save(cubes, netcdf_path)
            provenance_record = get_provenance_record(
                attributes, [feature, label],
                caption=get_caption(attributes, feature, label, group=group),
                plot_type='scatter',
                plot_file=plot_path)
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(netcdf_path, provenance_record)


def plot_merged_scatterplots(training_data, pred_input_data, attributes,
                             basename, cfg):
    """Plot merged scatterplots (all groups in one plot).

    Plot scatterplots for all pairs of ``(feature, label)`` data (all groups in
    one plot).

    Parameters
    ----------
    training_data : pandas.DataFrame
        Training data (features, label).
    pred_input_data : pandas.DataFrame
        Prediction input data (mean and error).
    attributes : dict
        Plot attributes for the different features and the label data.
    basename : str
        Basename for the name of the file.
    cfg : dict
        Recipe configuration.

    """
    logger.info("Plotting merged scatterplots")
    label = training_data.y.columns[0]
    groups = get_groups(training_data,
                        add_combined_group=cfg.get('combine_groups', False))
    numbers_as_markers = cfg.get('numbers_as_markers', False)
    plot_regression_line_mean = cfg.get('plot_regression_line_mean', False)
    colors = get_colors(cfg)

    # Iterate over features
    for feature in training_data.x.columns:
        (x_data, y_data) = get_xy_data_without_nans(training_data, feature,
                                                    label)
        (_, axes) = plt.subplots()
        if len(groups) > 1 and cfg.get('combine_groups', False):
            axes = _create_regplot(
                x_data,
                y_data,
                axes=axes,
                line_kwargs={
                    'color': COLOR_COMBINED_GROUPS,
                    'label': groups[0],
                    'linestyle': '-'
                },
                fill_between_kwargs={
                    'color': COLOR_COMBINED_GROUPS,
                    'alpha': 0.3
                },
                mean_kwargs=(None
                             if not cfg.get('plot_regression_line_mean') else {
                                 'color': COLOR_COMBINED_GROUPS,
                                 'marker': 'o'
                             }),
            )
            axes = _create_pred_output_plot(
                x_data,
                y_data,
                pred_input_data['mean'][feature].values,
                pred_input_data['error'][feature].values,
                axes,
                hline_kwargs={
                    'color': COLOR_COMBINED_GROUPS,
                    'linestyle': ':'
                },
            )
            for (idx, group) in enumerate(groups[1:]):
                axes = _create_scatterplot(
                    x_data.loc[group],
                    y_data.loc[group],
                    numbers_as_markers=numbers_as_markers,
                    plot_regression_line_mean=plot_regression_line_mean,
                    axes=axes,
                    color=colors[idx],
                    label=group,
                )
                axes = _create_pred_output_plot(
                    x_data.loc[group],
                    y_data.loc[group],
                    pred_input_data['mean'][feature].values,
                    pred_input_data['error'][feature].values,
                    axes,
                    hline_kwargs={
                        'color': colors[idx],
                        'linestyle': ':'
                    },
                )
        else:
            for (idx, group) in enumerate(groups):
                axes = _create_scatterplot(
                    x_data.loc[group],
                    y_data.loc[group],
                    numbers_as_markers=numbers_as_markers,
                    plot_regression_line_mean=plot_regression_line_mean,
                    axes=axes,
                    color=colors[idx],
                    label=group,
                )
                axes = _create_pred_output_plot(
                    x_data.loc[group],
                    y_data.loc[group],
                    pred_input_data['mean'][feature].values,
                    pred_input_data['error'][feature].values,
                    axes,
                    hline_kwargs={
                        'color': colors[idx],
                        'linestyle': ':'
                    },
                )
        axes = _create_pred_input_plot(
            pred_input_data['mean'][feature].values,
            pred_input_data['error'][feature].values, axes)
        set_plot_appearance(axes,
                            attributes,
                            plot_title=feature,
                            plot_xlabel=feature,
                            plot_ylabel=label,
                            plot_xlim=feature,
                            plot_ylim=label)
        plt.legend(**LEGEND_KWARGS)
        filename = f'scatterplot_merged_{basename}_{feature}'
        plot_path = get_plot_filename(filename, cfg)
        plt.savefig(plot_path,
                    **cfg.get('savefig_kwargs', {}))
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Provenance
        cubes = iris.cube.CubeList([
            pandas_object_to_cube(x_data,
                                  index_droplevel=[0, 2],
                                  var_name=feature,
                                  long_name=attributes[feature]['plot_xlabel'],
                                  units=attributes[feature]['units']),
            pandas_object_to_cube(y_data,
                                  index_droplevel=[0, 2],
                                  var_name=label,
                                  long_name=attributes[label]['plot_ylabel'],
                                  units=attributes[label]['units']),
        ])
        netcdf_path = get_diagnostic_filename(filename, cfg)
        io.iris_save(cubes, netcdf_path)
        provenance_record = get_provenance_record(attributes, [feature, label],
                                                  caption=get_caption(
                                                      attributes, feature,
                                                      label),
                                                  plot_type='scatter',
                                                  plot_file=plot_path)
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(netcdf_path, provenance_record)


def create_simple_scatterplot(x_data, y_data, obs_mean, obs_std):
    """Create simple scatterplot of an emergent relationship (without saving).

    Parameters
    ----------
    x_data : numpy.ndarray
        X data of the emergent constraint.
    y_data : numpy.ndarray
        Y data of the emergent constraint.
    obs_mean : float
        Mean of observational data.
    obs_std : float
        Standard deviation of observational data.

    """
    logger.debug("Plotting simple scatterplot")
    (fig, axes) = plt.subplots()
    axes.scatter(x_data, y_data, color='k', marker='o')
    line_kwargs = {'color': 'C1', 'linestyle': '-'}
    fill_between_kwargs = {**line_kwargs, 'alpha': 0.3}
    axes = _create_regplot(x_data, y_data, axes=axes, line_kwargs=line_kwargs,
                           fill_between_kwargs=fill_between_kwargs)
    axes = _create_pred_input_plot(obs_mean, obs_std, axes)
    return (fig, axes)


def plot_target_distributions(training_data, pred_input_data, attributes,
                              basename, cfg):
    """Plot distributions of target variable for every feature.

    Parameters
    ----------
    training_data : pandas.DataFrame
        Training data (features, label).
    pred_input_data : pandas.DataFrame
        Prediction input data (mean and error).
    attributes : dict
        Plot attributes for the different features and the label data.
    basename : str
        Basename for the name of the file.
    cfg : dict
        Recipe configuration.

    """
    logger.info("Plotting distributions of target variable")
    label = training_data.y.columns[0]
    groups = get_groups(training_data,
                        add_combined_group=cfg['combine_groups'])
    summary_columns = pd.MultiIndex.from_product(
        [groups, ['best estimate', 'range', 'min', 'max']])
    summary = pd.DataFrame(columns=summary_columns)

    # Iterate over features
    for feature in training_data.x.columns:
        (x_data, y_data) = get_xy_data_without_nans(training_data, feature,
                                                    label)
        colors = get_colors(cfg, groups=groups)
        summary_for_feature = pd.Series(index=summary_columns, name=feature)

        # Iterate over groups
        for (idx, group) in enumerate(groups):
            try:
                x_sub_data = x_data.loc[group]
                y_sub_data = y_data.loc[group]
            except KeyError:
                x_sub_data = x_data
                y_sub_data = y_data
            (y_lin, y_pdf) = target_pdf(
                x_sub_data,
                y_sub_data,
                pred_input_data['mean'][feature].values,
                pred_input_data['error'][feature].values,
            )

            # Plots
            axes = sns.histplot(y_sub_data,
                                bins=7,
                                stat='density',
                                color=colors[idx],
                                alpha=0.4)
            axes.plot(y_lin,
                      y_pdf,
                      color=colors[idx],
                      linestyle='-',
                      label=group)

            # Print results
            (y_min, y_mean, y_max) = get_constraint(
                x_sub_data,
                y_sub_data,
                pred_input_data['mean'][feature].values,
                pred_input_data['error'][feature].values,
                confidence_level=cfg['confidence_level'],
            )
            y_error = np.max([y_max - y_mean, y_mean - y_min])
            reg = linregress(x_sub_data.values, y_sub_data.values)
            logger.info(
                "Constrained %s for feature '%s' and group '%s': %.2f ± %.2f "
                "(%i%% confidence level), R2 = %f, p = %f", label,
                feature, group, y_mean, y_error,
                int(100.0 * cfg['confidence_level']), reg.rvalue**2,
                reg.pvalue)

            # Save results of group
            summary_for_feature[(group, 'best estimate')] = y_mean
            summary_for_feature[(group, 'range')] = y_max - y_min
            summary_for_feature[(group, 'min')] = y_min
            summary_for_feature[(group, 'max')] = y_max

        # Save results to feature
        summary = summary.append(summary_for_feature)

        # Plot appearance
        set_plot_appearance(axes, attributes, plot_title=feature)
        axes.set_xlabel(attributes[label]['plot_ylabel'])
        axes.set_ylabel('Probability density')
        if attributes[label]['plot_ylim'] is not None:
            axes.set_xlim(attributes[label]['plot_ylim'])
        axes.set_ylim([0.0, 1.0])
        plt.legend(loc='best')

        # Save plot
        plot_path = get_plot_filename(
            f'target_distribution_{basename}_{feature}', cfg)
        plt.savefig(plot_path,
                    **cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

    # Print mean results
    with pd.option_context(*PANDAS_PRINT_OPTIONS):
        logger.info("Constrained ranges:\n%s", summary)
        summary = summary.mean(axis=0)
        logger.info("Mean of constrained ranges:\n%s", summary)


def export_csv(data_frame, attributes, basename, cfg, tags=None):
    """Export CSV file.

    Parameters
    ----------
    data_frame : pandas.DataFrame
        Data to export.
    attributes : dict
        Plot attributes for the different features and the label data. Used to
        retrieve provenance information.
    basename : str
        Basename for the name of the file.
    cfg : dict
        Recipe configuration.
    tags : iterable of str, optional
        Tags for which provenance information should be retrieved (using
        ``attributes``). If not specified, use (last level of) columns of the
        given ``data_frame``.

    Returns
    -------
    str
        Path to the new CSV file.

    """
    logger.info("Exporting CSV file for '%s'", basename)
    csv_path = get_diagnostic_filename(basename, cfg).replace('.nc', '.csv')
    data_frame.to_csv(csv_path)
    logger.info("Wrote %s", csv_path)
    if tags is None:
        tags = data_frame.columns.get_level_values(-1)
    provenance_record = get_provenance_record(attributes,
                                              tags,
                                              caption=basename)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(csv_path, provenance_record)
    return csv_path


def standard_prediction_error(x_data, y_data):
    """Return a function to calculate standard prediction error.

    The standard prediction error of a linear regression is the error when
    predicting a data point which was not used to fit the regression line in
    the first place.

    Parameters
    ----------
    x_data : numpy.ndarray
        X data used to fit the linear regression.
    y_data : numpy.ndarray
        Y data used to fit the linear regression.

    Returns
    -------
    callable
        Function that takes a :obj:`float` as single argument (representing the
        X value of a new data point) and returns the standard prediction error
        for that.

    """
    (x_data, y_data) = _check_x_y_arrays(x_data, y_data)
    reg = linregress(x_data, y_data)
    y_pred = reg.slope * x_data + reg.intercept
    n_data = x_data.shape[0]
    see = np.sqrt(np.sum(np.square(y_data - y_pred)) / (n_data - 2))
    x_mean = np.mean(x_data)
    ssx = np.sum(np.square(x_data - x_mean))

    def spe(x_new):
        """Return standard prediction error."""
        return see * np.sqrt(1.0 + 1.0 / n_data + (x_new - x_mean)**2 / ssx)

    return spe


def regression_line(x_data, y_data, n_points=1000):
    """Return x and y coordinates of the regression line (mean and error).

    Parameters
    ----------
    x_data : numpy.ndarray
        X data used to fit the linear regression.
    y_data : numpy.ndarray
        Y data used to fit the linear regression.
    n_points : int, optional (default: 1000)
        Number of points for the regression lines.

    Returns
    -------
    dict
        :class:`numpy.ndarray` s for the keys ``'x'``, ``'y'``,
        ``'y_minus_err'``, ``'y_plus_err'``, ``'slope'``, ``'intercept'``,
        ``'pvalue'`` and ``'rvalue'``.

    """
    (x_data, y_data) = _check_x_y_arrays(x_data, y_data)
    spe = np.vectorize(standard_prediction_error(x_data, y_data))
    out = {}
    reg = linregress(x_data, y_data)
    x_range = np.max(x_data) - np.min(x_data)
    x_lin = np.linspace(
        np.min(x_data) - x_range,
        np.max(x_data) + x_range, n_points)
    out['x'] = x_lin
    out['y'] = reg.slope * x_lin + reg.intercept
    out['y_minus_err'] = out['y'] - spe(x_lin)
    out['y_plus_err'] = out['y'] + spe(x_lin)
    out['slope'] = reg.slope
    out['intercept'] = reg.intercept
    out['pvalue'] = reg.pvalue
    out['rvalue'] = reg.rvalue
    return out


def target_pdf(x_data,
               y_data,
               obs_mean,
               obs_std,
               n_points=1000,
               necessary_p_value=None):
    """Calculate probability density function (PDF) for target variable.

    Parameters
    ----------
    x_data : numpy.ndarray
        X data of the emergent constraint.
    y_data : numpy.ndarray
        Y data of the emergent constraint.
    obs_mean : float
        Mean of observational data.
    obs_std : float
        Standard deviation of observational data.
    n_points : int, optional (default: 1000)
        Number of sampled points for PDF of target variable.
    necessary_p_value : float, optional
        If given, return unconstrained PDF (using Gaussian distribution with
        unconstrained mean and standard deviation) when `p`-value of emergent
        relationship is greater than the given necessary `p`-value.

    Returns
    -------
    tuple of numpy.ndarray
        x and y values for the PDF.

    """
    (y_lin, y_pdf, _) = _get_target_pdf(x_data,
                                        y_data,
                                        obs_mean,
                                        obs_std,
                                        n_points=n_points,
                                        necessary_p_value=necessary_p_value)
    return (y_lin, y_pdf)


def cdf(data, pdf):
    """Calculate cumulative distribution function for a 1-dimensional PDF.

    Parameters
    ----------
    data : numpy.ndarray
        Data points (1D array).
    pdf : numpy.ndarray
        Corresponding probability density function (PDF).

    Returns
    -------
    numpy.ndarray
        Corresponding cumulative distribution function (CDF).

    """
    idx_range = range(1, len(data) + 1)
    cum_dens = [integrate.simps(pdf[:idx], data[:idx]) for idx in idx_range]
    return np.array(cum_dens)


def constraint_info_array(x_data,
                          y_data,
                          obs_mean,
                          obs_std,
                          n_points=1000,
                          necessary_p_value=None):
    """Get array with all relevant parameters of emergent constraint.

    Parameters
    ----------
    x_data : numpy.ndarray
        X data of the emergent constraint.
    y_data : numpy.ndarray
        Y data of the emergent constraint.
    obs_mean : float
        Mean of observational data.
    obs_std : float
        Standard deviation of observational data.
    n_points : int, optional (default: 1000)
        Number of sampled points for PDF of target variable.
    necessary_p_value : float, optional
        If given, replace constrained mean and standard deviation with
        unconstrained values when `p`-value of emergent relationship is greater
        than the given necessary `p`-value.

    Returns
    -------
    numpy.ndarray
        Array of shape (8,) with the elements:
            0. Constrained mean of target variable.
            1. Constrained standard deviation of target variable.
            2. Unconstrained mean of target variable.
            3. Unconstrained standard deviation of target variable.
            4. Slope of emergent relationship.
            5. Intercept of emergent relationship.
            6. Correlation coefficient `r` of emergent relationship.
            7. `p`-value of emergent relationship.

    """
    (y_lin, y_pdf, reg) = _get_target_pdf(x_data,
                                          y_data,
                                          obs_mean,
                                          obs_std,
                                          n_points=n_points,
                                          necessary_p_value=necessary_p_value)
    norm = np.sum(y_pdf)
    y_mean = np.sum(y_lin * y_pdf) / norm
    y_std = np.sqrt(np.sum((y_lin - y_mean)**2 * y_pdf) / norm)
    info = [
        y_mean, y_std,
        np.ma.mean(y_data),
        np.ma.std(y_data), reg.slope, reg.intercept, reg.rvalue, reg.pvalue
    ]
    return np.array(info)


def get_constraint(x_data, y_data, obs_mean, obs_std, confidence_level=0.66):
    """Get constraint on target variable.

    Parameters
    ----------
    x_data : numpy.ndarray
        X data of the emergent constraint.
    y_data : numpy.ndarray
        Y data of the emergent constraint.
    obs_mean : float
        Mean of observational data.
    obs_std : float
        Standard deviation of observational data.
    confidence_level : float, optional (default: 0.66)
        Confindence level to estimate the range of the target variable.

    Returns
    -------
    tuple of float
        Lower confidence limit, best estimate and upper confidence limit of
        target variable.

    """
    (x_data, y_data) = _check_x_y_arrays(x_data, y_data)
    (y_lin, y_pdf) = target_pdf(x_data, y_data, obs_mean, obs_std)
    y_mean = np.sum(y_lin * y_pdf) / np.sum(y_pdf)
    y_cdf = cdf(y_lin, y_pdf)
    y_index_range = np.nonzero((y_cdf >= (1.0 - confidence_level) / 2.0)
                               & (y_cdf <= (1.0 + confidence_level) / 2.0))
    y_range = y_lin[y_index_range]
    return (np.min(y_range), y_mean, np.max(y_range))


def get_constraint_from_df(training_data,
                           pred_input_data,
                           confidence_level=0.66):
    """Get constraint on target variable from :class:`pandas.DataFrame`.

    Parameters
    ----------
    training_data : pandas.DataFrame
        Training data (features, label).
    pred_input_data : pandas.DataFrame
        Prediction input data (mean and error).
    confidence_level : float, optional (default: 0.66)
        Confindence level to estimate the range of the target variable.

    Returns
    -------
    tuple of float
        Lower confidence limit, best estimate and upper confidence limit of
        target variable.

    """
    if len(training_data.columns) != 2:
        raise ValueError(
            f"Expected exactly two columns for training data (feature and "
            f"label), got {len(training_data.columns):d}")
    if len(pred_input_data.columns) != 2:
        raise ValueError(
            f"Expected exactly two columns for prediction input data (mean "
            f"and error, got {len(pred_input_data.columns):d}")

    # Extract data
    label = training_data.y.columns[0]
    feature = training_data.x.columns[0]
    (x_data, y_data) = get_xy_data_without_nans(training_data, feature, label)
    x_pred = pred_input_data['mean'][feature].values[0]
    x_pred_error = pred_input_data['error'][feature].values[0]

    # Calculate constraint
    constraint = get_constraint(x_data,
                                y_data,
                                x_pred,
                                x_pred_error,
                                confidence_level=confidence_level)
    return constraint
