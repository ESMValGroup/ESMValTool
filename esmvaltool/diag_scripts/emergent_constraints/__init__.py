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
from scipy.stats import linregress, multivariate_normal
from sklearn.linear_model import LinearRegression

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, io)

logger = logging.getLogger(__name__)


COLORS = sns.color_palette()
COLORS_ALL = ['gray'] + COLORS[1:]
LEGEND_KWARGS = {
    'loc': 'center left',
    'bbox_to_anchor': [1.05, 0.5],
    'borderaxespad': 0.0,
}
PANDAS_PRINT_OPTIONS = ['display.max_rows', None, 'display.max_colwidth', -1]


def _check_feature_array(x_array, single_sample=False):
    """Check X array."""
    x_array = np.array(x_array)
    if x_array.ndim == 0:
        x_array = x_array.reshape(1, 1)
    elif x_array.ndim == 1:
        if single_sample:
            x_array = x_array.reshape(1, -1)
        else:
            x_array = x_array.reshape(-1, 1)
    elif x_array.ndim > 2:
        raise ValueError(
            f"Expected at most 2D array for X, got {x_array.ndim:d}D array")
    if single_sample:
        if x_array.shape[0] > 1:
            raise ValueError(
                f"Expected single sample for X (i.e. one observation), got "
                f"{x_array.shape[0]:d}")
    return x_array


def _check_training_arrays(x_array, y_array=None):
    """Ensure that the training input arrays have correct shapes."""
    x_array = _check_feature_array(x_array, single_sample=False)
    if y_array is None:
        return x_array
    y_array = np.array(y_array)
    if y_array.ndim != 1:
        raise ValueError(
            f"Expected 1D array for Y, got {y_array.ndim:d}D array")
    if x_array.shape[0] != y_array.shape[0]:
        raise ValueError(
            f"Expected X and Y arrays with identical first dimension, got "
            f"shapes {x_array.shape} and {y_array.shape}")
    return (x_array, y_array)


def _check_prediction_arrays(x_new, x_train=None, x_cov=None):
    """Ensure that the prediction input arrays have correct shapes."""
    x_new = _check_feature_array(x_new, single_sample=True)
    if x_train is not None:
        _check_feature_array(x_train, single_sample=False)
        if x_new.shape[1] != x_train.shape[1]:
            raise ValueError(
                f"Expected identical number of predictors for training and "
                f"prediction data, got {x_train.shape[1]:d} and "
                f"{x_new.shape[1]:d}, respectively")
    if x_cov is None:
        return x_new
    x_cov = np.array(x_cov)
    if x_cov.ndim == 0:
        x_cov = x_cov.reshape(1, 1)
    elif x_cov.ndim != 2:
        raise ValueError(
            f"Expected 0D or 2D array for covariance matrix of observations, "
            f"got {x_cov.ndim:d}D array")
    if x_cov.shape != (x_new.shape[1], x_new.shape[1]):
        raise ValueError(
            f"Expected quadratic array for covariance matrix of observations "
            f"with shape {(x_new.shape[1], x_new.shape[1])} (number of "
            f"predictors), got {x_cov.shape}")
    return (x_new, x_cov)


def _get_x_ranges(obs_mean, obs_cov):
    """Get integration limits (8 sigma interval, includes > 99.99% of area)."""
    x_ranges = []
    for (idx, var) in enumerate(np.diag(obs_cov)):
        std = np.sqrt(var)
        x_range = [obs_mean[0][idx] - 4.0 * std, obs_mean[0][idx] + 4.0 * std]
        x_ranges.append(x_range)
    return x_ranges


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
    diff_not_in_data_frame = list(ref_data_frame.columns.difference(
        data_frame.columns))
    if diff_not_in_data_frame:
        raise ValueError(
            f"No '{data_name}' given for tags {diff_not_in_data_frame}")
    diff_not_in_ref = list(data_frame.columns.difference(
        ref_data_frame.columns))
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


def _get_attributes(cubes):
    """Extract attributes for features and labels."""
    attributes = {}

    # Extract attributes
    for cube in cubes:
        cube_attrs = cube.attributes
        tag = cube_attrs['tag']
        attributes.setdefault(tag, {})
        if cube_attrs['var_type'] in ('feature', 'label'):
            attributes[tag] = _combine_dicts(
                attributes[tag], _metadata_to_dict(cube.metadata))
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


def _get_cube_list(input_files, external_file=None,
                   merge_identical_pred_input=True):
    """Get :class:`iris.cube.CubeList` of input files."""
    cubes = iris.cube.CubeList()

    # Input files
    for filename in input_files:
        logger.info("Loading '%s'", filename)
        cube = _load_cube_with_dataset_coord(filename)
        cube.attributes['filename'] = filename
        (feature_cube, prediction_cube) = _split_cube(
            cube, merge_identical_pred_input)
        if feature_cube is not None:
            cubes.append(feature_cube)
        if prediction_cube is not None:
            cubes.append(prediction_cube)

    # External file
    cubes.extend(_get_external_cube_list(external_file))

    # Check cubes
    for cube in cubes:
        for attr in ('var_type', 'tag'):
            if attr not in cube.attributes:
                raise AttributeError(
                    f"Necessary attribute '{attr}' not given in cube "
                    f"'{cube.attributes['filename']}'")

    return cubes


def _get_external_cube_list(external_file):
    """Get external cube list."""
    if external_file is None:
        return iris.cube.CubeList()
    cubes = iris.cube.CubeList()
    with open(external_file, 'r') as infile:
        metadata_list = yaml.safe_load(infile)
    for metadata in metadata_list:
        for attr in ('data', 'dataset'):
            if attr not in metadata:
                raise AttributeError(
                    f"Entry {metadata} of external file '{external_file}' "
                    f"does not contain necessary attribute '{attr}'")
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
        cube.attributes['filename'] = external_file
        cubes.append(cube)
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
            index = pd.MultiIndex.from_product([[group],
                                                cube.coord('dataset').points],
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
                coord = cubes[0].coord(coord_name)
                returned_cube = cube
                returned_coord = coord
                break
            except iris.exceptions.CoordinateNotFoundError:
                pass
        if returned_cube is not None:
            break
    else:
        raise ValueError(
            f"No cube of {cubes} does contain coordinate 'dataset' (i.e. one "
            f"of {accepted_coord_names})")
    return (returned_cube, returned_coord)


def _load_cube_with_dataset_coord(filename):
    """Load cube with single ``dataset``-like coordinate.

    Files created by NCL cannot be read using a simple :func:`iris.load_cube`.

    """
    cubes = iris.load(filename)
    accepted_coord_names = ('dataset', 'model')

    # Handle single cube
    if len(cubes) == 1:
        (cube, coord) = _get_first_cube_with_coord(cubes, accepted_coord_names)
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
    if data_cube.ndim > 1:
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


def _create_scatterplot(x_data, y_data, numbers_as_markers=True, axes=None,
                        **kwargs):
    """Create single scatterplot."""
    if axes is None:
        (_, axes) = plt.subplots()

    # Scatterplots
    scatter_kwargs = dict(kwargs)
    scatter_kwargs.pop('label', None)
    for (idx, _) in enumerate(x_data):
        if numbers_as_markers:
            axes.text(x_data[idx], y_data[idx],
                      x_data.index.get_level_values(-1)[idx], size=7,
                      **scatter_kwargs)
        else:
            axes.scatter(x_data[idx], y_data[idx], **scatter_kwargs)

    # Regression line
    axes = _create_regplot(x_data, y_data, axes=axes, **kwargs)
    return axes


def _create_pred_input_plot(mean, error, axes, **kwargs):
    """Create plot for prediction input data (vertical lines)."""
    vspan_kwargs = dict(kwargs)
    vspan_kwargs['alpha'] = max([kwargs.get('alpha', 1.0) - 0.2, 0.1])
    vspan_kwargs.pop('label', None)
    axes.axvline(mean, **kwargs)
    axes.axvspan(mean - error, mean + error, **vspan_kwargs)
    return axes


def _create_regplot(x_data, y_data, axes=None, **kwargs):
    """Create single regression line plot."""
    if axes is None:
        (_, axes) = plt.subplots()

    # Label for plot
    kwargs = dict(kwargs)
    label = kwargs.pop('label')

    # Create regression line
    reg = linregress(x_data, y_data)
    x_range = np.max(x_data) - np.min(x_data)
    x_reg = np.linspace(np.min(x_data) - x_range, np.max(x_data) + x_range, 2)
    y_reg = reg.slope * x_reg + reg.intercept
    if label is not None:
        label += rf" ($R^2={reg.rvalue**2:.2f})$"
    else:
        label = rf"$R^2={reg.rvalue**2:.2f}$"

    # Plot
    axes.plot(x_reg, y_reg, linestyle='-', label=label, **kwargs)
    return axes


def _get_pandas_cube(pandas_object):
    """Convert :mod:`pandas` object to cube and fix coordinates."""
    cube = iris.pandas.as_cube(pandas_object)
    for coord_name in ('index', 'columns'):
        try:
            names = getattr(pandas_object, coord_name).names
        except AttributeError:
            continue
        names = [n for n in names if n is not None]
        if not names:
            continue
        new_coord_name = '-'.join(names)
        coord = cube.coord(coord_name)
        coord.var_name = new_coord_name
        coord.long_name = new_coord_name
        if not np.issubdtype(coord.dtype, np.number):
            coord.points = coord.points.astype(str)
            if coord.bounds is not None:
                coord.bounds = coord.bounds.astype(str)
    return cube


def get_input_files(cfg, patterns=None, ignore_patterns=None):
    """Get input files.

    Parameters
    ----------
    cfg : dict
        Recipe configuration.
    patterns : list of str, optional
        Use only files that match these patterns as input files.
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
    input_files = get_input_files(cfg, patterns=cfg.get('patterns'),
                                  ignore_patterns=cfg.get('ignore_patterns'))
    logger.debug("Found files:\n%s", pformat(input_files))

    # Get cubes
    external_file = _get_external_file(cfg.get('read_external_file'),
                                       cfg['auxiliary_data_dir'])
    cubes = _get_cube_list(
        input_files,
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
    for data_frame in (features, label, pred_input,
                       pred_input_err):
        data_frame.sort_index(axis=0, inplace=True)
        data_frame.sort_index(axis=1, inplace=True)

    # Check data
    (features, label, pred_input, pred_input_err) = _check_data_frames(
        features, label, pred_input, pred_input_err)
    training_data = pd.concat([features, label],
                              axis=1,
                              keys=['x', 'y'])
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


def pandas_object_to_cube(pandas_object, index_droplevel=None,
                          columns_droplevel=None, **kwargs):
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
            f"Attributes for feature '{feature}' does not inlcude necessary "
            f"key 'plot_title'")
    if 'plot_xlabel' not in feature_attrs:
        raise KeyError(
            f"Attributes for feature '{feature}' does not inlcude necessary "
            f"key 'plot_xlabel'")
    if 'plot_ylabel' not in label_attrs:
        raise KeyError(
            f"Attributes for label '{label}' does not inlcude necessary "
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


def get_groups(training_data, add_combined_group=True):
    """Extract groups from training data.

    Parameters
    ----------
    training_data : pandas.DataFrame
        Training data (features, label).
    add_combined_group : bool, optional (default: True)
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
    groups = get_groups(training_data)

    # Iterate over features
    for feature in training_data.x.columns:
        (x_data, y_data) = get_xy_data_without_nans(training_data, feature,
                                                    label)

        # Individual plots
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
                x_sub_data, y_sub_data,
                numbers_as_markers=cfg.get('numbers_as_markers', False),
                color=COLORS_ALL[idx], label=group)
            axes = _create_pred_input_plot(
                pred_input_data['mean'][feature].values,
                pred_input_data['error'][feature].values, axes, alpha=0.4,
                color=COLORS[0], label='Observation')
            set_plot_appearance(axes, attributes, plot_title=feature,
                                plot_xlabel=feature, plot_ylabel=label,
                                plot_xlim=feature, plot_ylim=label)
            legend = plt.legend(**LEGEND_KWARGS)
            filename = (f"scatterplot_{basename}_{feature}_"
                        f"{group.replace(', ', '-')}")
            plot_path = get_plot_filename(filename, cfg)
            plt.savefig(plot_path, additional_artists=[legend],
                        **cfg.get('savefig_kwargs', {}))
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Provenance
            cubes = iris.cube.CubeList([
                pandas_object_to_cube(
                    x_sub_data, index_droplevel=index_droplevel,
                    var_name=feature,
                    long_name=attributes[feature]['plot_xlabel'],
                    units=attributes[feature]['units']),
                pandas_object_to_cube(
                    y_sub_data, index_droplevel=index_droplevel,
                    var_name=label, long_name=attributes[label]['plot_ylabel'],
                    units=attributes[label]['units']),
            ])
            netcdf_path = get_diagnostic_filename(filename, cfg)
            io.iris_save(cubes, netcdf_path)
            provenance_record = get_provenance_record(
                attributes, [feature, label],
                caption=get_caption(attributes, feature, label, group=group),
                plot_type='scatter', plot_file=plot_path)
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
    groups = get_groups(training_data)
    numbers_as_markers = cfg.get('numbers_as_markers', False)

    # Iterate over features
    for feature in training_data.x.columns:
        (x_data, y_data) = get_xy_data_without_nans(training_data, feature,
                                                    label)
        (_, axes) = plt.subplots()
        if len(groups) > 1:
            axes = _create_regplot(x_data, y_data, axes=axes,
                                   color=COLORS_ALL[0],
                                   label=groups[0])
            for (idx, group) in enumerate(groups[1:]):
                axes = _create_scatterplot(
                    x_data.loc[group],
                    y_data.loc[group],
                    numbers_as_markers=numbers_as_markers,
                    axes=axes,
                    color=COLORS_ALL[idx + 1],
                    label=group,
                )
        else:
            axes = _create_scatterplot(
                x_data,
                y_data,
                numbers_as_markers=numbers_as_markers,
                axes=axes,
                color=COLORS_ALL[0],
                label=groups[0],
            )
        axes = _create_pred_input_plot(
            pred_input_data['mean'][feature].values,
            pred_input_data['error'][feature].values, axes, alpha=0.4,
            color=COLORS[0], label='Observation')
        set_plot_appearance(axes, attributes, plot_title=feature,
                            plot_xlabel=feature, plot_ylabel=label,
                            plot_xlim=feature, plot_ylim=label)
        legend = plt.legend(**LEGEND_KWARGS)
        filename = f'scatterplot_merged_{basename}_{feature}'
        plot_path = get_plot_filename(filename, cfg)
        plt.savefig(plot_path, additional_artists=[legend],
                    **cfg.get('savefig_kwargs', {}))
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Provenance
        cubes = iris.cube.CubeList([
            pandas_object_to_cube(
                x_data, index_droplevel=[0, 2], var_name=feature,
                long_name=attributes[feature]['plot_xlabel'],
                units=attributes[feature]['units']),
            pandas_object_to_cube(
                y_data, index_droplevel=[0, 2], var_name=label,
                long_name=attributes[label]['plot_ylabel'],
                units=attributes[label]['units']),
        ])
        netcdf_path = get_diagnostic_filename(filename, cfg)
        io.iris_save(cubes, netcdf_path)
        provenance_record = get_provenance_record(
            attributes, [feature, label],
            caption=get_caption(attributes, feature, label),
            plot_type='scatter', plot_file=plot_path)
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(netcdf_path, provenance_record)


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
    provenance_record = get_provenance_record(attributes, tags,
                                              caption=basename)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(csv_path, provenance_record)
    return csv_path


def standard_prediction_error(x_data, y_data):
    """Return function to calculate standard prediction error.

    The standard prediction error of a (multivariate) linear regression is the
    error when predicting a new value which is not in the original data.

    Parameters
    ----------
    x_data : numpy.ndarray
        Independent observations of predictors.
    y_data : numpy.ndarray
        Independent observations of the target variable.

    Returns
    -------
    callable
        Standard prediction error function for new observation of the
        predictors.

    """
    (x_data, y_data) = _check_training_arrays(x_data, y_data)
    lin = LinearRegression()
    lin.fit(x_data, y_data)

    # Standard error of estimates
    dof = x_data.shape[0] - x_data.shape[1]
    y_pred = lin.predict(x_data)
    see = np.sqrt(np.sum(np.square(y_data - y_pred)) / dof)

    # Get design matrix
    ones = np.ones((x_data.shape[0], 1), dtype=x_data.dtype)
    x_design = np.hstack([ones, x_data])

    # Standard prediction error for new input
    def spe(x_new):
        """Return standard prediction error."""
        x_new = _check_prediction_arrays(x_new, x_train=x_data)
        one = np.ones((1, 1), dtype=x_new.dtype)
        x_new = np.hstack([one, x_new])
        return see * (1.0 +
                      x_new @ np.linalg.inv(x_design.T @ x_design) @ x_new.T)

    return np.vectorize(spe, signature='(n)->()')


def regression_surface(x_data, y_data, n_points=50):
    """Return points of the regression surface (mean and error).

    Parameters
    ----------
    x_data : numpy.ndarray
        Independent observations of predictors.
    y_data : numpy.ndarray
        Independent observations of the target variable.
    n_points : int, optional (default: 50)
        Number of sampled points per predictor.

    Returns
    -------
    dict
        :class:`numpy.ndarray`s for the keys ``'x'``, ``'y'``,
        ``'y_minus_err'``, ``'y_plus_err'``, ``'rvalue'``, ``'slope'`` and
        ``'intercept'``.

    """
    (x_data, y_data) = _check_training_arrays(x_data, y_data)
    out = {}
    lin = LinearRegression()
    lin.fit(x_data, y_data)
    spe = standard_prediction_error(x_data, y_data)
    x_max = np.max(x_data, axis=0)
    x_min = np.min(x_data, axis=0)
    x_range = x_max - x_min
    slices = [
        slice(x_min[i] - d, x_max[i] + d, complex(0, n_points))
        for (i, d) in enumerate(x_range)
    ]
    x_lin = np.array(np.mgrid[slices])
    x_lin = x_lin.reshape(-1, np.prod(x_lin.shape[1:], dtype=int)).T
    out['x'] = x_lin
    out['y'] = lin.predict(x_lin)
    out['y_minus_err'] = out['y'] - spe(x_lin)
    out['y_plus_err'] = out['y'] + spe(x_lin)
    out['x'] = x_lin
    out['coef'] = lin.coef_
    out['intercept'] = lin.intercept_
    out['R2'] = lin.score(x_data, y_data)

    return out


def gaussian_pdf(x_data, y_data, obs_mean, obs_cov, n_points=1000):
    """Calculate Gaussian probability densitiy function for target variable.

    Parameters
    ----------
    x_data : numpy.ndarray
        Independent observations of predictors.
    y_data : numpy.ndarray
        Independent observations of the target variable.
    obs_mean : numpy.ndarray
        Mean of observational data.
    obs_cov : numpy.ndarray
        Covariance matrix of observational data.
    n_points : int, optional (default: 1000)
        Number of sampled points for target variable for PDF.

    Returns
    -------
    tuple of numpy.ndarray
        x and y values for the PDF.

    """
    (x_data, y_data) = _check_training_arrays(x_data, y_data)
    (obs_mean, obs_cov) = _check_prediction_arrays(obs_mean,
                                                   x_train=x_data,
                                                   x_cov=obs_cov)
    lin = LinearRegression()
    lin.fit(x_data, y_data)
    spe = standard_prediction_error(x_data, y_data)

    def obs_pdf(x_new):
        """Return PDF of observations P(x)."""
        gaussian = multivariate_normal(mean=obs_mean.squeeze(),
                                       cov=obs_cov.squeeze())
        return gaussian.pdf(x_new)

    def cond_pdf(x_new, y_new):
        """Return conditional PDF P(y|x)."""
        y_pred = lin.predict(x_new.reshape(1, -1))
        gaussian = multivariate_normal(mean=y_pred, cov=spe(x_new)**2)
        return gaussian.pdf(y_new)

    def comb_pdf(*args):
        """Return combined PDF P(y,x)."""
        x_new = np.array(args[:-1])
        y_new = args[-1]
        return obs_pdf(x_new) * cond_pdf(x_new, y_new)

    # Calculate PDF of target variable P(y)
    x_ranges = _get_x_ranges(obs_mean, obs_cov)
    y_range = max(y_data) - min(y_data)
    y_lin = np.linspace(min(y_data) - y_range, max(y_data) + y_range, n_points)
    y_pdf = [integrate.nquad(comb_pdf, x_ranges, args=(y, ))[0] for y in y_lin]
    return (y_lin, np.array(y_pdf))


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


def get_constraint(training_data, pred_input_data, confidence_level):
    """Get constraint on target variable.

    Parameters
    ----------
    training_data : pandas.DataFrame
        Training data (features, label).
    pred_input_data : pandas.DataFrame
        Prediction input data (mean and error).
    confidence_level : float
        Confindence level to estimate the range of the target variable.

    Returns
    -------
    tuple of float
        Lower confidence limit, best estimate and upper confidence limit of
        target variable.

    Raises
    ------
    ValueError
        Input data has not the correct shape.

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
    (x_data, y_data) = _check_training_arrays(x_data, y_data)
    pred_input_mean = pred_input_data['mean'][feature].values[0]
    pred_input_error = pred_input_data['error'][feature].values[0]
    (pred_input_mean, pred_input_error) = _check_prediction_arrays(
        pred_input_mean, x_data, pred_input_error)

    # Calculate constraint
    logger.info(
        "Calculating constraint on '%s' using %.2f%% confindence level",
        label, 100.0 * confidence_level)
    (y_lin, y_pdf) = gaussian_pdf(x_data, y_data, pred_input_mean,
                                  pred_input_error**2)
    y_cdf = cdf(y_lin, y_pdf)
    y_mean = y_lin[np.argmax(y_pdf)]
    y_index_range = np.nonzero((y_cdf >= (1.0 - confidence_level) / 2.0) &
                               (y_cdf <= (1.0 + confidence_level) / 2.0))
    y_range = y_lin[y_index_range]
    return (y_range.min(), y_mean, y_range.max())
