"""Convenience functions for writing netcdf files."""
import fnmatch
import logging
import os
from pprint import pformat

import iris
import numpy as np

from .iris_helpers import unify_1d_cubes

logger = logging.getLogger(__name__)

VAR_KEYS = [
    'long_name',
    'units',
]
NECESSARY_KEYS = VAR_KEYS + [
    'dataset',
    'filename',
    'project',
    'short_name',
]


def _has_necessary_attributes(metadata,
                              only_var_attrs=False,
                              log_level='debug'):
    """Check if dataset metadata has necessary attributes."""
    output = True
    keys_to_check = (VAR_KEYS +
                     ['short_name'] if only_var_attrs else NECESSARY_KEYS)
    for dataset in metadata:
        for key in keys_to_check:
            if key not in dataset:
                getattr(logger, log_level)(
                    "Dataset '%s' does not have necessary attribute '%s'",
                    dataset, key)
                output = False
    return output


def get_all_ancestor_files(cfg, pattern=None):
    """Return a list of all files in the ancestor directories.

    Parameters
    ----------
    cfg : dict
        Diagnostic script configuration.
    pattern : str, optional
        Only return files which match a certain pattern.

    Returns
    -------
    list of str
        Full paths to the ancestor files.

    """
    ancestor_files = []
    input_dirs = [
        d for d in cfg['input_files'] if not d.endswith('metadata.yml')
    ]
    for input_dir in input_dirs:
        for (root, _, files) in os.walk(input_dir):
            if pattern is not None:
                files = fnmatch.filter(files, pattern)
            files = [os.path.join(root, f) for f in files]
            ancestor_files.extend(files)
    return ancestor_files


def get_ancestor_file(cfg, pattern):
    """Return a desired file in the ancestor directories.

    Parameters
    ----------
    cfg : dict
        Diagnostic script configuration.
    pattern : str
        Pattern which specifies the name of the file.

    Returns
    -------
    str
        Full path to the file.

    Raises
    ------
    ValueError
        No or more than one file found.

    """
    files = get_all_ancestor_files(cfg, pattern=pattern)
    if len(files) != 1:
        raise ValueError(
            f"Expected to find exactly one ancestor file for pattern "
            f"'{pattern}', got {len(files):d}:\n{pformat(files)}")
    return files[0]


def netcdf_to_metadata(cfg, pattern=None, root=None):
    """Convert attributes of netcdf files to list of metadata.

    Parameters
    ----------
    cfg : dict
        Diagnostic script configuration.
    pattern : str, optional
        Only consider files which match a certain pattern.
    root : str, optional (default: ancestor directories)
        Root directory for the search.

    Returns
    -------
    list of dict
        List of dataset metadata.

    Raises
    ------
    ValueError
        Necessary attributes are missing.

    """
    if root is None:
        all_files = get_all_ancestor_files(cfg, pattern)
    else:
        all_files = []
        for (base, _, files) in os.walk(root):
            if pattern is not None:
                files = fnmatch.filter(files, pattern)
            files = [os.path.join(base, f) for f in files]
            all_files.extend(files)
    all_files = fnmatch.filter(all_files, '*.nc')

    # Iterate over netcdf files
    metadata = []
    for path in all_files:
        cube = iris.load_cube(path)
        dataset_info = dict(cube.attributes)
        for var_key in VAR_KEYS:
            dataset_info[var_key] = str(getattr(cube, var_key))
        dataset_info['short_name'] = cube.var_name
        dataset_info['standard_name'] = cube.standard_name
        dataset_info['filename'] = path
        metadata.append(dataset_info)

    # Check if necessary keys are available
    if not _has_necessary_attributes(metadata, log_level='error'):
        raise ValueError("Necessary attributes are missing for metadata")

    return metadata


def metadata_to_netcdf(cube, metadata):
    """Convert single metadata dictionary to netcdf file.

    Parameters
    ----------
    cube : iris.cube.Cube
        Cube to be written.
    metadata : dict
        Metadata for the cube.

    Raises
    ------
    ValueError
        Saving of cube not possible because of invalid metadata.

    """
    metadata = dict(metadata)
    if not _has_necessary_attributes([metadata], log_level='error'):
        raise ValueError(f"Cannot save cube {cube.summary(shorten=True)}")
    for var_key in VAR_KEYS:
        setattr(cube, var_key, metadata.pop(var_key))
    cube.var_name = metadata.pop('short_name')
    cube.standard_name = None
    if 'standard_name' in metadata:
        standard_name = metadata.pop('standard_name')
        try:
            cube.standard_name = standard_name
        except ValueError:
            logger.warning(
                "Got invalid standard_name '%s', setting it to 'None'",
                standard_name)
            cube.attributes['invalid_standard_name'] = standard_name
    for (attr, val) in metadata.items():
        if isinstance(val, bool):
            metadata[attr] = str(val)
    cube.attributes.update(metadata)
    iris_save(cube, metadata['filename'])


def iris_save(source, path):
    """Save :mod:`iris` objects with correct attributes.

    Parameters
    ----------
    source : iris.cube.Cube or iterable of iris.cube.Cube
        Cube(s) to be saved.
    path : str
        Path to the new file.

    """
    if isinstance(source, iris.cube.Cube):
        source.attributes['filename'] = path
    else:
        for cube in source:
            cube.attributes['filename'] = path
    iris.save(source, path)
    logger.info("Wrote %s", path)


def save_1d_data(cubes, path, coord_name, var_attrs, attributes=None):
    """Save 1D data for multiple datasets.

    Create 2D cube with the dimensionsal coordinate ``coord_name`` and the
    auxiliary coordinate ``dataset`` and save 1D data for every dataset given.
    The cube is filled with missing values where no data exists for a dataset
    at a certain point.

    Note
    ----
    Does not check metadata of the ``cubes``, i.e. different names or units
    will be ignored.

    Parameters
    ----------
    cubes : dict of iris.cube.Cube
        1D :class.:`iris.cube.Cube`s (values) and corresponding dataset names
        (keys).
    path : str
        Path to the new file.
    coord_name : str
        Name of the coordinate.
    var_attrs : dict
        Attributes for the variable (``short_name``, ``long_name``, and
        ``units``).
    attributes : dict, optional
        Additional attributes for the cube.

    Raises
    ------
    ValueError
        Empty list of cubes given or necessary variable attributes are missing.

    """
    var_attrs = dict(var_attrs)
    if not cubes:
        raise ValueError("Cannot save 1D data, no cubes given")
    if not _has_necessary_attributes(
            [var_attrs], only_var_attrs=True, log_level='error'):
        raise ValueError(
            f"Cannot save 1D data to {path} because necessary variable "
            f"attributes are missing")
    datasets = list(cubes.keys())
    cube_list = iris.cube.CubeList(list(cubes.values()))
    cube_list = unify_1d_cubes(cube_list, coord_name)
    data = [c.data for c in cube_list]
    dataset_coord = iris.coords.AuxCoord(datasets, long_name='dataset')
    coord = cube_list[0].coord(coord_name)
    if attributes is None:
        attributes = {}
    var_attrs['var_name'] = var_attrs.pop('short_name')

    # Create new cube
    cube = iris.cube.Cube(np.ma.array(data),
                          aux_coords_and_dims=[(dataset_coord, 0), (coord, 1)],
                          attributes=attributes,
                          **var_attrs)
    iris_save(cube, path)


def save_scalar_data(data, path, var_attrs, aux_coord=None, attributes=None):
    """Save scalar data for multiple datasets.

    Create 1D cube with the auxiliary dimension ``dataset`` and save scalar
    data for every dataset given.

    Note
    ----
    Missing values can be added by :obj:`numpy.nan`.

    Parameters
    ----------
    data : dict
        Scalar data (values) and corresponding datasets (keys).
    path : str
        Path to the new file.
    var_attrs : dict
        Attributes for the variable (``short_name``, ``long_name``, and
        ``units``).
    aux_coord : iris.coords.AuxCoord, optional
        Optional auxiliary coordinate.
    attributes : dict, optional
        Additional attributes for the cube.

    Raises
    ------
    ValueError
        No data given or necessary variable attributes are missing.

    """
    var_attrs = dict(var_attrs)
    if not data:
        raise ValueError("Cannot save scalar data, no data given")
    if not _has_necessary_attributes(
            [var_attrs], only_var_attrs=True, log_level='error'):
        raise ValueError(
            f"Cannot save scalar data to {path} because necessary variable "
            f"attributes are missing")
    dataset_coord = iris.coords.AuxCoord(list(data), long_name='dataset')
    if attributes is None:
        attributes = {}
    var_attrs['var_name'] = var_attrs.pop('short_name')
    coords = [(dataset_coord, 0)]
    if aux_coord is not None:
        coords.append((aux_coord, 0))
    cube = iris.cube.Cube(np.ma.masked_invalid(list(data.values())),
                          aux_coords_and_dims=coords,
                          attributes=attributes,
                          **var_attrs)
    iris_save(cube, path)
