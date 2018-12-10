"""Convenience functions for writing netcdf files."""
import fnmatch
import logging
import os
from datetime import datetime

import iris

logger = logging.getLogger(__name__)

VAR_KEYS = [
    'long_name',
    'standard_name',
    'units',
]
NECESSARY_KEYS = VAR_KEYS + [
    'dataset',
    'filename',
    'project',
    'short_name',
]


def _add_standard_name_to_iris(standard_name, units, cube=None):
    """Add (invalid) `standard_name` to iris."""
    if standard_name is None:
        return
    if standard_name not in iris.std_names.STD_NAMES:
        iris.std_names.STD_NAMES[standard_name] = {
            'canonical_units': units,
        }
    if cube is not None:
        cube.standard_name = standard_name


def _has_necessary_attributes(metadata, log_level='debug'):
    """Check if dataset metadata has necessary attributes."""
    for dataset in metadata:
        for key in NECESSARY_KEYS:
            if key not in dataset:
                getattr(logger, log_level)("Dataset '%s' does not have "
                                           "necessary attribute '%s'", dataset,
                                           key)
                return False
    return True


def _has_necessary_var_attributes(metadata, log_level='debug'):
    """Check if dataset metadata has necessary variable attributes."""
    for var_key in VAR_KEYS + ['short_name']:
        if var_key not in metadata:
            getattr(logger, log_level)("Metadata '%s' do not have "
                                       "necessary variable attribute '%s'",
                                       metadata, var_key)
            return False
    return True


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
    str or None
        Full path to the file or `None` if file not found.

    """
    files = get_all_ancestor_files(cfg, pattern=pattern)
    if not files:
        logger.warning(
            "No file with requested name %s found in ancestor "
            "directories", pattern)
        return None
    if len(files) != 1:
        logger.warning(
            "Multiple files with requested pattern %s found (%s), returning "
            "first appearance", pattern, files)
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

    """
    if root is None:
        all_files = get_all_ancestor_files(cfg, pattern)
    else:
        all_files = []
        for (base, _, files) in os.walk(root):
            if pattern is not None:
                files = fnmatch.filter(files, pattern)
            files = fnmatch.filter(files, '*.nc')
            files = [os.path.join(base, f) for f in files]
            all_files.extend(files)

    # Iterate over netcdf files
    metadata = []
    for path in all_files:
        cube = iris.load_cube(path)
        _add_standard_name_to_iris(
            cube.attributes.pop('invalid_standard_name'), cube.units, cube)
        dataset_info = dict(cube.attributes)
        for var_key in VAR_KEYS:
            dataset_info[var_key] = getattr(cube, var_key)
        dataset_info['short_name'] = getattr(cube, 'var_name')
        dataset_info['filename'] = path

        # Check if necessary keys are available
        if _has_necessary_attributes([dataset_info], log_level='warning'):
            metadata.append(dataset_info)
        else:
            logger.warning("Skipping '%s'", path)

    return metadata


def metadata_to_netcdf(cube, metadata, cfg):
    """Convert list of metadata to netcdf files.

    Parameters
    ----------
    cube : iris.cube.Cube
        Cube to be written.
    metadata : dict
        Metadata for the cube.
    cfg : dict
        Diagnostic script configuration.

    """
    if not _has_necessary_attributes([metadata], 'error'):
        logger.error("Cannot save cube %s", cube)
        return
    _add_standard_name_to_iris(metadata['standard_name'], metadata['units'],
                               cube)
    for var_key in VAR_KEYS:
        setattr(cube, var_key, metadata.pop(var_key))
    setattr(cube, 'var_name', metadata.pop('short_name'))
    for (attr, val) in metadata.items():
        if isinstance(val, bool):
            metadata[attr] = str(val)
    cube.attributes.update(metadata)
    save_iris_cube(cube, metadata['filename'], cfg)


def save_iris_cube(cube, path, cfg):
    """Save `iris.cube.Cube` and append ESMValTool information.

    Parameters
    ----------
    cube : iris.cube.Cube
        Cube to be saved.
    path : str
        Desired path.
    cfg : dict
        Diagnostic script configuration.

    """
    if not cfg['write_netcdf']:
        logger.warning(
            "Could not write netcdf file '%s', 'write_netcdf' is "
            "set to 'False' in user configuration file.", path)
        return
    attr = {
        'created_by':
        'ESMValTool version {}'.format(cfg['version']) +
        ', diagnostic {}'.format(cfg['script']),
        'creation_date':
        datetime.utcnow().isoformat(' ') + ' UTC',
    }
    cube.attributes.update(attr)
    iris.save(cube, path)
    logger.info("Wrote %s", path)


def save_scalar_data(data, path, cfg, var_attrs):
    """Save scalar data for multiple datasets.

    Create 1D cube with the auxiliary dimension `dataset` and save scalar data
    for every appearing dataset.

    Parameters
    ----------
    data : dict
        Scalar data (values) and corresponding datasets (keys).
    path : str
        Desired path.
    cfg : dict
        Diagnostic script configuration.
    var_attrs : dict
        Attributes for the variable (`var_name`, `standard_name`, `long_name`
        or `units`).

    """
    if not _has_necessary_var_attributes(var_attrs, 'error'):
        logger.error("Cannot write file '%s'", path)
        return
    dataset_coord = iris.coords.AuxCoord(list(data), long_name='dataset')
    _add_standard_name_to_iris(var_attrs['standard_name'], var_attrs['units'])
    if 'short_name' in var_attrs:
        var_attrs['var_name'] = var_attrs.pop('short_name')
    cube = iris.cube.Cube(
        list(data.values()),
        aux_coords_and_dims=[(dataset_coord, 0)],
        **var_attrs)
    cube.attributes['filename'] = path
    save_iris_cube(cube, path, cfg)
