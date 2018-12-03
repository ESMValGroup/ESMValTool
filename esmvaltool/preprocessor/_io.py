"""Functions for loading and saving cubes."""
import copy
import logging
import os
import shutil
from itertools import groupby

import iris
import iris.exceptions
import yaml

from .._config import use_legacy_iris
from .._task import write_ncl_settings

logger = logging.getLogger(__name__)

GLOBAL_FILL_VALUE = 1e+20

DATASET_KEYS = {
    'mip',
}
VARIABLE_KEYS = {
    'reference_dataset',
    'alternative_dataset',
}


def _get_attr_from_field_coord(ncfield, coord_name, attr):
    if coord_name is not None:
        attrs = ncfield.cf_group[coord_name].cf_attrs()
        attr_val = [value for (key, value) in attrs if key == attr]
        if attr_val:
            return attr_val[0]
    return None


def concatenate_callback(raw_cube, field, _):
    """Use this callback to fix anything Iris tries to break."""
    # Remove attributes that cause issues with merging and concatenation
    for attr in ['creation_date', 'tracking_id', 'history']:
        if attr in raw_cube.attributes:
            del raw_cube.attributes[attr]
    for coord in raw_cube.coords():
        # Iris chooses to change longitude and latitude units to degrees
        # regardless of value in file, so reinstating file value
        if coord.standard_name in ['longitude', 'latitude']:
            units = _get_attr_from_field_coord(field, coord.var_name, 'units')
            if units is not None:
                coord.units = units


def load_cubes(files, filename, metadata, constraints=None, callback=None):
    """Load iris cubes from files."""
    logger.debug("Loading:\n%s", "\n".join(files))
    if constraints is not None:
        cubes = iris.load_raw(files, constraints=constraints['standard_name'],
                              callback=callback)
        # test for long name
        if not cubes:
            iris_constraint = iris.Constraint(
                cube_func=(lambda c: c.long_name == constraints['long_name']))
            cubes = iris.load(files, constraints=iris_constraint,
                              callback=callback)
        # last resort - test for short name
        if not cubes:
            iris_constraint = iris.Constraint(
                cube_func=(lambda c: c.var_name == constraints['short_name']))
            cubes = iris.load(files, constraints=iris_constraint,
                              callback=callback)
    else:
        cubes = iris.load_raw(files, constraints=constraints,
                              callback=callback)
    iris.util.unify_time_units(cubes)
    if not cubes:
        raise Exception('Can not load cubes from {0}'.format(files))

    for cube in cubes:
        cube.attributes['_filename'] = filename
        cube.attributes['metadata'] = yaml.safe_dump(metadata)

    return cubes


def concatenate(cubes):
    """Concatenate all cubes after fixing metadata."""
    try:
        cube = iris.cube.CubeList(cubes).concatenate_cube()
        return cube
    except iris.exceptions.ConcatenateError as ex:
        logger.error('Can not concatenate cubes: %s', ex)
        logger.error('Differences: %s', ex.differences)
        logger.error('Cubes:')
        for cube in cubes:
            logger.error(cube)
        raise ConcatenationError('Can not concatenate cubes {0}'.format(cubes))


def _save_cubes(cubes, **args):
    """Save iris cube to file."""
    filename = args['target']
    optimize_accesss = args.pop('optimize_access')

    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    if (os.path.exists(filename)
            and all(cube.has_lazy_data() for cube in cubes)):
        logger.debug(
            "Not saving cubes %s to %s to avoid data loss. "
            "The cube is probably unchanged.", cubes, filename)
    else:
        logger.debug("Saving cubes %s to %s", cubes, filename)
        if optimize_accesss:
            cube = cubes[0]
            if optimize_accesss == 'map':
                dims = set(
                    cube.coord_dims('latitude') + cube.coord_dims('longitude'))
            elif optimize_accesss == 'timeseries':
                dims = set(cube.coord_dims('time'))
            else:
                dims = tuple()
                for coord_dims in (
                        cube.coord_dims(dimension)
                        for dimension in optimize_accesss.split(' ')):
                    dims += coord_dims
                dims = set(dims)

            args['chunksizes'] = tuple(
                length if index in dims else 1
                for index, length in enumerate(cube.shape))
        iris.save(cubes, **args)

    return filename


def save(cubes, optimize_access=None, compress=False, debug=False, step=None):
    """
    Save iris cubes to file.

    Path is taken from the _filename attributte in the code.

    Parameters
    ----------
    cubes: iterable of iris.cube.Cube
        Data cubes to be saved

    optimize_access: str
        Set internal NetCDF chunking to favour a reading scheme

        Values can be map or timeseries, which improve performance when
        reading the file one map or time series at a time.
        Users can also provide a coordinate or a list of coordinates. In that
        case the better performance will be avhieved by loading all the values
        in that coordinate at a time

    compress: bool, optional
        Use NetCDF internal compression.

    debug: bool, optional
        Inform the function if this save is an intermediate save

    step: int, optional
        Number of the preprocessor step.

        Only used if debug is True

    Returns
    -------
    list
        List of paths
    """
    paths = {}
    for cube in cubes:
        if '_filename' not in cube.attributes:
            raise ValueError("No filename specified in cube {}".format(cube))
        if debug:
            dirname = os.path.splitext(cube.attributes.get('_filename'))[0]
            if os.path.exists(dirname) and os.listdir(dirname):
                num = int(sorted(os.listdir(dirname)).pop()[:2]) + 1
            else:
                num = 0
            filename = os.path.join(dirname, '{:02}_{}.nc'.format(num, step))
        else:
            filename = cube.attributes.pop('_filename')
        if filename not in paths:
            paths[filename] = []
        paths[filename].append(cube)

    for filename in paths:
        if use_legacy_iris():
            _save_cubes(
                cubes=paths[filename],
                target=filename,
                zlib=compress,
                optimize_access=optimize_access)
        else:
            _save_cubes(
                cubes=paths[filename],
                target=filename,
                optimize_access=optimize_access,
                fill_value=GLOBAL_FILL_VALUE)

    return list(paths)


def cleanup(files, remove=None):
    """Clean up after running the preprocessor."""
    if remove is None:
        remove = []

    for path in remove:
        if os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.isfile(path):
            os.remove(path)

    return files


def extract_metadata(files, write_ncl=False):
    """Extract the metadata attribute from cubes and write to file."""
    output_files = []
    for output_dir, filenames in groupby(files, os.path.dirname):
        metadata = {}
        for filename in filenames:
            cube = iris.load_cube(filename)
            raw_cube_metadata = cube.attributes.get('metadata')
            if raw_cube_metadata:
                cube_metadata = yaml.safe_load(raw_cube_metadata)
                metadata[filename] = cube_metadata

        output_filename = os.path.join(output_dir, 'metadata.yml')
        output_files.append(output_filename)
        with open(output_filename, 'w') as file:
            yaml.safe_dump(metadata, file)
        if write_ncl:
            output_files.append(_write_ncl_metadata(output_dir, metadata))

    return output_files


def _write_ncl_metadata(output_dir, metadata):
    """Write NCL metadata files to output_dir."""
    variables = copy.deepcopy(list(metadata.values()))

    for variable in variables:
        fx_files = variable.pop('fx_files', {})
        for fx_type in fx_files:
            variable[fx_type] = fx_files[fx_type]

    info = {'input_file_info': variables}

    # Split input_file_info into dataset and variable properties
    # dataset keys and keys with non-identical values will be stored
    # in dataset_info, the rest in variable_info
    variable_info = {}
    info['variable_info'] = [variable_info]
    info['dataset_info'] = []
    for variable in variables:
        dataset_info = {}
        info['dataset_info'].append(dataset_info)
        for key in variable:
            dataset_specific = any(
                variable[key] != var.get(key, object()) for var in variables)
            if ((dataset_specific or key in DATASET_KEYS)
                    and key not in VARIABLE_KEYS):
                dataset_info[key] = variable[key]
            else:
                variable_info[key] = variable[key]

    filename = os.path.join(output_dir,
                            variable_info['short_name'] + '_info.ncl')
    write_ncl_settings(info, filename)

    return filename


class ConcatenationError(Exception):
    """Exception class for concatenation errors."""
