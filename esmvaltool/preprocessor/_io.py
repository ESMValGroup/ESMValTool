"""Functions for loading and saving cubes"""
import logging
import os
import shutil
from itertools import groupby

import iris
import iris.exceptions
import yaml

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
    cubes = iris.load_raw(files, constraints=constraints, callback=callback)
    iris.util.unify_time_units(cubes)
    if not cubes:
        raise Exception('Can not load cubes from {0}'.format(files))

    for cube in cubes:
        cube.attributes['_filename'] = filename
        cube.attributes['metadata'] = yaml.safe_dump(metadata)
        # TODO add block below when using iris 2.0
        # always set fillvalue to 1e+20
        # if np.ma.is_masked(cube.data):
        #     np.ma.set_fill_value(cube.data, GLOBAL_FILL_VALUE)

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

    # TODO replace block when using iris 2.0
    for filename in paths:
        # _save_cubes(cubes=paths[filename], target=filename,
        #             fill_value=GLOBAL_FILL_VALUE)
        _save_cubes(
            cubes=paths[filename],
            target=filename,
            zlib=compress,
            optimize_access=optimize_access)

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
    variables = list(metadata.values())
    # 'variables' is a list of dicts, but NCL does not support nested
    # dicts, so convert to dict of lists.
    keys = sorted({k for v in variables for k in v})
    input_file_info = {k: [v.get(k) for v in variables] for k in keys}
    fx_file_list = input_file_info.pop('fx_files', None)
    if fx_file_list:
        for fx_files in fx_file_list:
            for key in fx_files:
                if key not in input_file_info:
                    input_file_info[key] = []
                input_file_info[key].append(fx_files[key])
    # NCL cannot handle nested arrays so delete for now
    # TODO: switch to NCL list type
    input_file_info.pop('institute', None)
    input_file_info.pop('modeling_realm', None)
    info = {
        'input_file_info': input_file_info,
        'dataset_info': {},
        'variable_info': {}
    }

    # Split input_file_info into dataset and variable properties
    # dataset keys and keys with non-identical values will be stored
    # in dataset_info, the rest in variable_info
    for key, values in input_file_info.items():
        dataset_specific = any(values[0] != v for v in values)
        if (dataset_specific or key in DATASET_KEYS) and \
                key not in VARIABLE_KEYS:
            info['dataset_info'][key] = values
        else:
            # Select a value that is filled
            attribute_value = None
            for value in values:
                if value is not None:
                    attribute_value = value
                    break
            info['variable_info'][key] = attribute_value

    short_name = info['variable_info']['short_name']
    filename = os.path.join(output_dir, short_name + '_info.ncl')
    write_ncl_settings(info, filename)

    return filename


class ConcatenationError(Exception):
    """Exception class for concatenation errors."""
