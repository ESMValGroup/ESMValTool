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

MODEL_KEYS = {
    'mip',
}
VARIABLE_KEYS = {
    'reference_model',
    'alternative_model',
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
    """Load iris cubes from files"""
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
    """Concatenate all cubes after fixing metadata"""
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

    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    if (os.path.exists(filename)
            and all(cube.has_lazy_data() for cube in cubes)):
        logger.debug("Not saving cubes %s to %s to avoid data loss. "
                     "The cube is probably unchanged.", cubes, filename)
    else:
        logger.debug("Saving cubes %s to %s", cubes, filename)
        iris.save(cubes, **args)

    return filename


def save_cubes(cubes, debug=False, step=None):
    """Save iris cubes to the file specified in the _filename attribute."""
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
        _save_cubes(cubes=paths[filename], target=filename)

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
    """Write NCL metadata files to output_dir"""
    variables = list(metadata.values())
    # 'variables' is a list of dicts, but NCL does not support nested
    # dicts, so convert to dict of lists.
    keys = sorted({k for v in variables for k in v})
    input_file_info = {k: [v.get(k) for v in variables] for k in keys}
    info = {
        'input_file_info': input_file_info,
        'model_info': {},
        'variable_info': {}
    }

    # Split input_file_info into model and variable properties
    # model keys and keys with non-identical values will be stored
    # in model_info, the rest in variable_info
    for key, values in input_file_info.items():
        model_specific = any(values[0] != v for v in values)
        if (model_specific or key in MODEL_KEYS) and key not in VARIABLE_KEYS:
            info['model_info'][key] = values
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
    """Exception class for concatenation errors"""

    pass
