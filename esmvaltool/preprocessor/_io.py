"""Functions for loading and saving cubes"""
import logging
import os

import iris

iris.FUTURE.netcdf_promote = True
iris.FUTURE.netcdf_no_unlimited = True

logger = logging.getLogger(__name__)


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


def load_cubes(files, filename, constraints=None, callback=None):
    """Load iris cubes from files"""
    logger.debug("Loading and concatenating:\n%s", "\n".join(files))
    cubes = iris.load_raw(files, constraints=constraints, callback=callback)
    iris.util.unify_time_units(cubes)
    cubes = cubes.concatenate()
    if not cubes:
        raise Exception('Can not load cubes from {0}'.format(files))
    for cube in cubes:
        cube.attributes['_filename'] = filename
    return cubes


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
            filename = cube.attributes.get('_filename')
            filename = os.path.splitext(filename)[0]
            filename = os.path.join(filename, step + '.nc')
        else:
            filename = cube.attributes.pop('_filename')
        if filename not in paths:
            paths[filename] = []
        paths[filename].append(cube)

    for filename in paths:
        _save_cubes(cubes=paths[filename], target=filename)

    return list(paths)
