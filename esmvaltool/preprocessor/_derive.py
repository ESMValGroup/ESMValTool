"""Automatically derive variables."""


import logging

import iris
import yaml

from .derived_variables._derived_variable import DerivedVariable

logger = logging.getLogger(__name__)


def get_required(short_name, field=None):
    """Get variable short_name and field pairs required to derive variable."""
    frequency = field[2] if field else 'M'
    derived_var = DerivedVariable.get_derived_variable(short_name)
    return derived_var.get_required(frequency)


def derive(cubes, variable):
    """Derive variable."""
    short_name = variable['short_name']

    # Do nothing if variable is already available
    if short_name == cubes[0].var_name:
        return cubes[0]

    # Preprare input cubes and derive correct variable
    cubes = iris.cube.CubeList(cubes)
    derived_var = DerivedVariable.get_derived_variable(short_name)
    cube = derived_var.calculate(cubes)

    # Set standard attributes
    cube.var_name = short_name
    if variable['standard_name'] not in iris.std_names.STD_NAMES:
        iris.std_names.STD_NAMES[variable['standard_name']] = {
            'canonical_units': variable['units']
        }
    for attribute in ('standard_name', 'long_name', 'units'):
        setattr(cube, attribute, variable[attribute])

    # Set attributes required by preprocessor
    cube.attributes['_filename'] = variable['filename']
    cube.attributes['metadata'] = yaml.safe_dump(variable)

    return cube
