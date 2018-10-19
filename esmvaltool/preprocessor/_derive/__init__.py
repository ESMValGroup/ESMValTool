"""Automatically derive variables."""


import logging

import iris
import yaml

from ._derived_variable_base import DerivedVariableBase

logger = logging.getLogger(__name__)


def get_required(short_name, field=None):
    """Get variable short_name and field pairs required to derive variable.

    It is also possible to process fx variables using the tuple ('fx_files',
    [...]), e.g. ('fx_files', ['sftlf, 'orog']).

    Parameters
    ----------
    short_name : str
        `short_name` of the derived variable.
    field : str, optional
        `field_type` of the derived variable.

    Returns
    -------
    list of tuples
        List of tuples `(short_name, field)` of all variables required for
        derivation, in case of fx variables also the tuple `('fx_files',
        [...]).

    """
    frequency = field[2] if field else 'M'
    derived_var = DerivedVariableBase.get_derived_variable(short_name)
    return derived_var.get_required(frequency)


def derive(cubes, variable):
    """Derive variable.

    Parameters
    ----------
    cubes : iris.cube.CubeList
        Includes all the needed variables for derivation defined in
        :func:`get_required`.
    variable : dict
        All information of the derived variable.

    Returns
    -------
    iris.cube.Cube
        The new derived variable.

    """
    short_name = variable['short_name']

    # Do nothing if variable is already available
    if short_name == cubes[0].var_name:
        return cubes[0]

    # Preprare input cubes and derive correct variable
    cubes = iris.cube.CubeList(cubes)
    derived_var = DerivedVariableBase.get_derived_variable(short_name)
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
