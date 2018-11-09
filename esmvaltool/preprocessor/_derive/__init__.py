"""Automatically derive variables."""

import importlib
import logging
import os

import iris
import yaml

from ._derived_variable_base import DerivedVariableBase

logger = logging.getLogger(__name__)

ALL_DERIVED_VARIABLES = {}


def get_required(short_name, field=None):
    """Return all required variables for derivation.

    Get all information (at least `short_name`) required for derivation and
    optionally a list of needed fx files.

    Parameters
    ----------
    short_name : str
        `short_name` of the derived variable.
    field : str, optional
        `field_type` of the derived variable.

    Returns
    -------
    dict
        Dictionary containing a :obj:`list` of dictionaries (including at least
        the key `short_name`) with the key `vars` and optionally a :obj:`list`
        of fx variables with the key `fx_files`.

    """
    frequency = field[2] if field else 'M'
    derived_var = DerivedVariableBase.get_derived_variable(short_name)
    return derived_var.get_required(frequency)


def derive(cubes, variable, fx_files=None):
    """Derive variable.

    Parameters
    ----------
    cubes : iris.cube.CubeList
        Includes all the needed variables for derivation defined in
        :func:`get_required`.
    variable : dict
        All information of the derived variable. Required keys are
        `short_name`, `standard_name`, `long_name`, and `units`.
    fx_files : dict, optional
        If required, dictionary containing fx files  with `short_name`
        (keys) and path (values) of the fx variable.

    Returns
    -------
    iris.cube.Cube
        The new derived variable.

    """
    short_name = variable['short_name']

    # Do nothing if variable is already available
    if short_name == cubes[0].var_name:
        return cubes[0]

    # Preprare input cubes and add fx files if necessary
    cubes = iris.cube.CubeList(cubes)
    if fx_files:
        for (fx_var, fx_path) in fx_files.items():
            if fx_path is not None:
                cubes.append(iris.load_cube(fx_path))
            else:
                logger.debug(
                    "Requested fx variable '%s' for derivation of "
                    "'%s' not found", fx_var, short_name)

    # Derive correct variable
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

    return cube


def get_all_derived_variables():
    """Get all possible derived variables.

    Returns
    -------
    dict
        All derived variables with `short_name` (keys) and the associated
        python classes (values).

    """
    if ALL_DERIVED_VARIABLES:
        return ALL_DERIVED_VARIABLES
    current_path = os.path.dirname(os.path.realpath(__file__))
    for var_file in os.listdir(current_path):
        var_name = os.path.splitext(var_file)[0]
        try:
            var_module = importlib.import_module(
                'esmvaltool.preprocessor._derive.{}'.format(var_name))
            try:
                derived_var = getattr(var_module, 'DerivedVariable')(var_name)
                ALL_DERIVED_VARIABLES[var_name] = derived_var
            except AttributeError:
                pass
        except ImportError:
            pass

    return ALL_DERIVED_VARIABLES
