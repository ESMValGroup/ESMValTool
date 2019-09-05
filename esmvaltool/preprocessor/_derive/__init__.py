"""Automatically derive variables."""

import importlib
import logging
from copy import deepcopy
from pathlib import Path

import iris

logger = logging.getLogger(__name__)


def _get_all_derived_variables():
    """Get all possible derived variables.

    Returns
    -------
    dict
        All derived variables with `short_name` (keys) and the associated
        python classes (values).

    """
    derivers = {}
    for path in Path(__file__).parent.glob('[a-z]*.py'):
        short_name = path.stem
        module = importlib.import_module(
            f'esmvaltool.preprocessor._derive.{short_name}')
        derivers[short_name] = getattr(module, 'DerivedVariable')
    return derivers


ALL_DERIVED_VARIABLES = _get_all_derived_variables()

__all__ = list(ALL_DERIVED_VARIABLES)


def get_required(short_name):
    """Return all required variables for derivation.

    Get all information (at least `short_name`) required for derivation and
    optionally a list of needed fx files.

    Parameters
    ----------
    short_name : str
        `short_name` of the variable to derive.

    Returns
    -------
    list
        List of dictionaries (including at least the key `short_name`)
        and occasionally mip or fx_files.

    """
    DerivedVariable = ALL_DERIVED_VARIABLES[short_name]
    variables = deepcopy(DerivedVariable().required)
    return variables


def derive(cubes,
           short_name,
           long_name,
           units,
           standard_name=None,
           fx_files=None):
    """Derive variable.

    Parameters
    ----------
    cubes: iris.cube.CubeList
        Includes all the needed variables for derivation defined in
        :func:`get_required`.
    short_name: str
        short_name
    long_name: str
        long_name
    units: str
        units
    standard_name: str, optional
        standard_name
    fx_files: dict, optional
        If required, dictionary containing fx files  with `short_name`
        (keys) and path (values) of the fx variable.

    Returns
    -------
    iris.cube.Cube
        The new derived variable.

    """
    if short_name == cubes[0].var_name:

        # FIXME (not necessary after reformat OBS PR)
        if variable['standard_name'] not in iris.std_names.STD_NAMES:
            iris.std_names.STD_NAMES[variable['standard_name']] = {
                'canonical_units': variable['units']
            }
        cubes[0].standard_name = variable['standard_name']

        return cubes[0]

    cubes = iris.cube.CubeList(cubes)
    # Preprare input cubes and add fx files if necessary
    if fx_files:
        for (fx_var, fx_path) in fx_files.items():
            if fx_path is not None:
                fx_cube = iris.load_cube(
                    fx_path,
                    constraint=iris.Constraint(
                        cube_func=lambda c, var=fx_var: c.var_name == var))
                cubes.append(fx_cube)
            else:
                logger.debug(
                    "Requested fx variable '%s' for derivation of "
                    "'%s' not found", fx_var, short_name)

    # Derive variable
    DerivedVariable = ALL_DERIVED_VARIABLES[short_name]
    cube = DerivedVariable().calculate(cubes)

    # Set standard attributes
    cube.var_name = short_name
    cube.standard_name = standard_name if standard_name else None
    cube.long_name = long_name
    cube.units = units
    for temp in cubes:
        if 'source_file' in temp.attributes:
            cube.attributes['source_file'] = temp.attributes['source_file']

    return cube
