"""
Apply automatic fixes for known errors in cmorized data

All functions in this module will work even if no fixes are available
for the given dataset. Therefore is recommended to apply them to all
variables to be sure that all known errors are
fixed.

"""
from ._fixes.fix import Fix
from .check import _get_cmor_checker


def fix_file(filename, short_name, project, dataset, output_dir):
    """
    Apply fixes to the netCDF files

    This allow fixing errors that prevent loading or can not be fixed
    in the cube.

    Parameters:
    -----------
    cube: iris.cube.Cube
        Data cube to fix
    short_name: basestring
        Short name of the variable to fix
    project: basestring
    dataset: basestring
    cmor_table: basestring or None
    mip: basestring or None

    Returns:
    --------
    path:
        Filename to the fixed file. If no fix has been applied
        it will be the original

    """
    for fix in Fix.get_fixes(
            project=project, dataset=dataset, variable=short_name):
        filename = fix.fix_file(filename, output_dir)
    return filename


def fix_metadata(cube, short_name, project, dataset, cmor_table=None,
                 mip=None):
    """
    Apply fixes to the metadata of the cube.

    This fixes will not cause the data to be loaded on memory

    Parameters:
    -----------
    cube: iris.cube.Cube
        Data cube to fix
    short_name: basestring
        Short name of the variable to fix
    project: basestring
    dataset: basestring
    cmor_table: basestring or None
    mip: basestring or None

    Returns:
    --------
    iris.cube.Cube:
        Fixed cube. If no fixes were applied, returns the original cube

    """
    for fix in Fix.get_fixes(
            project=project, dataset=dataset, variable=short_name):
        cube = fix.fix_metadata(cube)
    if cmor_table and mip:
        checker = _get_cmor_checker(
            table=cmor_table,
            mip=mip,
            short_name=short_name,
            fail_on_error=False,
            automatic_fixes=True)
        checker(cube).check_metadata()
    return cube


def fix_data(cube, short_name, project, dataset, cmor_table=None, mip=None):
    """
    Apply fixes to the metadata of the cube.

    Fixes at this step require the data to be loaded onto memory,
    but the returned cube can be lazy loaded if no fixes were applied.

    Parameters:
    -----------
    cube: iris.cube.Cube
        Data cube to fix
    short_name: basestring
        Short name of the variable to fix
    project: basestring
    dataset: basestring
    cmor_table: basestring or None
    mip: basestring or None

    Returns:
    --------
    iris.cube.Cube
        Fixed cube. If no fixes were applied, returns the original cube
    """
    for fix in Fix.get_fixes(
            project=project, dataset=dataset, variable=short_name):
        cube = fix.fix_data(cube)
    if cmor_table and mip:
        checker = _get_cmor_checker(
            table=cmor_table,
            mip=mip,
            short_name=short_name,
            fail_on_error=False,
            automatic_fixes=True)
        checker(cube).check_data()
    return cube
