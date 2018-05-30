"""CMOR fixer for Iris cubes."""
from ._fixes.fix import Fix
from .check import _get_cmor_checker


def fix_file(filename, short_name, project, model, output_dir):
    """Fix errors that prevent loading or can not be fixed in the cube."""
    for fix in Fix.get_fixes(
            project=project, model=model, variable=short_name):
        filename = fix.fix_file(filename, output_dir)
    return filename


def fix_metadata(cube, short_name, project, model, cmor_table=None, mip=None):
    """Apply fixes to the metadata of the cube."""
    for fix in Fix.get_fixes(
            project=project, model=model, variable=short_name):
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


def fix_data(cube, short_name, project, model, cmor_table=None, mip=None):
    """Apply fixes to the data of the cube."""
    for fix in Fix.get_fixes(
            project=project, model=model, variable=short_name):
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
