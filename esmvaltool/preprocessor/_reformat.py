"""Simple interface to reformat and CMORize functions."""
from ..interface_scripts.cmor_check import CMORCheck
from ..interface_scripts.fixes.fix import Fix
from ..interface_scripts.variable_info import CMIP5Info, CMIP6Info

CMOR_TABLES = {
    'CMIP5': CMIP5Info(),
    'CMIP6': CMIP6Info(),
}


def _get_cmor_checker(table,
                      mip,
                      short_name,
                      fail_on_error=True,
                      automatic_fixes=False):
    """Get a CMOR checker/fixer."""
    if table not in CMOR_TABLES:
        raise NotImplementedError("No CMOR checker implemented for table {}"
                                  .format(table))

    cmor_table = CMOR_TABLES[table]
    var_info = cmor_table.get_variable(mip, short_name)

    def _checker(cube):
        return CMORCheck(
            cube,
            var_info,
            fail_on_error=fail_on_error,
            automatic_fixes=automatic_fixes)

    return _checker


def fix_file(filename, short_name, project, model):
    """Fix errors that prevent loading or can not be fixed in the cube."""
    for fix in Fix.get_fixes(
            project=project, model=model, variable=short_name):
        fix.fix_file(filename)
    # TODO: create a copy if file needs to be changed and return name to copy
    return filename


def fix_metadata(cube, short_name, project, model, cmor_table=None, mip=None):
    """Apply fixes to the metadata of the cube."""
    for fix in Fix.get_fixes(
            project=project, model=model, variable=short_name):
        fix.fix_metadata(cube)
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
        fix.fix_data(cube)
    if cmor_table and mip:
        checker = _get_cmor_checker(
            table=cmor_table,
            mip=mip,
            short_name=short_name,
            fail_on_error=False,
            automatic_fixes=True)
        checker(cube).check_data()
    return cube


def cmor_check_metadata(cube, cmor_table, mip, short_name):
    """Check if metadata conforms to CMOR."""
    checker = _get_cmor_checker(cmor_table, mip, short_name)
    checker(cube).check_metadata()
    return cube


def cmor_check_data(cube, cmor_table, mip, short_name):
    """Check if data conforms to CMOR."""
    checker = _get_cmor_checker(cmor_table, mip, short_name)
    checker(cube).check_data()
    return cube


def cmor_check(cube, cmor_table, mip, short_name):
    """Check if cube conforms to CMOR."""
    cmor_check_metadata(cube, cmor_table, mip, short_name)
    cmor_check_data(cube, cmor_table, mip, short_name)
    return cube
