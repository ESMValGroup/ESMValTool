"""Simple interface to reformat and CMORize functions."""
from ..interface_scripts.cmor_check import CMORCheck
from ..interface_scripts.fixes.fix import Fix
from ..interface_scripts.variable_info import CMIP5Info, CMIP6Info

CMOR_TABLES = {
    'CMIP5': CMIP5Info(),
    'CMIP6': CMIP6Info(),
}


def _get_cmor_checker(short_name,
                      project,
                      mip,
                      fail_on_error=True,
                      automatic_fixes=False):
    """Get a CMOR checker/fixer."""
    if project in CMOR_TABLES:
        variables_info = CMOR_TABLES[project]
    else:
        raise NotImplementedError

    var_info = variables_info.get_variable(mip, short_name)

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


def fix_metadata(cube, short_name, project, model, mip=None):
    """Apply fixes to the metadata of the cube."""
    for fix in Fix.get_fixes(
            project=project, model=model, variable=short_name):
        fix.fix_metadata(cube)
    if mip:
        checker = _get_cmor_checker(
            short_name=short_name,
            project=project,
            mip=mip,
            fail_on_error=False,
            automatic_fixes=True)
        checker(cube).check_metadata()


def fix_data(cube, short_name, project, model, mip=None):
    """Apply fixes to the data of the cube."""
    for fix in Fix.get_fixes(
            project=project, model=model, variable=short_name):
        fix.fix_data(cube)
    if mip:
        checker = _get_cmor_checker(
            short_name=short_name,
            project=project,
            mip=mip,
            fail_on_error=False,
            automatic_fixes=True)
        checker(cube).check_data()


def cmor_check_metadata(cube, short_name, project, mip):
    """Check if metadata conforms to CMOR."""
    checker = _get_cmor_checker(short_name, project, mip)
    checker(cube).check_metadata()


def cmor_check_data(cube, short_name, project, mip):
    """Check if data conforms to CMOR."""
    checker = _get_cmor_checker(short_name, project, mip)
    checker(cube).check_data()


def cmor_check(cube, short_name, project, mip):
    """Check if cube conforms to CMOR."""
    cmor_check_metadata(cube, short_name, project, mip)
    cmor_check_data(cube, short_name, project, mip)
