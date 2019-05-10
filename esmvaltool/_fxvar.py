"""Module to handle fx variables."""
import logging
from copy import deepcopy

from ._data_finder import get_output_file

logger = logging.getLogger(__name__)


def _add_fxvar_keys(fx_var_dict, variable):
    """Add a couple keys specific to fx variable."""
    fx_variable = deepcopy(variable)

    # add internal recognition flag
    fx_variable['fxvar'] = True
    fx_variable['variable_group'] = fx_var_dict['short_name']
    fx_variable['short_name'] = fx_var_dict['short_name']

    # specificities of project
    if fx_variable['project'] == 'CMIP5':
        fx_variable['mip'] = 'fx'
    elif fx_variable['project'] == 'CMIP6':
        fx_variable['grid'] = variable['grid']
        if 'mip' in fx_var_dict:
            fx_variable['mip'] = fx_var_dict['mip']

    return fx_variable


def _update_fx_files(fx_varlist, config_user, parent_variable):
    """Get the fx files dict for a list of fx variables."""
    fx_files_dict = {}
    for fx_variable in fx_varlist:
        fx_files_dict[fx_variable['short_name']] = get_output_file(
            fx_variable,
            config_user['preproc_dir'],
            parent_variable)
    return fx_files_dict
