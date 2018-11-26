"""Code that is shared between multiple diagnostic scripts."""
from . import names, plot
from ._base import (get_cfg, group_metadata, run_diagnostic, select_metadata,
                    sorted_group_metadata, sorted_metadata, extract_variables,
                    variables_available)
from ._write_netcdf import save_iris_cube, save_scalar_data
from ._validation import (get_control_exper_obs, apply_supermeans)
from ._diag import Datasets, Variable, Variables

__all__ = [
    'names',
    'get_cfg',
    'plot',
    'run_diagnostic',
    'Variable',
    'Variables',
    'Datasets',
    'select_metadata',
    'sorted_metadata',
    'group_metadata',
    'sorted_group_metadata',
    'extract_variables',
    'variables_available',
    'save_iris_cube',
    'save_scalar_data',
    'get_control_exper_obs',
    'apply_supermeans',
]
