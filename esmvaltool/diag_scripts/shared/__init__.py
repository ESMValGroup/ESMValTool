"""Code that is shared between multiple diagnostic scripts."""
from . import plot
from ._base import (get_cfg, group_metadata, run_diagnostic, select_metadata,
                    sorted_group_metadata, sorted_metadata,
                    variables_available, save_iris_cube)
from ._validation import (get_control_exper_obs, apply_supermeans)

__all__ = [
    'get_cfg',
    'plot',
    'run_diagnostic',
    'select_metadata',
    'sorted_metadata',
    'group_metadata',
    'sorted_group_metadata',
    'variables_available',
    'save_iris_cube',
    'get_control_exper_obs',
    'apply_supermeans',
]
