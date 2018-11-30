"""Code that is shared between multiple diagnostic scripts."""
from . import names, plot
from ._base import (ProvenanceLogger, extract_variables, get_cfg,
                    get_diagnostic_filename, get_plot_filename, group_metadata,
                    run_diagnostic, save_iris_cube, select_metadata,
                    sorted_group_metadata, sorted_metadata,
                    variables_available)
from ._diag import Datasets, Variable, Variables
from ._validation import apply_supermeans, get_control_exper_obs

__all__ = [
    # Main entry point for diagnostics
    'run_diagnostic',
    # Define output filenames
    'get_diagnostic_filename',
    'get_plot_filename',
    # Log provenance
    'ProvenanceLogger',
    # Select and sort input metadata
    'select_metadata',
    'sorted_metadata',
    'group_metadata',
    'sorted_group_metadata',
    'extract_variables',
    'variables_available',
    'save_iris_cube',
    'names',
    'Variable',
    'Variables',
    'Datasets',
    'get_cfg',
    # Plotting module
    'plot',
    # Validation module
    'get_control_exper_obs',
    'apply_supermeans',
]
