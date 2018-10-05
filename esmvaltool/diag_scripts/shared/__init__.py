"""Code that is shared between multiple diagnostic scripts."""
from . import names, plot
from ._base import (ProvenanceLogger, get_cfg, get_diagnostic_filename,
                    get_plot_filename, group_metadata, run_diagnostic,
                    select_metadata, sorted_group_metadata, sorted_metadata)
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
