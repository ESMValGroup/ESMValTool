"""Code that is shared between multiple diagnostic scripts."""
from . import io, iris_helpers, names, plot
from ._base import (
    ProvenanceLogger,
    extract_variables,
    get_cfg,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_group_metadata,
    sorted_metadata,
    variables_available,
)
from ._diag import Datasets, Variable, Variables
from ._validation import apply_supermeans, get_control_exper_obs

__all__ = [
    # Main entry point for diagnostics
    'run_diagnostic',
    # Define and write output files
    'save_figure',
    'save_data',
    'get_plot_filename',
    'get_diagnostic_filename',
    # Log provenance
    'ProvenanceLogger',
    # Select and sort input metadata
    'select_metadata',
    'sorted_metadata',
    'group_metadata',
    'sorted_group_metadata',
    'extract_variables',
    'variables_available',
    'names',
    'Variable',
    'Variables',
    'Datasets',
    'get_cfg',
    # IO module
    'io',
    # Iris helpers module
    'iris_helpers',
    # Plotting module
    'plot',
    # Validation module
    'get_control_exper_obs',
    'apply_supermeans',
]
