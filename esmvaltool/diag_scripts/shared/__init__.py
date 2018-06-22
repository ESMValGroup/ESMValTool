"""Code that is shared between multiple diagnostic scripts."""
from . import names
from . import plot
from ._base import get_cfg, run_diagnostic
from ._diag import Variable, Variables, Datasets

__all__ = [
    'names',
    'get_cfg',
    'plot',
    'run_diagnostic',
    'Variable',
    'Variables',
    'Datasets',
]
