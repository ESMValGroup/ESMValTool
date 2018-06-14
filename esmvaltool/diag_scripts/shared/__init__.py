"""Code that is shared between multiple diagnostic scripts."""
from . import plot
from ._base import get_cfg, run_diagnostic
from .python_diag import *

__all__ = [
    'Variable',
    'Variables',
    'Models',
    'get_cfg',
    'plot',
    'run_diagnostic',
]
