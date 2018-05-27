"""Code that is shared between multiple diagnostic scripts."""
from . import plot
from ._base import *

__all__ = [
    'Variable',
    'Variables',
    'get_cfg',
    'plot',
    'run_diagnostic',
]
