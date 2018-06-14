"""Code that is shared between multiple diagnostic scripts."""
from . import plot
from ._base import get_cfg, run_diagnostic
from .python_diag import *

__all__ = [
    'TIME',
    'YEAR',
    'MONTH',
    'DAY_Y',
    'DAY_M',
    'LAT',
    'LON',
    'HEIGHT',
    'Variable',
    'Variables',
    'Models',
    'get_cfg',
    'plot',
    'run_diagnostic',
]
