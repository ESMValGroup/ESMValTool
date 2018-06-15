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
    'EXP_STR',
    'MODEL_STR',
    'OBS_STR',
    'PROJECT_STR',
    'SHORT_NAME_STR',
    'Variable',
    'Variables',
    'Models',
    'get_cfg',
    'plot',
    'run_diagnostic',
]
