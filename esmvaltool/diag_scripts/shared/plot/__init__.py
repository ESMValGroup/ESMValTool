"""Module that provides common plot functions."""
# set matplotlib non-interactive backend
import matplotlib
matplotlib.use('Agg')  # noqa

from ._plot import *

__all__ = [
    'get_path_to_mpl_style',
    'quickplot'
]
