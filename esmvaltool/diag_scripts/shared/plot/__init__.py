"""Module that provides common plot functions."""
# set matplotlib non-interactive backend
import matplotlib
matplotlib.use('Agg')  # noqa

from ._plot import example_map_plot

__all__ = [
    'example_map_plot'
]
