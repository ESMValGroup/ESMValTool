"""Module that provides common plot functions."""
# set matplotlib non-interactive backend
import matplotlib
matplotlib.use('Agg')  # noqa

from ._plot import quickplot

__all__ = ['quickplot']
