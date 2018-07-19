"""Module that provides common plot functions."""
# set matplotlib non-interactive backend
import matplotlib
matplotlib.use('Agg')  # noqa

from ._plot import (
    get_path_to_mpl_style,
    get_dataset_style,
    quickplot,
    multi_dataset_scatterplot,
    scatterplot,
)

__all__ = [
    'get_path_to_mpl_style',
    'get_dataset_style',
    'quickplot',
    'multi_dataset_scatterplot',
    'scatterplot',
]
