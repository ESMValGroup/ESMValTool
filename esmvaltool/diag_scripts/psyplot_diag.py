#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create arbitrary Psyplot plots.

Description
-----------
This diagnostic provides a high-level interface to Psyplot.

Author
------
Manuel Schlund (DLR, Germany)

Notes
-----
For each input dataset, an individual plot is created. This diagnostic supports
arbitrary variables of arbitrary datasets.

Configuration options in recipe
-------------------------------
psyplot_func: str
    Function used to plot the data. Must be a function of
    :mod:`psyplot.project.plot`. Run ``python -c "from psyplot.project import
    plot; print(plot.show_plot_methods())"`` to get a list of all currently
    supported plotting functions (make sure to run this command in your
    ESMValTool environment).
psyplot_kwargs: dict, optional
    Optional keyword arguments for the plotting function given by
    ``psyplot_func``. String arguments can include facets in curly brackets
    which will be derived from the corresponding dataset, e.g., ``clabel:
    '{long_name} [{units}]'``, ``title: '{long_name} Climatology of {dataset}
    ({start_year}-{end_year})'``.
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`. By
    default, uses ``bbox_inches: tight, dpi: 300, orientation: landscape``.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set` (affects all plots).

"""
import logging
from contextlib import redirect_stdout
from copy import deepcopy
from io import StringIO
from pathlib import Path
from pprint import pformat

import matplotlib.pyplot as plt
import psyplot.project as psy
import seaborn as sns

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)


def _get_default_cfg(cfg):
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)
    cfg.setdefault('psyplot_kwargs', {})
    cfg.setdefault('savefig_kwargs', {
        'bbox_inches': 'tight',
        'dpi': 300,
        'orientation': 'landscape',
    })
    cfg.setdefault('seaborn_settings', {})
    return cfg


def _get_plot_func(cfg):
    """Get psyplot plot function."""
    if 'psyplot_func' not in cfg:
        raise ValueError("Necessary option 'psyplot_func' missing")
    if not hasattr(psy.plot, cfg['psyplot_func']):
        with redirect_stdout(StringIO()) as str_in:
            psy.plot.show_plot_methods()
        all_plot_funcs = str_in.getvalue()
        raise AttributeError(
            f"Invalid psyplot_func '{cfg['psyplot_func']}' (must be a "
            f"function of the module psyplot.project.plot). Currently "
            f"supported:\n{all_plot_funcs}")
    logger.info(
        "Using plotting function psyplot.project.plot.%s", cfg['psyplot_func'])
    return getattr(psy.plot, cfg['psyplot_func'])


def _get_psyplot_kwargs(cfg, dataset):
    """Get keyword arguments for psyplot plotting function."""
    psyplot_kwargs = deepcopy(cfg['psyplot_kwargs'])
    for (key, val) in psyplot_kwargs.items():
        if isinstance(val, str):
            try:
                val = val.format(**dataset)
            except KeyError as exc:
                raise ValueError(
                    f"Not all necessary facets psyplot_kwargs '{key}: {val}' "
                    f"available for dataset" f"\n{pformat(dataset)}") from exc
            psyplot_kwargs[key] = val
    return psyplot_kwargs


def main(cfg):
    """Run diagnostic."""
    cfg = _get_default_cfg(cfg)
    sns.set(**cfg['seaborn_settings'])
    plot_func = _get_plot_func(cfg)

    # Create individual plots for each dataset
    input_data = list(cfg['input_data'].values())
    for dataset in input_data:
        filename = dataset['filename']
        logger.info("Creating plot '%s' for %s", cfg['psyplot_func'], filename)

        # Create plot
        psyplot_kwargs = _get_psyplot_kwargs(cfg, dataset)
        plot_func(filename, **psyplot_kwargs)

        # Save plot
        basename = Path(filename).stem
        plot_path = get_plot_filename(basename, cfg)
        plt.savefig(plot_path, **cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Provenance tracking
        caption = (f"Plot {cfg['psyplot_func']} of {dataset['long_name']} of "
                   f"dataset {dataset['dataset']} ({dataset['start_year']}-"
                   f"{dataset['end_year']}).")
        provenance_record = {
            'ancestors': [filename],
            'authors': ['schlund_manuel'],
            'caption': caption,
        }
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
