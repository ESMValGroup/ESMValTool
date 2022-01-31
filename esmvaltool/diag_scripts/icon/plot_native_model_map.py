#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot native ICON model output (maps).

Description
-----------
This diagnostic creates map plots of native ICON model output.

Author
------
Manuel Schlund (DLR, Germany)

Notes
-----
Input data must be 1D (the single dimensional coordinate must be the spatial
coordinate of the unstructured grid) and must contain the coordinates
``latitude`` and ``longitude``.

Configuration options in recipe
-------------------------------
mapplot_kwargs:
    Keyword arguments for :func:`psyplot.project.plot.mapplot`.

"""
import logging
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import psyplot.project as psy

from esmvaltool.diag_scripts.shared import get_plot_filename, run_diagnostic

logger = logging.getLogger(Path(__file__).stem)


def get_default_cfg(cfg):
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)
    cfg.setdefault('mapplot_kwargs', {})
    return cfg


def main(cfg):
    """Run diagnostic."""
    cfg = get_default_cfg(cfg)
    input_data = list(cfg['input_data'].values())

    # Create individual plots for each dataset
    for dataset in input_data:
        filename = dataset['filename']

        # Create plot
        psy.plot.mapplot(filename, **cfg['mapplot_kwargs'])

        # Save plot
        basename = Path(filename).stem
        plot_path = get_plot_filename(basename, cfg)
        plt.savefig(plot_path, bbox_inches='tight', orientation='landscape',
                    dpi=200)
        logger.info("Wrote %s", plot_path)
        plt.close()


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
