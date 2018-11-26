#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to temperature variablity metric psi (Cox et al., 2018).

Description
-----------
Calculate global temperature variablity metric psi following Cox et al. (2018).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recips
-------------------------------

"""

import logging
import os

import iris
import numpy as np
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    group_metadata, plot, run_diagnostic, save_iris_cube, select_metadata,
    variables_available, extract_variables)

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Run the diagnostic."""
    input_data = cfg['input_data'].values()
    logger.info(input_data)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
