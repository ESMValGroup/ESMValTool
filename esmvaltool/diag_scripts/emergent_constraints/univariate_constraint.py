#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to evaluate a single emergent constraint.

Description
-----------
Establish emergent constraint for an arbitrary input variable and an arbitrary
target variable. All input data must be marked with a ``var_type`` (either
``feature``, ``label``, ``prediction_input`` or ``prediction_input_error``) and
a ``tag``, which describes the data. This diagnostic supports only a single
``tag`` for ``label`` and ``feature``. For every ``tag``, a
``'reference_dataset'`` can be specified, which will be automatically
considered as ``prediction_input``. If ``reference_dataset`` contains ``'|'``
(e.g. ``'OBS1|OBS2'``), mutliple datasets are considered as
``prediction_input`` (in this case ``'OBS1'`` and ``'OBS2'``).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
all_data_label : str, optional (default: 'all')
    Label used in plots when all input data is considered. Only relevant if
    ``group_by`` is not used.
confidence_level : float, optional (default: 0.66)
    Confidence level for estimation of target variable.
group_by : str, optional
    Group input data by an attribute (e.g. produces separate plots for the
    individual groups, etc.).
ignore_patterns : list of str, optional
    Patterns matched against ancestor files. Those files are ignored.
merge_identical_pred_input : bool, optional (default: True)
    Use identical prediction_input values as single value.
numbers_as_markers : bool, optional (default: False)
    Use numbers as markers in scatterplots.
patterns : list of str, optional
    Patterns matched against ancestor files.
read_external_file : str, optional
    Read input datasets from external file given as absolute path or relative
    path. In the latter case, ``'auxiliary_data_dir'`` from the user
    configuration file is used as base directory.
savefig_kwargs : dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots), see
    <https://seaborn.pydata.org/generated/seaborn.set.html>.

"""

import logging
import os

import pandas as pd
import seaborn as sns

import esmvaltool.diag_scripts.emergent_constraints as ec
from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def check_training_data(training_data):
    """Check training data."""
    features = training_data.x
    if len(features.columns) != 1:
        raise ValueError(
            f"Expected exactly 1 'feature' variable, got "
            f"{len(features.columns):d}")


def get_default_settings(cfg):
    """Get default configuration settings."""
    cfg.setdefault('all_data_label', 'all')
    cfg.setdefault('confidence_level', 0.66)
    cfg.setdefault('merge_identical_pred_input', True)
    cfg.setdefault('savefig_kwargs', {
        'bbox_inches': 'tight',
        'dpi': 600,
        'orientation': 'landscape',
    })
    cfg.setdefault('seaborn_settings', {})
    return cfg


def main(cfg):
    """Run the diagnostic."""
    cfg = get_default_settings(cfg)
    sns.set(**cfg['seaborn_settings'])

    # Load data
    (training_data, prediction_data, attributes) = ec.get_input_data(cfg)
    check_training_data(training_data)

    # Plots
    with pd.option_context(*ec.PANDAS_PRINT_OPTIONS):
        logger.info(
            "Correlation of training data (considering all available data):\n"
            "%s", training_data.corr())
    ec.plot_individual_scatterplots(training_data,
                                    prediction_data,
                                    attributes,
                                    'training_data',
                                    cfg)
    ec.plot_merged_scatterplots(training_data,
                                prediction_data,
                                attributes,
                                'training_data',
                                cfg)

    # Export CSV
    ec.export_csv(training_data, attributes, 'training_data', cfg)
    ec.export_csv(prediction_data, attributes, 'prediction_data', cfg)

    # Print constraint
    label = training_data.y.columns[0]
    units = attributes[label]['units']
    constrained_target = ec.get_constraint(training_data, prediction_data,
                                           cfg['confidence_level'])
    logger.info(
        "Constraint on target variable '%s': [%.2f, %.2f] %s with best "
        "estimate %.2f %s", label, constrained_target[0],
        constrained_target[2], units, constrained_target[1], units)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
