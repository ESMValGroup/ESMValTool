#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to evaluate multiple emergent constraints simultaneously.

Description
-----------
Establish multiple emergent constraints for arbitrary input variables and an
arbitrary target variable. All input datasets need to be one-dimensional and
must include a coordinate ``'dataset'`` or ``'model'`` (thus, the data
describes a single scalar value for each dataset). All input datasets must be
marked with a ``var_type`` (either ``feature``, ``label``, ``prediction_input``
or ``prediction_input_error``) and a ``tag``, which describes the type of data.
This diagnostic supports only a single ``tag`` for ``label`` and an arbitrary
number of ``tag`` s for ``feature``. For every ``tag``, a
``'reference_dataset'`` can be specified, which will be automatically
considered as ``prediction_input``.  If ``reference_dataset`` contains ``'|'``
(e.g. ``'OBS1|OBS2'``), multiple datasets are considered as
``prediction_input`` (in this case ``'OBS1'`` and ``'OBS2'``).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
additional_data: list of dict, optional
    Additional datasets given as list of metadata.
all_data_label: str, optional (default: 'all')
    Label used in plots when all input data is considered. Only relevant if
    ``group_by`` is not used.
combine_groups: bool, optional (default: False)
    Add results to plots for data generated by combining the data of all
    individual groups.
confidence_level: float, optional (default: 0.66)
    Confidence level for estimation of target variable.
group_by: str, optional
    Group input data by an attribute (e.g. produces separate plots for the
    individual groups, etc.).
ignore_patterns: list of str, optional
    Patterns matched against ancestor files. Those files are ignored.
merge_identical_pred_input: bool, optional (default: True)
    Use identical prediction_input values as single value.
numbers_as_markers: bool, optional (default: False)
    Use numbers as markers in scatterplots.
patterns: list of str, optional
    Patterns matched against ancestor files.
plot_regression_line_mean: bool, optional (default: False)
    Plot means of regression lines in scatterplots.
read_external_file: str, optional
    Read input datasets from external file given as absolute path or relative
    path. In the latter case, ``'auxiliary_data_dir'`` from the user
    configuration file is used as base directory.
savefig_kwargs: dict
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings: dict
    Options for :func:`seaborn.set` (affects all plots), see
    `<https://seaborn.pydata.org/generated/seaborn.set.html>`_.

"""

import logging
import os
from copy import deepcopy

import pandas as pd
import seaborn as sns

import esmvaltool.diag_scripts.emergent_constraints as ec
from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def get_default_settings(cfg):
    """Get default configuration settings."""
    cfg = deepcopy(cfg)
    cfg.setdefault('all_data_label', 'all')
    cfg.setdefault('combine_groups', False)
    cfg.setdefault('confidence_level', 0.66)
    cfg.setdefault('merge_identical_pred_input', True)
    cfg.setdefault('patterns', [])
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

    # Load data and perform PCA
    (training_data, prediction_data, attributes) = ec.get_input_data(cfg)
    training_data_no_nans = training_data.dropna()

    # Plots
    with pd.option_context(*ec.PANDAS_PRINT_OPTIONS):
        logger.info(
            "Correlation of training data (considering all available data):\n"
            "%s", training_data.corr())
        logger.info(
            "Correlation of training data (considering only climate models "
            "where data for all constraints is available):\n%s",
            training_data_no_nans.corr())
    ec.plot_individual_scatterplots(training_data,
                                    prediction_data,
                                    attributes,
                                    'training_data',
                                    cfg)
    ec.plot_merged_scatterplots(training_data, prediction_data, attributes,
                                'training_data', cfg)
    ec.plot_target_distributions(training_data, prediction_data, attributes,
                                 'training_data', cfg)

    # Export CSV
    ec.export_csv(training_data, attributes, 'training_data', cfg)
    ec.export_csv(training_data_no_nans, attributes, 'training_data_no_nans',
                  cfg)
    ec.export_csv(prediction_data, attributes, 'prediction_data', cfg)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
