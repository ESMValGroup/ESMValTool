#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generate config files which can be used as input for tests."""

import os
import random
from unittest import mock

import yaml

from tests.integration.diag_scripts.mlr.test_read_input import (
    SimplifiedMLRModel,
    get_logger_msg,
)

FUNCTION = '_load_input_datasets'
N_CFG = 20
N_MAX_DATA = 10
OUTFILE = os.path.expanduser(os.path.join('~', 'outfile.yml'))
DATASET = {
    'dataset': ['dataset_1', 'dataset_2', 'dataset_3'],
    'project': ['project_1', 'project_2', 'project_3'],
    'standard_name': ['std_name_1', 'std_name_2', 'std_name_3'],
    'short_name': ['short_name_1', 'short_name_2', 'short_name_3'],
    'long_name': ['long_name_1', 'long_name_2', 'long_name_3'],
    'units': ['units_1', 'units_2', 'units_3'],
    'filename': ['/path/1', '/path/2', 'path/3'],
    'var_type': ['feature', 'label', 'prediction_input', 'wrong_var_type'],
    'tag': ['tag_1', 'tag_2', 'tag_3'],
    'prediction_name': [None, 'pred_name_1', 'pred_name_2'],
    'broadcast_from': [None, None, None, None, 0, [0], [0, 1], [4, 5]],
}
CFG = {
    'input_files': [[], ['/ancestor/1'], ['/ancestor/1', '/ancestor/2']],
    'output_file_type': ['png', 'pdf', 'ps'],
    'plot_dir': ['/plot/dir/1', '/plot/dir/2'],
    'work_dir': ['/work/dir/1', '/work/dir_2'],
    'write_plots': [True, False],
    'accept_only_scalar_data': [None, True, False],
    'allow_missing_features': [None, True, False],
    'grid_search_cv_param_grid': [
        None,
        None,
        {
            'a': [1, 2, 3]
        },
        [{
            'a': [1, 2],
            'b': [3.14, 2.71]
        }, {
            'b': [6.28, -1]
        }],
    ],
    'group_datasets_by_attributes': [
        None,
        None,
        ['dataset'],
        ['dataset', 'project'],
    ],
    'imputation_strategy': [None, None, 'remove', 'mean'],
    'mlr_model_type': [None, None, 'gbr', 'gpr'],
    'parameters': [
        None,
        None,
        {
            'a': 1
        },
        {
            'a': 1,
            'b': 3.1415
        },
    ],
    'predict_kwargs': [None, None, {
        'return_var': True
    }],
    'test_size': [None, None, 0.25, -1.0, 2.5],
    'coords_as_features': [
        None,
        None,
        None,
        ['latitude'],
        ['air_pressure', 'latitude'],
    ],
}

random.seed(42)


def generate_random_dict(source, remove_prob=0.0):
    """Generate random dict_ using `source`."""
    dict_ = {}
    for (attr, values) in source.items():
        value = values[random.randrange(len(values))]
        if value is not None:
            dict_[attr] = value
    if random.random() < remove_prob:
        rand_attr = list(dict_.keys())[random.randrange(len(dict_))]
        dict_.pop(rand_attr)
    return dict_


if __name__ == '__main__':
    CFGS = []

    # Generate data
    for idx_cfg in range(N_CFG):
        key_cfg = 'cfg_{:d}'.format(idx_cfg)
        cfg = generate_random_dict(CFG)
        datasets = {}
        for idx_data in range(random.randrange(N_MAX_DATA)):
            key_data = 'dataset_{:d}'.format(idx_data)
            dataset = generate_random_dict(DATASET, remove_prob=0.1)
            datasets[key_data] = dataset
        cfg['input_data'] = datasets
        input_datasets = list(datasets.values())

        # Output
        mlr_model = SimplifiedMLRModel(cfg)
        logger_calls = []
        with mock.patch('esmvaltool.diag_scripts.mlr.logger',
                        autospec=True) as mlr_logger:
            with mock.patch('esmvaltool.diag_scripts.mlr.models.logger',
                            autospec=True) as models_logger:
                try:
                    getattr(mlr_model, FUNCTION)(input_datasets)
                except Exception as exc:
                    output = {'EXCEPTION': {}}
                    output['EXCEPTION']['type'] = type(exc).__name__
                    output['EXCEPTION']['value'] = exc.args[0]
                else:
                    output = mlr_model._datasets
                logger_calls.extend(models_logger.method_calls)
            logger_calls.extend(mlr_logger.method_calls)
        CFGS.append({
            'cfg': cfg,
            'output': output,
            'logger': get_logger_msg(logger_calls),
        })

    # Write data
    NoAnchorDumper = yaml.dumper.SafeDumper
    NoAnchorDumper.ignore_aliases = lambda self, data: True
    with open(OUTFILE, 'w') as outfile:
        yaml.dump(CFGS,
                  outfile,
                  default_flow_style=False,
                  Dumper=NoAnchorDumper)
        print("Wrote '{}'".format(OUTFILE))
