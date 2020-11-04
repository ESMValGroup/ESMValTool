"""Plotting boxplots to show weighted ensemble results for ECS."""
import logging
import os
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import requests
import xarray as xr
from _boxplot import boxplot
from climwip import get_diagnostic_filename, log_provenance

from esmvaltool.diag_scripts.shared import get_plot_filename, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def read_weights(filename: str) -> dict:
    """Read a `.nc` file into a weights DataArray."""
    weights_ds = xr.open_dataset(filename)
    return weights_ds.to_dataframe().to_dict()['weight']


def read_metadata(cfg: dict) -> dict:
    """Read the metadata from the configure file."""
    datasets = defaultdict(list)

    metadata = cfg['input_data'].values()

    for item in metadata:
        variable = item['variable_group']

        datasets[variable].append(item)

    return datasets


def reformat_model_strings(weights):
    """Reformat the model strings to compare to model strings in ECS data."""
    models_weights = list(weights.keys())
    new_models = list()
    for model in models_weights:
        # remove experiment from identifiers
        string1 = model.split("_")[0]
        string2 = model.split("_")[1]
        model_id = '%s_%s' % (string1, string2)
        new_models.append(model_id)
        weights[model_id] = weights.pop(model)

    return new_models, weights


def read_ecs(models, cmip):
    """Read the ECS values from file on GitHub.

    Only include models which are present in ECS data and given by
    recipe.
    """
    r_get = requests.get(
        'https://raw.githubusercontent.com/mzelinka/'
        'cmip56_forcing_feedback_ecs/master/cmip56_forcing_feedback_ecs.json')
    data = r_get.json()
    ecs_data = dict()
    models_zelinka = data[cmip].keys()
    for model in models_zelinka:
        ensembles = list(data[cmip][model].keys())
        for ensemble in ensembles:
            ecs_value = data[cmip][model][ensemble]['ECS']
            new_key = model + '_' + ensemble
            ecs_data[new_key] = ecs_value

    models_zelinka_new = set(list(ecs_data.keys()))
    models_intercept = models_zelinka_new.intersection(models)

    ecs_data_return = [ecs_data[member] for member in models_intercept]
    da_ecs = xr.DataArray(ecs_data_return,
                          coords=[list(models_intercept)],
                          dims=["model_ensemble"])
    return da_ecs


def _visualize_and_save_percentiles(data: 'xr.DataArray', weights: dict,
                                    cfg: dict, ancestors: list):

    # ensure weights are sorted the same way as data and convert to list
    weights = [weights[member] for member in data.model_ensemble.values]

    figure, axes = plt.subplots(1, 1, figsize=(4, 10), dpi=300)
    figure.subplots_adjust(left=.1, right=.99, bottom=.22, top=.91)

    box_h1 = boxplot(
        axes,
        0,
        median=data,
        mean=data,
        box=data,
        whisk=data,
        width=.8,
        color='darkgrey',
        alpha=.3,
    )

    box_h2 = boxplot(
        axes,
        0,
        median=data,
        mean=data,
        box=data,
        whisk=data,
        dots=data,
        dots_sizes=(1., 6.),
        weights=weights,
        width=.6,
        color='g',
        alpha=1,
    )

    axes.set_ylabel('Effective Climate Sensitivity')
    axes.grid(axis='y')

    axes.set_xticklabels([])

    plt.ylim(0, 6)
    plt.legend((box_h1, box_h2), ('unweighted', 'weighted'))

    plt.title('Weighted ECS')

    filename_plot = get_plot_filename('boxplot_ecs', cfg)
    figure.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(figure)

    filename_data = get_diagnostic_filename('ecs_percentiles',
                                            cfg,
                                            extension='nc')
    data.to_netcdf(filename_data)

    caption = 'Temperature change boxplot'
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data, cfg, ancestors)


def main(cfg):
    """Plot weighted boxplot."""
    input_files = cfg['input_files']
    filename = cfg['weights']
    weights_path = Path(input_files[0]) / filename
    weights = read_weights(weights_path)

    new_models, new_weights = reformat_model_strings(weights)

    models = read_metadata(cfg)['tas']
    data_files = []
    for info in models:
        filename = info['filename']
        data_files.append(filename)

    any_model = next(iter(cfg['input_data'].values()))
    cmip = any_model['project']
    ecs = read_ecs(new_models, cmip)

    _visualize_and_save_percentiles(
        ecs,
        new_weights,
        cfg,
        data_files,
    )


if __name__ == '__main__':
    with run_diagnostic() as configure:
        main(configure)
