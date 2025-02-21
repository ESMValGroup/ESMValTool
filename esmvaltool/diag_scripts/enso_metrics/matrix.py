"""Diagnostic script to plot ENSO metrics portrait matrix."""

import os
import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import (run_diagnostic,
                                            save_figure)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def plot_matrix(diag_path):
    """Read output values from csv and plot."""
    metric_df = pd.read_csv(diag_path, header=None)
    # run normalisation on all these values
    # metric_df[2] = (metric_df[2] - metric_df[2].mean()) / metric_df[2].std()

    t_list = []
    for mod in metric_df[0].unique():  # iterate model, translate metrics
        mod_df = metric_df.loc[metric_df[0] == mod, :]
        t_list.append(mod_df[[1, 2]].set_index(1).T.rename(index={2: mod}))

    matrixdf = pd.concat(t_list)
    #normalise column by column
    for col in matrixdf.columns:
        matrixdf[col] = (matrixdf[col] - matrixdf[col].mean()) / matrixdf[col].std()

    figure = plt.figure(dpi=300)
    plt.imshow(matrixdf, cmap='coolwarm')
    plt.colorbar()
    plt.xticks(range(len(matrixdf.columns)), matrixdf.columns,
               rotation=45, ha='right')
    plt.yticks(range(len(matrixdf.index)), matrixdf.index, wrap=True)
    plt.xticks(np.arange(matrixdf.shape[1] + 1) - 0.5, minor=True)
    plt.yticks(np.arange(matrixdf.shape[0] + 1) - 0.5, minor=True)
    plt.tick_params(which="both", bottom=False, left=False)
    plt.grid(which="minor", color="black", linestyle="-", linewidth=0.5)

    return figure


def main(cfg):
    """Read metrics and plot matrix."""
    provenance_record = {
        'caption': "ENSO metrics",
        'authors': [
            'chun_felicity',
        ],
        'references': [''],
        'ancestors': cfg['diag_metrics']
    }

    metrics = cfg['diag_metrics']
    diag_path = '/'.join(cfg['work_dir'].split('/')[:-2])
    diag_path = '/'.join([diag_path, metrics, 'matrix.csv'])
    logger.info(diag_path)

    figure = plot_matrix(diag_path)

    save_figure('plot_matrix', provenance_record, cfg,
                figure=figure, bbox_inches='tight')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
