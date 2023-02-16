"""Python example diagnostic."""
import logging
import os
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
from iris.analysis.stats import pearsonr
import iris.plot as iplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import seaborn as sns
from scipy.stats import bootstrap

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import (
    extract_variables,
    group_metadata,
    run_diagnostic,
    get_diagnostic_filename,
    get_plot_filename,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
    io,
)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(Path(__file__).stem)

VAR_NAMES = {
    'clt': 'Total Cloud Fraction',
    'lwp': 'Liquid Water Path',
    'clivi': 'Ice Water Path',
    'netcre': 'net_cre',
    'swcre': 'sw_cre',
    'lwcre': 'lw_cre',
}
PALETTE = {
    'ECS_high': 'royalblue',
    'ECS_med': 'green',
    'ECS_low': 'orange',
}


def create_data_frame(input_files):
    """Collect all correlation values from the input_files."""

    corr_df = None
    
    for ifile in input_files:
        df = pd.read_csv(ifile)
        select_corr = df.loc[(df['Statistic'] == 'Corr') & (df['Group'] != 'OBS')] 
        var = ifile.split("_")[-1].split(".")[0]
        select_corr['Variable'] = var
        print(var)
        if corr_df is None:
            corr_df = select_corr
        else:
            corr_df = pd.concat([corr_df, select_corr], axis = 0)
            print("Append")
        print(corr_df.tail())

    return corr_df


def plot_corr(data_frame, cfg):
    """Plot correlation figure"""

    plt.figure(constrained_layout=True, figsize=(12, 8))

    sns.set_style('darkgrid')
    sns.set(font_scale=2)
    #sns.boxplot(data=data_frame, x='Variable', y='Value', hue='Group', palette=PALETTE)
    sns.stripplot(data=data_frame, x='Variable', y='Value', 
                  order = cfg['diag_order'],
                  marker = '_', s = 20, jitter=False, 
                  hue='Group', hue_order = ['ECS_low', 'ECS_med', 'ECS_high'],
                  dodge=True, palette=PALETTE)
    plt.ylabel('Correlation')

    # Save plot
    plot_path = get_plot_filename('pattern_corr', cfg)
    plt.savefig(plot_path)
    logger.info("Wrote %s", plot_path)
    plt.close()



def main(cfg):
    """Run diagnostic."""

    # Get input files
    pattern = cfg.get('patterns')
    
    input_files = [i for i in io.get_all_ancestor_files(cfg) if pattern in i]
    if not input_files:
        raise ValueError("No input files found")
    logger.info("Found input files:\n%s", pformat(input_files))

    # Create data frame
    data_frame = create_data_frame(input_files)

    # Plot correlation
    plot_corr(data_frame, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
