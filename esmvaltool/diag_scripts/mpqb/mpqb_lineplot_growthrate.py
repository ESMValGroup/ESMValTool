#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os
from pprint import pformat

import iris
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from mpqb_utils import get_mpqb_cfg

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger,
    get_plot_filename,
)

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(caption):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['diff'],
        'domains': ['user'],
        'plot_type': 'lineplot',
        'authors': [
            'mueller_benjamin',
            'crezee_bas',
            'hassler_birgit',
        ],
        'projects': ['cmug'],
        'references': [
            'acknow_project',
        ],
    }
    return record


def main(cfg):
    """Create lineplot."""
    # Limits for the y-axis
    ylims = [-5., 15.]

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(input_data, 'alias', sort='alias')

    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))

    plt.clf()
    fig, (ax1, lax) = plt.subplots(nrows=2,
                                   gridspec_kw={"height_ratios": [10, 1]},
                                   figsize=(10, 5))

    plt.sca(ax1)
    for dataset in grouped_input_data:
        dataset_cfg = grouped_input_data[dataset][0]
        alias = dataset_cfg['alias']

        logger.info("Opening dataset: %s", dataset)
        cube = iris.load_cube(dataset_cfg['filename'])
        if cube.coords('time', dim_coords=True):
            ih.unify_time_coord(cube)

        # calculate the growth rates
        diff_data = np.diff(cube.data)
        cube = cube[1:]
        cube.data = diff_data

        # Set default if not defined.
        label = get_mpqb_cfg('datasetname', alias)
        color = get_mpqb_cfg('datasetcolor', alias)

        # iris.quickplot.plot(cube, label=label, color=color,
        #                     linestyle='dotted')
        # iris.quickplot.plot(cube, label=label, color=color,
        #                     linestyle='dashed')
        iris.quickplot.plot(cube, label=label, color=color)
    plt.xticks(rotation=90)
    # Add the zero line
    plt.axhline(y=0, linestyle=':', color='k')
    plt.tight_layout()

    # Time axis formatting
    years = mdates.YearLocator()  # every year
    years_fmt = mdates.DateFormatter('%Y')
    ax1 = plt.gca()
    ax1.xaxis.set_major_locator(years)
    ax1.xaxis.set_major_formatter(years_fmt)
    ax1.grid(True, which='major', axis='x')

    ax1.set_ylim(ylims)
    ax1.set_ylabel(f"Change in {cube.var_name.upper()} ({cube.units})")
    ax1.set_title(f"Time series of annual {cube.var_name.upper()} growth rate")

    h1, l1 = ax1.get_legend_handles_labels()
    leg = lax.legend(h1, l1, borderaxespad=0, ncol=4, loc='center')
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    lax.axis("off")

    baseplotname = (f"lineplot_{dataset_cfg['variable_group']}_"
                    f"{dataset_cfg['start_year']}-{dataset_cfg['end_year']}")

    filename = get_plot_filename(baseplotname, cfg)
    logger.info("Saving as %s", filename)
    fig.savefig(filename, bbox_inches='tight')

    caption = (
        f"Domain average annual growth rates of "
        f"{dataset_cfg['variable_group']} "
        f"between {dataset_cfg['start_year']} and {dataset_cfg['end_year']}")

    provenance_record = get_provenance_record(caption)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    plt.close(fig)
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
