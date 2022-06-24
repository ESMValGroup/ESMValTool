#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os
from pprint import pformat

import iris
import matplotlib.pyplot as plt
from mpqb_utils import get_mpqb_cfg

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
        'statistics': ['clim'],
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
    ylims = [cfg.pop('y0', None), cfg.pop('y1', None)]

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(input_data, 'alias', sort='alias')

    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))

    # In order to get the right line colors for MPQB soil moisture
    # here we put ERA-Interim-Land at the end of the dictionary if
    # it is included.
    if 'ERA-Interim-Land' in grouped_input_data.keys():
        grouped_input_data.move_to_end('ERA-Interim-Land')

    plt.clf()
    # fig = plt.figure(figsize=(10, 4))
    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot()

    for dataset in grouped_input_data:
        dataset_cfg = grouped_input_data[dataset][0]
        alias = dataset_cfg['alias']

        logger.info("Opening dataset: %s", dataset)
        cube = iris.load_cube(dataset_cfg['filename'])

        # Set default if not defined.
        label = get_mpqb_cfg('datasetname', alias)
        color = get_mpqb_cfg('datasetcolor', alias)

        # iris.quickplot.plot(cube, label=label, color=color,
        #                     linestyle='dotted')
        # iris.quickplot.plot(cube, label=label, color=color,
        #                     linestyle='dashed')
        iris.quickplot.plot(cube, label=label, color=color)
    plt.legend()
    plt.tight_layout()
    # Time axis formatting
    ax1 = plt.gca()
    ax1.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    ax1.set_xticklabels(['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                         'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'])
    ax1.set_xlabel('month')
    ax1.grid(True, which='major', axis='x')

    ax1.set_ylim(ylims)
    ax1.set_ylabel(f"{cube.var_name.upper()} ({cube.units})")
    ax1.set_title(f"{cube.var_name.upper()} annual cycle")

    baseplotname = (f"lineplot_{dataset_cfg['variable_group']}_"
                    f"{dataset_cfg['start_year']}-{dataset_cfg['end_year']}")

    filename = get_plot_filename(baseplotname, cfg)
    logger.info("Saving as %s", filename)
    fig.savefig(filename, bbox_inches='tight')

    # Provenance
    caption = (
        f"Domain average annual cycle of {dataset_cfg['variable_group']} "
        f"between {dataset_cfg['start_year']} and {dataset_cfg['end_year']}")

    provenance_record = get_provenance_record(caption)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    plt.close(fig)
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
