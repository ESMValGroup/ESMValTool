"""Hydro forcing diagnostic."""
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(caption: str, ancestors: list):
    """Create a provenance record describing the diagnostic data and plots."""
    record = {
        'caption': caption,
        'domains': ['global'],
        'authors': [
            'smeets_stef',
            'aerts_jerom',
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestors,
    }
    return record


def log_provenance(caption: str, filename: str, cfg: dict, ancestors: list):
    """Log provenance info."""
    provenance_record = get_provenance_record(caption, ancestors=ancestors)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    logger.info('Output stored as %s', filename)


def plot_data(
    *,
    cfg: dict,
    datasets: dict,
    xaxis: str,
    yaxis: str,
    xlabel: str,
    ylabel: str,
    caption: str,
    filename: str,
    ancestors: list,
):
    """Plot data."""
    figure, axes = plt.subplots(dpi=300)

    for label, dataset in datasets.items():
        plt.plot(dataset[xaxis], dataset[yaxis], label=label)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(caption)
    plt.legend()
    plt.show()

    filename_plot = get_plot_filename(filename, cfg)
    figure.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(figure)

    # Store provenance
    log_provenance(caption, filename_plot, cfg, ancestors)


def plot_timeseries_data(cfg, time_unit: str = 'D', title: str = 'plot'):
    """Plot timeseries data."""
    input_data = cfg['input_data']

    datasets = {}
    for filename, metadata in input_data.items():
        alias = metadata['alias']
        dataset = xr.load_dataset(filename)
        datasets[alias] = dataset

    start_date = np.datetime_as_string(dataset.time.min(), unit=time_unit)
    end_date = np.datetime_as_string(dataset.time.max(), unit=time_unit)

    plot_data(
        cfg=cfg,
        datasets=datasets,
        xaxis='time',
        yaxis='pr',
        xlabel='time / {mip}'.format(**metadata),
        ylabel='{long_name} / {units}'.format(**metadata),
        caption=f"{title} for {start_date}:{end_date}",
        filename=title.lower().replace(' ', '_') + '_plot',
        ancestors=list(input_data.keys()),
    )


def plot_climatology(cfg, title: str = 'plot'):
    """Plot climatology data."""
    input_data = cfg['input_data']

    datasets = {}
    for filename, metadata in input_data.items():
        alias = metadata['alias']
        datasets[alias] = xr.load_dataset(filename)

    plot_data(
        cfg=cfg,
        datasets=datasets,
        xaxis='month_number',
        yaxis='pr',
        xlabel='Month number',
        ylabel='{long_name} / {units}'.format(**metadata),
        caption=f"{title} per month",
        filename=title.lower().replace(' ', '_') + '_plot',
        ancestors=list(input_data.keys()),
    )


def main(cfg):
    """Main function."""
    script = cfg['script']

    if script == 'plot_sample_year':
        plot_timeseries_data(
            cfg,
            time_unit='D',
            title='Daily precipitation',
        )
    elif script == 'plot_total_precip':
        plot_timeseries_data(
            cfg,
            time_unit='M',
            title='Monthly total precipitation',
        )
    elif script == 'plot_climatology':
        plot_climatology(
            cfg,
            title='Precipitation per month',
        )
    else:
        raise ValueError(f'Unknown script: {script}')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
