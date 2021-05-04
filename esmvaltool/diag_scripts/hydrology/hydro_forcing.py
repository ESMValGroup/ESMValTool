"""Hydro forcing diagnostic."""
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
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


def plot_data(*, cfg: dict, datasets: dict, xaxis: str, yaxis: str,
              xlabel: str, ylabel: str, caption: str, name: str,
              ancestors: list):
    """Plot data."""
    figure, _ = plt.subplots(dpi=300)

    for label in datasets.dataset:
        label = str(label.data)
        dataset = datasets.sel(dataset=label)
        if 'time' in dataset:
            dataset = dataset.dropna(dim='time')  # remove nan
            figure.autofmt_xdate()  # rotate date labels
        plt.plot(dataset[xaxis], dataset[yaxis], label=label)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(caption)
    plt.legend()
    plt.show()

    filename_plot = get_plot_filename(name + '_plot', cfg)
    figure.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(figure)

    # Store provenance
    log_provenance(caption, filename_plot, cfg, ancestors)


def plot_timeseries(cfg, metadata):
    """Plot timeseries data."""
    short_name = 'pr'
    xaxis = 'time'

    datasets = read_input_data(metadata)
    ancestors = [info['filename'] for info in metadata]

    time_period = cfg['time_period']

    var = datasets[short_name]

    time_unit = time_period[0].upper()
    start_date = np.datetime_as_string(datasets.time.min(), unit=time_unit)
    end_date = np.datetime_as_string(datasets.time.max(), unit=time_unit)

    name = f'{var.long_name}_{time_period}'
    caption = f"{var.long_name} per {time_period} for {start_date}:{end_date}"

    plot_data(
        cfg=cfg,
        datasets=datasets,
        xaxis=xaxis,
        yaxis=short_name,
        xlabel=f'{xaxis.capitalize()} / {time_period}',
        ylabel=f'{var.long_name} / {var.units}',
        caption=caption,
        name=name,
        ancestors=ancestors,
    )

    filename_data = get_diagnostic_filename(name, cfg, extension='nc')
    datasets.to_netcdf(filename_data)
    log_provenance(caption, filename_data, cfg, ancestors)


def plot_climatology(cfg, metadata):
    """Plot climatology data."""
    short_name = 'pr'

    datasets = read_input_data(metadata)
    var = datasets[short_name]

    xaxis = var.dims[-1]  # i.e. month_number / day_of_year
    xlabel = xaxis.replace('_', ' ')
    caption = f'{var.long_name} climatology statistics per {xlabel}'

    ancestors = [info['filename'] for info in metadata]

    name = f'{var.long_name}_climatology_{xaxis}'

    plot_data(
        cfg=cfg,
        datasets=datasets,
        xaxis=xaxis,
        yaxis=short_name,
        xlabel=xlabel.capitalize(),
        ylabel=f'{var.long_name} / {var.units}',
        caption=caption,
        name=name,
        ancestors=ancestors,
    )

    filename_data = get_diagnostic_filename(name, cfg, extension='nc')
    datasets.to_netcdf(filename_data)
    log_provenance(caption, filename_data, cfg, ancestors)


def read_input_data(metadata: list, dim: str = 'dataset'):
    """Load data from metadata.

    Read the input data from the list of given data sets. `metadata` is
    a list of metadata containing the filenames to load. The datasets
    are stacked along the `dim` dimension. Returns an xarray.DataArray.
    """
    identifiers = []
    datasets = []
    for info in metadata:
        filename = info['filename']
        dataset = xr.load_dataset(filename)
        datasets.append(dataset)

        identifier = info[dim]
        identifiers.append(identifier)

    stacked_datasets = xr.concat(datasets, dim=dim)
    stacked_datasets[dim] = identifiers

    return stacked_datasets


def main(cfg):
    """Load and plot hydro forcing data."""
    plot_type = cfg['plot_type']

    input_data = cfg['input_data'].values()
    variable_groups = group_metadata(input_data, 'variable_group')

    plot_func_mapping = {
        'climatology': plot_climatology,
        'timeseries': plot_timeseries,
    }

    for metadata in variable_groups.values():
        try:
            plot_func = plot_func_mapping[plot_type]
        except KeyError as err:
            raise ValueError(f'Unknown plot_type: {plot_type!r}') from err

        plot_func(cfg, metadata=metadata)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
