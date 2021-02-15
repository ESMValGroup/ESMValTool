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


def plot_data(
    *,
    cfg: dict,
    datasets: dict,
    xaxis: str,
    yaxis: str,
    xlabel: str,
    ylabel: str,
    caption: str,
    name: str,
    ancestors: list,
):
    """Plot data."""
    figure, axes = plt.subplots(dpi=300)

    for label in datasets.dataset:
        label = str(label.data)
        dataset = datasets.sel(dataset=label)
        if 'time' in dataset:
            dataset = dataset.dropna(dim='time')  # remove nan
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


def plot_timeseries_data(cfg,
                         metadata,
                         time_unit: str = 'D',
                         title: str = 'plot'):
    """Plot timeseries data."""
    datasets = read_input_data(metadata)
    ancestors = [info['filename'] for info in metadata]

    var = datasets.pr

    start_date = np.datetime_as_string(datasets.time.min(), unit=time_unit)
    end_date = np.datetime_as_string(datasets.time.max(), unit=time_unit)

    name = title.lower().replace(' ', '_')
    caption = f"{title} for {start_date}:{end_date}"

    plot_data(
        cfg=cfg,
        datasets=datasets,
        xaxis='time',
        yaxis='pr',
        xlabel='time / {mip}'.format(**datasets.attrs),
        ylabel=f'{var.long_name} / {var.units}',
        caption=caption,
        name=name,
        ancestors=ancestors,
    )

    filename_data = get_diagnostic_filename(name, cfg, extension='nc')
    datasets.to_netcdf(filename_data)
    log_provenance(caption, filename_data, cfg, ancestors)


def plot_climatology(cfg, metadata, title: str = 'plot'):
    """Plot climatology data."""
    datasets = read_input_data(metadata)
    var = datasets.pr

    ancestors = [info['filename'] for info in metadata]

    name = title.lower().replace(' ', '_')
    caption = f"{title}"

    plot_data(
        cfg=cfg,
        datasets=datasets,
        xaxis='month_number',
        yaxis='pr',
        xlabel='Month number',
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


def get_diagnostic_data(cfg):
    """Get diagnostic name and data."""
    input_data = cfg['input_data'].values()
    diagnostic_data = group_metadata(input_data, 'diagnostic')
    return list(diagnostic_data.items())[0]


def main(cfg):
    """Main function."""
    diagnostic, metadata = get_diagnostic_data(cfg)

    if diagnostic == 'sample_year':
        plot_timeseries_data(
            cfg,
            metadata=metadata,
            time_unit='D',
            title='Daily precipitation',
        )
    elif diagnostic == 'total_precipitation':
        plot_timeseries_data(
            cfg,
            metadata=metadata,
            time_unit='M',
            title='Monthly total precipitation',
        )
    elif diagnostic == 'climatology':
        plot_climatology(
            cfg,
            metadata=metadata,
            title='Precipitation per month',
        )
    else:
        raise ValueError(f'Unknown diagnostic: {diagnostic}')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
