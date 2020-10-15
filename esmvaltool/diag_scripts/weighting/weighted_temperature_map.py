"""Implementation of a mapplot for the climwip weighting scheme.

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
from collections import defaultdict
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import numpy as np
import xarray as xr
from climwip import log_provenance, read_model_data
from climwip import get_diagnostic_filename

from esmvaltool.diag_scripts.shared import get_plot_filename, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def read_weights(filename: str) -> dict:
    """Read a `.nc` file into a weights DataArray."""
    weights_ds = xr.open_dataset(filename)
    return weights_ds['weight']


def read_metadata(cfg: dict) -> dict:
    """Read the metadata from the config file."""
    datasets = defaultdict(list)

    metadata = cfg['input_data'].values()

    for item in metadata:
        variable = item['variable_group']

        datasets[variable].append(item)

    return datasets


def weighted_quantile(values: list,
                      quantile: float,
                      weights: list = None) -> 'np.array':
    """Calculate weighted quantiles.

    Analogous to np.quantile, but supports weights.

    Based on: https://stackoverflow.com/a/29677616/6012085

    Parameters
    ----------
    values: array_like
        List of input values.
    quantile: array_like
        List of quantiles between 0.0 and 1.0.
    weights: array_like
        List with same length as `values` containing the weights.

    Returns
    -------
    np.array
        Numpy array with computed quantile.
    """
    values = np.array(values)
    if weights is None:
        weights = np.ones(len(values))
    weights = np.array(weights)

    if not ((quantile >= 0) & (quantile <= 1)):
        raise ValueError('Quantile should be between 0.0 and 1.0')

    idx = np.argsort(values)
    values = values[idx]
    weights = weights[idx]

    weighted_quantiles = np.cumsum(weights) - 0.5 * weights

    # Cast weighted quantiles to 0-1 To be consistent with np.quantile
    min_val = weighted_quantiles.min()
    max_val = weighted_quantiles.max()
    weighted_quantiles = (weighted_quantiles - min_val) / max_val

    return np.interp(quantile, weighted_quantiles, values)


def visualize_and_save_temperature(temperature: 'xr.DataArray',
                                   cfg: dict,
                                   ancestors: list):
    """Visualize weighted temperature."""
    period = '{start_year}-{end_year}'.format(**read_metadata(cfg)['tas'][0])
    if (meta := read_metadata(cfg).get('tas_reference', None)) is not None:
        period = '{} minus {start_year}-{end_year}'.format(period, **meta[0])
    metric = cfg['model_aggregation']
    if isinstance(metric, int):
        metric = f'{metric}perc'

    proj = ccrs.PlateCarree(central_longitude=0)
    fig, ax = plt.subplots(subplot_kw={'projection': proj})

    temperature.plot.pcolormesh(
        ax=ax,
        transform=ccrs.PlateCarree(),
        levels=9,
        robust=True,
        extend='both',
        cmap='Reds',
    )

    ax.coastlines()
    lons = temperature.lon.values
    lats = temperature.lat.values
    longitude_formatter = LongitudeFormatter()
    latitude_formatter = LatitudeFormatter()
    if (xticks := cfg.get('xticks', None)) is not None:
        ax.set_xticks(xticks)
    else:
        ax.set_xticks(np.arange(np.floor(lons.min()), np.ceil(lons.max()), 10), crs=proj)
    if (yticks := cfg.get('yticks', None)) is not None:
        ax.set_yticks(yticks)
    else:
        ax.set_yticks(np.arange(np.floor(lats.min()), np.ceil(lats.max()), 10), crs=proj)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_major_formatter(longitude_formatter)
    ax.yaxis.set_major_formatter(latitude_formatter)
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax.set_title(f'Weighted {metric} temperature \n{period} ($\degree$C)')

    filename_plot = get_plot_filename('temperature_change_weighted_map', cfg)
    fig.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(fig)

    filename_data = get_diagnostic_filename(
        'temperature_map_weighted', cfg, extension='nc')
    temperature.to_netcdf(filename_data)

    caption = f'Temperature map {period}'
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data, cfg, ancestors)


def visualize_and_save_temperature_difference(temperature: 'xr.DataArray',
                                              cfg: dict,
                                              ancestors: list):
    """Visualize weighted temperature."""
    period = '{start_year}-{end_year}'.format(**read_metadata(cfg)['tas'][0])
    if (meta := read_metadata(cfg).get('tas_reference', None)) is not None:
        period = '{} minus {start_year}-{end_year}'.format(period, **meta[0])
    metric = cfg['model_aggregation']
    if isinstance(metric, int):
        metric = f'{metric}perc'

    proj = ccrs.PlateCarree(central_longitude=0)
    fig, ax = plt.subplots(subplot_kw={'projection': proj})

    temperature.plot.pcolormesh(
        ax=ax,
        transform=ccrs.PlateCarree(),
        levels=9,
        center=0,
        robust=True,
        extend='both',
    )

    ax.coastlines()
    lons = temperature.lon.values
    lats = temperature.lat.values
    longitude_formatter = LongitudeFormatter()
    latitude_formatter = LatitudeFormatter()
    if (xticks := cfg.get('xticks', None)) is not None:
        ax.set_xticks(xticks)
    else:
        ax.set_xticks(np.arange(np.floor(lons.min()), np.ceil(lons.max()), 10), crs=proj)
    if (yticks := cfg.get('yticks', None)) is not None:
        ax.set_yticks(yticks)
    else:
        ax.set_yticks(np.arange(np.floor(lats.min()), np.ceil(lats.max()), 10), crs=proj)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_major_formatter(longitude_formatter)
    ax.yaxis.set_major_formatter(latitude_formatter)
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax.set_title(f'Weighted - unweighted {metric} temperature\n{period} ($\degree$C)')

    filename_plot = get_plot_filename('temperature_map_difference', cfg)
    fig.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(fig)

    filename_data = get_diagnostic_filename(
        'temperature_map_difference', cfg, extension='nc')
    temperature.to_netcdf(filename_data)

    caption = 'Temperature map difference'
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data, cfg, ancestors)



def calculate_percentiles(data: 'xr.DataArray',
                          percentile: float,
                          weights: dict = None):
    """Calculate (weighted) percentiles.

    Calculate the (weighted) percentiles for the given data.

    Percentiles is a list of values between 0 and 100.

    Weights is a dictionary where the keys represent the model members,
    and the values the weights.
    If `weights` is not specified, the non-weighted percentiles are calculated.

    Returns a DataArray with 'percentiles' as the dimension.
    """
    if weights is not None:
        weights = weights.sel(model_ensemble=data.model_ensemble)

    output = xr.apply_ufunc(weighted_quantile,
                            data,
                            input_core_dims=[['model_ensemble']],
                            kwargs={
                                'weights': weights,
                                'quantile': percentile / 100
                            },
                            vectorize=True)

    return output


def model_aggregation(dataset, metric, weights=None):
    if isinstance(metric, int):
        return calculate_percentiles(dataset, metric, weights)
    elif metric.lower() == 'mean':
        if weights is not None:
            dataset = dataset.weighted(weights)
        return dataset.mean('model_ensemble')
    elif metric.lower() == 'median':
        return calculate_percentiles(dataset, 50, weights)
    else:
        errmsg = f'model_aggregation {metric} is not implemented!'
        raise NotImplementedError(errmsg)


def main(cfg):
    """Plot weighted temperature graph."""
    input_files = cfg['input_files']
    filename = cfg['weights']
    weights_path = Path(input_files[0]) / filename
    weights = read_weights(weights_path)

    models = read_metadata(cfg)['tas']
    model_data, model_data_files = read_model_data(models)

    # if a historical period is given calculate the change
    if (models_hist := read_metadata(cfg).get('tas_reference', None)) is not None:
        model_data_hist, model_data_files_hist = read_model_data(models_hist)
        model_data_files += model_data_files_hist
        model_data = model_data - model_data_hist

    metric = cfg.get('model_aggregation', 'mean')
    unweighted_mean = model_aggregation(model_data, metric)
    weighted_mean = model_aggregation(model_data, metric, weights)

    visualize_and_save_temperature(
        weighted_mean,
        cfg,
        model_data_files,
    )

    visualize_and_save_temperature_difference(
        weighted_mean - unweighted_mean,
        cfg,
        model_data_files,
    )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
