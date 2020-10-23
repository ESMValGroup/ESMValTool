"""Implementation of a mapplot for the climwip weighting scheme.

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
import logging
import os
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

from climwip import (
    get_diagnostic_filename,
    get_plot_filename,
    log_provenance,
    read_model_data,
)
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.weighting.plot_utilities import (
    calculate_percentiles,
    read_metadata,
    read_weights,
)

logger = logging.getLogger(os.path.basename(__file__))


def mapplot(dataarray, cfg, title_pattern, filename_part, ancestors, **colormesh_args):
    """Visualize weighted temperature."""
    period = '{start_year}-{end_year}'.format(**read_metadata(cfg)['tas'][0])
    if 'tas_reference' in read_metadata(cfg).keys():
        meta = read_metadata(cfg)['tas_reference']
        period = 'change: {} minus {start_year}-{end_year}'.format(
            period, **meta[0])
    metric = cfg['model_aggregation']
    if isinstance(metric, int):
        metric = f'{metric}perc'
    proj = ccrs.PlateCarree(central_longitude=0)
    fig, ax = plt.subplots(subplot_kw={'projection': proj})

    dataarray.plot.pcolormesh(
        ax=ax,
        transform=ccrs.PlateCarree(),
        levels=9,
        robust=True,
        extend='both',
        **colormesh_args
        # colorbar size often does not fit nicely
        # https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
        # cbar_kwargs={'fraction': .021}
    )

    lons = dataarray.lon.values
    lats = dataarray.lat.values
    longitude_formatter = LongitudeFormatter()
    latitude_formatter = LatitudeFormatter()
    default_xticks = np.arange(np.floor(lons.min()), np.ceil(lons.max()), 10)
    default_yticks = np.arange(np.floor(lats.min()), np.ceil(lats.max()), 10)

    ax.coastlines()
    ax.set_xticks(cfg.get('xticks', default_xticks), crs=proj)
    ax.set_yticks(cfg.get('yticks', default_yticks), crs=proj)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_major_formatter(longitude_formatter)
    ax.yaxis.set_major_formatter(latitude_formatter)
    ax.set_xlabel('')
    ax.set_ylabel('')

    title = title_pattern.format(metric=metric, period=period)
    ax.set_title(title)

    filename_plot = get_plot_filename(filename_part, cfg)
    fig.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(fig)

    filename_data = get_diagnostic_filename(filename_part, cfg, extension='nc')
    dataarray.to_netcdf(filename_data)

    log_provenance(title, filename_plot, cfg, ancestors)
    log_provenance(title, filename_data, cfg, ancestors)


def visualize_and_save_temperature(temperature: 'xr.DataArray', cfg: dict,
                                   ancestors: list):
    """Wrapper for mapplot: absolute temperature."""
    title_pattern = 'Weighted {metric} temperature \n{period} ($\degree$C)'
    filename_part = 'temperature_change_weighted_map'
    mapplot(temperature,
            cfg,
            title_pattern,
            filename_part,
            ancestors,
            cmap='Reds')


def visualize_and_save_temperature_difference(
        temperature_difference: 'xr.DataArray', cfg: dict, ancestors: list):
    """Wrapper for mapplot: temperature difference between weighted and unweighted."""
    title_pattern = 'Difference: weighted minus unweighted {metric} temperature\n{period} ($\degree$C)'
    filename_part = 'temperature_change_difference_map'
    mapplot(temperature_difference,
            cfg,
            title_pattern,
            filename_part,
            ancestors,
            center=0)


def model_aggregation(dataset, metric, weights=None):
    """Call mean or percentile calculation."""
    if isinstance(metric, int):
        return calculate_percentiles(dataset, [metric], weights).squeeze(
            'percentile', drop=True)
    elif metric.lower() == 'mean':
        if weights is not None:
            dataset = dataset.weighted(weights)
        return dataset.mean('model_ensemble')
    elif metric.lower() == 'median':
        return calculate_percentiles(dataset, [50], weights).squeeze(
            'percentile', drop=True)
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
    if 'tas_reference' in read_metadata(cfg).keys():
        models_hist = read_metadata(cfg)['tas_reference']
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
