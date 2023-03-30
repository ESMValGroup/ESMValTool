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
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from climwip.io_functions import log_provenance, read_model_data

from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)
from esmvaltool.diag_scripts.weighting.plot_utilities import (
    calculate_percentiles,
    read_metadata,
    read_weights,
)

logger = logging.getLogger(os.path.basename(__file__))


def set_antimeridian(dataarray, to: str):
    """Flip the antimeridian (i.e. longitude discontinuity) between Europe
    (i.e., [0, 360)) and the Pacific (i.e., [-180, 180)).

    Parameters
    ----------
    - dataarray : xarray.DataArray
    - to : string, {'pacific', 'europe'}
      * 'europe': Longitude will be in [0, 360)
      * 'pacific': Longitude will be in [-180, 180)

    Returns
    -------
    dataarray : xarray.DataArray
    """
    lon = dataarray['lon']

    if to.lower() == 'europe':
        dataarray = dataarray.assign_coords(lon=lon % 360)
    elif to.lower() == 'pacific':
        dataarray = dataarray.assign_coords(lon=((lon + 180) % 360) - 180)
    else:
        errmsg = "to has to be one of ['europe', 'pacific'] not {}".format(to)
        raise ValueError(errmsg)

    idx = np.argmin(dataarray['lon'].values)
    dataarray = dataarray.roll(lon=-idx, roll_coords=True)
    dataarray['lon'].attrs = lon.attrs
    return dataarray


def mapplot(dataarray, cfg, title_pattern, filename_part, ancestors,
            **colormesh_args):
    """Visualize weighted temperature."""
    metadata = read_metadata(cfg)
    metadata_future = metadata['tas_CLIM_future']
    start_year = metadata_future[0]['start_year']
    end_year = metadata_future[0]['end_year']

    period = f'{start_year}-{end_year}'
    if 'tas_CLIM_reference' in metadata:
        metadata_reference = metadata['tas_CLIM_reference']
        start_year_ref = metadata_reference[0]['start_year']
        end_year_ref = metadata_reference[0]['end_year']
        period = f'change: {period} minus {start_year_ref}-{end_year_ref}'
    metric = cfg['model_aggregation']
    if isinstance(metric, int):
        metric = f'{metric}perc'
    proj = ccrs.PlateCarree(central_longitude=0)
    figure, axes = plt.subplots(subplot_kw={'projection': proj})

    dataarray = set_antimeridian(dataarray, cfg.get('antimeridian', 'pacific'))
    dataarray = dataarray.dropna('lon', how='all').dropna('lat', how='all')

    dataarray.plot.pcolormesh(
        ax=axes,
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

    axes.coastlines()
    axes.set_xticks(cfg.get('xticks', default_xticks), crs=proj)
    axes.set_yticks(cfg.get('yticks', default_yticks), crs=proj)
    axes.xaxis.set_ticks_position('both')
    axes.yaxis.set_ticks_position('both')
    axes.xaxis.set_major_formatter(longitude_formatter)
    axes.yaxis.set_major_formatter(latitude_formatter)
    axes.set_xlabel('')
    axes.set_ylabel('')

    title = title_pattern.format(metric=metric, period=period)
    axes.set_title(title)

    filename_plot = get_plot_filename(filename_part, cfg)
    figure.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(figure)

    filename_data = get_diagnostic_filename(filename_part, cfg, extension='nc')
    dataarray.to_netcdf(filename_data)

    log_provenance(title, filename_plot, cfg, ancestors)
    log_provenance(title, filename_data, cfg, ancestors)


def visualize_and_save_temperature(temperature, cfg: dict, ancestors: list):
    """Wrap mapplot: absolute temperature."""
    title_pattern = ('Weighted {metric} temperature\n'
                     r'{period} ($\degree$C)')
    filename_part = 'temperature_change_weighted_map'
    mapplot(temperature,
            cfg,
            title_pattern,
            filename_part,
            ancestors,
            cmap='Reds')


def visualize_and_save_difference(temperature_difference, cfg: dict,
                                  ancestors: list):
    """Wrap mapplot: temperature difference between weighted and unweighted."""
    title_pattern = '\n'.join([
        'Weighted minus unweighted {metric} temperature',
        r'{period} ($\degree$C)',
    ])
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
        return calculate_percentiles(dataset, [metric],
                                     weights).squeeze('percentile', drop=True)
    if metric.lower() == 'mean':
        if weights is not None:
            dataset = dataset.weighted(weights)
        return dataset.mean('model_ensemble')

    if metric.lower() == 'median':
        return calculate_percentiles(dataset, [50],
                                     weights).squeeze('percentile', drop=True)

    errmsg = f'model_aggregation {metric} is not implemented!'
    raise NotImplementedError(errmsg)


def main(cfg):
    """Plot weighted temperature graph."""
    input_files = cfg['input_files']
    filename = cfg['weights']
    weights_path = Path(input_files[0]) / filename
    weights = read_weights(weights_path)

    metadata = read_metadata(cfg)
    models = metadata['tas_CLIM_future']
    model_data, model_data_files = read_model_data(models)

    # if a historical period is given calculate the change
    if 'tas_CLIM_reference' in metadata:
        models_hist = metadata['tas_CLIM_reference']
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

    visualize_and_save_difference(
        weighted_mean - unweighted_mean,
        cfg,
        model_data_files,
    )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
