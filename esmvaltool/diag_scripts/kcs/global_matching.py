"""Align the target model with the CMIP ensemble."""
import logging
from pathlib import Path
from itertools import product

import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

from esmvaltool.diag_scripts.shared import (get_diagnostic_filename,
                                            get_plot_filename, run_diagnostic,
                                            select_metadata, ProvenanceLogger)

logger = logging.getLogger(Path(__file__).name)


def create_provenance_record(ancestor_files):
    """Create a provenance record."""
    record = {
        'caption':
        "Match temperature anomaly in target model to CMIP ensemble",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            'alidoost_sarah',
            'rol_evert',
        ],
        'ancestors': ancestor_files,
    }
    return record


def mean_of_target_models(metadata):
    """Get the average delta T of the target model ensemble members."""
    target_model_data = select_metadata(metadata, variable_group='tas_target')
    files = [
        tmd['filename'] for tmd in target_model_data
        if 'MultiModel' not in tmd['dataset']
    ]
    datasets = xr.open_mfdataset(files, combine='nested', concat_dim='ens')
    provenance = create_provenance_record(files)
    return datasets.tas.mean(dim='ens'), provenance


def get_cmip_dt(metadata, year, percentile):
    """Compute target delta T for KNMI scenarios."""
    attribute = f'MultiModel{percentile}'
    multimodelstat = select_metadata(metadata, alias=attribute)[0]
    dataset = xr.open_dataset(multimodelstat['filename'])
    return dataset.tas.sel(time=str(year)).values[0]


def get_resampling_period(target_dts, cmip_dt):
    """Return 30-year time bounds of the resampling period.

    This is the period for which the target model delta T
    matches the cmip delta T for a specific year.
    Uses a 30-year rolling window to get the best match.
    """
    target_dts = target_dts.rolling(time=30, center=True,
                                    min_periods=30).mean()
    time_idx = abs(target_dts - cmip_dt).argmin(dim='time').values
    year = target_dts.isel(time=time_idx).year.values.astype(int)
    target_dt = target_dts.isel(time=time_idx).values.astype(float)
    return [year - 14, year + 15], target_dt


def _timeline(axes, yloc, interval):
    """Plot an interval near the bottom of the plot."""
    xmin, xmax = interval

    # Later years should be located slightly higher:
    # yloc is relative to the axes, not in data coordinates.
    yloc = 0.05 + yloc / 20

    plot_args = dict(transform=axes.get_xaxis_transform(),
                     linewidth=2,
                     color='red')

    axes.plot([xmin, xmax], [yloc] * 2, **plot_args, label='Selected periods')
    axes.plot([xmin] * 2, [yloc - 0.01, yloc + 0.01], **plot_args)
    axes.plot([xmax] * 2, [yloc - 0.01, yloc + 0.01], **plot_args)


def make_plot(metadata, scenarios, cfg, provenance):
    """Make figure 3, left graph.

    Multimodel values as line, reference value in black square,
    steering variables in dark dots.
    """
    fig, axes = plt.subplots()
    for member in select_metadata(metadata, variable_group='tas_cmip'):
        filename = member['filename']
        dataset = xr.open_dataset(filename)
        if 'MultiModel' not in filename:
            axes.plot(dataset.time.dt.year,
                      dataset.tas.values,
                      c='grey',
                      alpha=0.3,
                      lw=.5,
                      label='CMIP members')
        else:
            # Only display stats for the future period:
            dataset = dataset.sel(time=slice('2010', None, None))
            axes.plot(dataset.time.dt.year,
                      dataset.tas.values,
                      color='k',
                      linewidth=2,
                      label='CMIP ' + Path(filename).stem.split('_')[0][10:])

    for member in select_metadata(metadata, variable_group='tas_target'):
        filename = member['filename']
        dataset = xr.open_dataset(filename)
        if 'MultiModel' not in filename:
            axes.plot(dataset.time.dt.year,
                      dataset.tas.values,
                      color='blue',
                      linewidth=1,
                      label=member['dataset'])

    # Add the scenario's with dots at the cmip dt and bars for the periods
    for i, scenario in enumerate(scenarios):
        axes.scatter(scenario['year'],
                     scenario['cmip_dt'],
                     s=50,
                     zorder=10,
                     color='r',
                     label=r"Scenarios' steering $\Delta T_{CMIP}$")
        _timeline(axes, i, scenario['period_bounds'])

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # dict removes dupes
    axes.legend(by_label.values(), by_label.keys())
    axes.set_xlabel('Year')
    axes.set_ylabel(r'Global mean $\Delta T$ (K) w.r.t. reference period')

    # Save figure
    filename = get_plot_filename('global_matching', cfg)
    fig.savefig(filename, bbox_inches='tight', dpi=300)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance)


def save(output, cfg, provenance):
    """Save the output as csv file."""
    scenarios = pd.DataFrame(output)
    filename = get_diagnostic_filename('scenarios', cfg, extension='csv')
    scenarios.to_csv(filename)
    print(scenarios.round(2))
    print(f"Output written to {filename}")
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance)


def main(cfg):
    """Return scenarios tables."""
    # A list of dictionaries describing all datasets passed on to the recipe
    metadata = cfg['input_data'].values()

    # Get the average delta T of the target model
    target_dts, provenance = mean_of_target_models(metadata)

    # Define the different scenario's
    scenarios = []
    combinations = product(cfg['scenario_years'], cfg['scenario_percentiles'])
    for year, percentile in combinations:
        cmip_dt = get_cmip_dt(metadata, year, percentile)
        bounds, target_dt = get_resampling_period(target_dts, cmip_dt)

        scenario = {
            'year': year,
            'percentile': percentile,
            'cmip_dt': cmip_dt,
            'period_bounds': bounds,
            'target_dt': float(target_dt),
            'pattern_scaling_factor': cmip_dt / target_dt
        }
        scenarios.append(scenario)

    # Save scenarios tables as csv file
    save(scenarios, cfg, provenance)

    # Plot the results
    if cfg['write_plots']:
        make_plot(metadata, scenarios, cfg, provenance)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
