"""Align the target model with the CMIP ensemble.

1. Get the global mean temperature change for specified years and
specified percentiles (CMIP Delta T). These define the scenarios.

- Select 30-year periods from the target model (all ensemble members)
where they match the CMIP Delta T for each scenario.

- Make a plot of the global mean temperature change
according to all datasets, and indicate the scenario's.

"""
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
        "Resampling of local climate model.",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            'alidoost_sarah',
            # 'rol_evert',
        ],
        'ancestors': ancestor_files,
    }
    return record


def mean_of_target_models(metadata, model):
    """Get the average delta T of the target model ensemble members."""
    target_model_data = select_metadata(metadata, variable_group='tas_target')
    files = [tmd['filename'] for tmd in target_model_data
             if not 'MultiModel' in tmd['dataset']]
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
    target_dts = target_dts.rolling(time=30, center=True, min_periods=30).mean()
    time_idx = abs(target_dts - cmip_dt).argmin(dim='time').values
    year = target_dts.isel(time=time_idx).year.values.astype(int)
    target_dt = target_dts.isel(time=time_idx).values.astype(float)
    return [year - 14, year + 15], target_dt


def _timeline(ax, yloc, interval, text):
    """Plot an interval near the bottom right of the plot."""
    xmin, xmax = interval
    xcenter = (xmin + xmax) / 2
    xspan = (xmax - xmin) / 2  # use 'error bar: extends to both sides.

    # Later years should be located slightly higher:
    yloc = 0.05 + yloc/20
    ax.errorbar(xcenter, yloc, xerr=xspan, capsize=5, lw=2, capthick=2,
                transform=ax.get_xaxis_transform(), c='r', label='Resampling periods')


def make_plot(metadata, scenarios, cfg):
    """Make figure 3, left graph.

    Multimodel values as line, reference value in black square,
    steering variables in dark dots.
    """
    fig, ax = plt.subplots()
    for member in select_metadata(metadata, variable_group='tas_cmip'):
        filename = member['filename']
        dataset = xr.open_dataset(filename)
        if not 'MultiModel' in filename:
            ax.plot(dataset.time.dt.year, dataset.tas.values,
                    c='grey', alpha=0.3, lw=.5, label='CMIP members')
        else:
            statistic = 'CMIP ' + Path(filename).stem.split('_')[0][10:]
            # e.g. get "CMIP Mean" from "path_to/file/MultiModelMean_Amon_tas....nc"

            # Only display stats for the future period:
            dataset = dataset.sel(time=slice('2010', None, None))
            ax.plot(dataset.time.dt.year, dataset.tas.values,
                    c='k', lw=2, label=statistic)

    for member in select_metadata(metadata, variable_group='tas_target'):
        filename = member['filename']
        dataset = xr.open_dataset(filename)
        if not 'MultiModel' in filename:
            ax.plot(dataset.time.dt.year, dataset.tas.values,
                    c='blue', lw=1, label=member['dataset'])

    # Add the scenario's with dots at the cmip dt and bars for the periods
    dotlabel = r"Scenarios' steering $\Delta T_{CMIP}$"
    for i, scenario in enumerate(scenarios):
        ax.scatter(scenario['year'], scenario['cmip_dt'], s=50,
                         zorder=10, color='r', label=dotlabel)
        _timeline(ax, i, scenario['period_bounds'], scenario)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # dict removes dupes
    ax.legend(by_label.values(), by_label.keys())
    ax.set_xlabel('Year')
    ax.set_ylabel(r'Global mean $\Delta T$ (K) w.r.t. reference period')

    # Save figure
    filename = get_plot_filename('global_matching', cfg)
    fig.savefig(filename, bbox_inches='tight', dpi=300)


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
    target_dts, provenance = mean_of_target_models(
        metadata, cfg['target_model']
    )

    # Define the different scenario's
    scenarios = []
    combinations = product(cfg['scenario_years'], cfg['scenario_percentiles'])
    for year, percentile in combinations:
        cmip_dt = get_cmip_dt(metadata, year, percentile)
        bounds, target_dt = get_resampling_period(
            target_dts, cmip_dt
        )

        scenario = {
            'year': year,
            'percentile': percentile,
            'cmip_dt': cmip_dt,
            'period_bounds': bounds,
            'target_dt': float(target_dt),
            'pattern_scaling_factor': cmip_dt / target_dt
        }
        scenarios.append(scenario)

    save(scenarios, cfg, provenance)

    if cfg['write_plots']:
        make_plot(metadata, scenarios, cfg)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
