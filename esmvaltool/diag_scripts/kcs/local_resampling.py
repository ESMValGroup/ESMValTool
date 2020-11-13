"""Resample the target model for the selected time periods."""
import logging
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, run_diagnostic,
                                            select_metadata)

LOGGER = logging.getLogger(Path(__file__).name)


def _create_provenance_record(ancestor_files):
    """Create a provenance record."""
    record = {
        'caption': "Resampling of local climate model.",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            'alidoost_sarah',
            'rol_evert',
        ],
        'ancestors': ancestor_files,
    }
    return record


def _get_data_target_model(cfg):
    """Load ensembles for each variable and merge."""
    LOGGER.info("Reading input data for target model")
    dataset_dicts = cfg['input_data'].values()
    dataset = []
    ancestor_files = []
    for short_name in "pr", "tas":
        group = f'{short_name}_target'
        var = select_metadata(dataset_dicts, variable_group=group)
        files = [metadata['filename'] for metadata in var]
        dataset.append(
            xr.open_mfdataset(files,
                              concat_dim='ensemble_member',
                              combine='nested'))
        ancestor_files.extend(files)
    provenance = _create_provenance_record(ancestor_files)
    return xr.merge(dataset).load(), provenance


def _segment(dataset, period, step=5):
    """Compute season means for each n-year segment of input dataset."""
    segments = []
    for year in range(*period, step):
        segments.append(
            dataset.sel(time=slice(str(year), str(year + step - 1))))
    segmented_dataset = xr.concat(segments, dim='segment')
    return segmented_dataset


def get_segment_season_means(cfg):
    """Return a dict with segment season means for control and all futures.

    Read input data for the target model; extract segmented subsets for both
    the control and future periods; precompute seasonal means for each segment.

    Store intermediate results; if files already exist, load them instead.
    """
    dataset, provenance = _get_data_target_model(cfg)

    # Combine the future scenarios (name and resampling period)
    # with the control period in a single dictionary
    periods = {
        name: info['resampling_period']
        for name, info in cfg['scenarios'].items()
    }
    periods['control'] = cfg['control_period']

    # Get the segment season means for each of these periods
    segments_season_means = {}
    for name, period in periods.items():
        LOGGER.info("Get segment season means for %s", name)
        filename = f"{cfg['run_dir']}/season_means_{name}.nc"
        if Path(filename).exists():
            LOGGER.info("Found intermediate file %s", filename)
        else:
            LOGGER.info("Computing seasonal means for segmented dataset")
            means = _segment(dataset, period).groupby('time.season').mean()
            means.to_netcdf(filename)
            LOGGER.info("Intermediate results stored as %s", filename)
        segments_season_means[name] = xr.open_dataset(filename)
    return segments_season_means, provenance


def _find_single_top1000(segment_means, target):
    """Select n_combinations that are closest to the target.

    First, get all possible combinations, then select n_combinations.
    Store the result under 'filename'
    """
    n_segments = len(segment_means.segment)
    n_members = len(segment_means.ensemble_member)
    segment_indices = range(n_segments)
    segment_means = segment_means.values  # much faster indexing

    all_possible_combinations = product(range(n_members), repeat=n_segments)
    results = []

    for combination in all_possible_combinations:
        results.append(
            list(combination) +
            [abs(segment_means[segment_indices, combination].mean() - target)])

    # Create a pandas dataframe with the combinations and distance to target
    dataframe = pd.DataFrame(results,
                             columns=list(segment_indices) + ['distance'])
    top1000 = dataframe.sort_values('distance').head(1000)
    return top1000


def get_all_top1000s(cfg, segment_season_means):
    """Return a dict with 1000 combinations for control and all futures.

    For control, these samples should have the same mean winter
    precipitation as the overall mean of the x ensemble members.

    For the future periods, the target value is a relative change
    with respect to the overall mean of the control period.
    """
    # Create a dict of target values for control and all futures
    control_segments = segment_season_means['control'].pr.sel(season='DJF')
    control_mean = control_segments.mean().values
    target_values = {'control': control_mean}
    for name, info in cfg['scenarios'].items():
        target_values[name] = control_mean * (1 + info['dpr_winter'] / 100)

    # Find the 1000 recombinations that are closest to the target values
    top1000s = {}
    for name, target in target_values.items():
        LOGGER.info('Get 1000 recombinations for %s', name)
        filename = f"{cfg['run_dir']}/top1000_{name}.csv"
        if Path(filename).exists():
            LOGGER.info("Found intermediate file %s", filename)
        else:
            segments = segment_season_means[name].pr.sel(season='DJF')
            top1000 = _find_single_top1000(segments, target)
            top1000.to_csv(filename, index=False)
            LOGGER.info("Intermediate results stored as %s.", filename)
        top1000s[name] = pd.read_csv(filename)
    return top1000s


def _index_with_xarray(combinations):
    """Generate indexer for all selected combinations of segmented dataset.

    combinations: numpy 2d array with shape (n_combinations, n_segments)

    Note that numpy indexing is much faster than xarray labelled indexing,
    but in this case it's nice to keep working with xarray labelled arrays.
    """
    n_segments = len(combinations[0])

    # Create a DataArray once...
    indices = xr.DataArray(data=np.zeros(n_segments),
                           dims=['segment'],
                           coords={'segment': np.arange(n_segments)})

    # ...and update its values for each selected combination
    for combination in combinations:
        indices.values = combination
        yield indices


def _season_means(combinations, segment_means):
    """Compute summer pr,and summer and winter tas for recombined climates."""
    interesting_variables = []
    columns = ['combination', 'pr_summer', 'tas_winter', 'tas_summer']
    for combination in _index_with_xarray(combinations):
        recombined_segments = segment_means.sel(ensemble_member=combination)
        season_means = recombined_segments.mean('segment')

        interesting_variables.append([
            combination.values,
            season_means.pr.sel(season='JJA').values,
            season_means.tas.sel(season='DJF').values,
            season_means.tas.sel(season='JJA').values
        ])
    return pd.DataFrame(interesting_variables, columns=columns)


def _within_bounds(values, bounds):
    """Return true if value is within percentile bounds."""
    low, high = np.percentile(values, bounds)
    return values.between(low, high)


def _get_subset(top1000, info, period):
    """Select samples based on the percentile bounds.

    Select samples for which summer pr, and summer and winter tas
    are within the percentile bounds specified in the recipe.
    """
    pr_summer = _within_bounds(top1000['pr_summer'],
                               info[f'pr_summer_{period}'])
    tas_winter = _within_bounds(top1000['tas_winter'],
                                info[f'tas_winter_{period}'])
    tas_summer = _within_bounds(top1000['tas_summer'],
                                info[f'tas_summer_{period}'])
    subset = top1000[np.all([pr_summer, tas_winter, tas_summer], axis=0)]
    return subset


def get_percentile_subsets(cfg, segment_season_means, top1000s):
    """Get subsets based on percentile ranges.

    For each set of 1000 samples, compute summer precipitation,
    and summer and winter temperature.
    Then, for each scenario, select samples for which
    summer precipitation, and summer and winter temperature are
    within the percentile bounds specified in the recipe.
    """
    # Overwrite top1000s with seasonal mean characteristics
    for name, dataframe in top1000s.items():
        LOGGER.info(
            "Compute summer mean pr and summer and winter "
            "mean tas for 1000 selected combinations for %s", name)
        segment_means = segment_season_means[name]
        combinations = dataframe.drop('distance', axis=1).values
        top1000s[name] = _season_means(combinations, segment_means)

    # For each scenario, get a subset for the control and future period.
    subsets = {}
    for scenario, info in cfg['scenarios'].items():
        LOGGER.info("Get percentile-based subsets for scenario %s", scenario)
        subsets[scenario] = {
            'control': _get_subset(top1000s['control'], info, 'control'),
            'future': _get_subset(top1000s[scenario], info, 'future')
        }
    return subsets


def _penalties(overlap):
    """Determine penalties dependent on the number of overlaps."""
    return np.piecewise(
        overlap,
        condlist=[overlap < 3, overlap == 3, overlap == 4, overlap > 4],
        funclist=[0, 1, 5, 100])


def _best_subset(combinations, n_sample=8):
    """Find n samples with minimal reuse of ensemble members per segment.

    combinations: a pandas series with the remaining candidates
    n: the final number of samples drawn from the remaining set.
    """
    # Convert series of 1d arrays to 2d array (much faster!)
    combinations = np.array(
        [list(combination) for combination in combinations])

    # Store the indices in a nice dataframe
    n_segments = combinations.shape[1]
    best_subset = pd.DataFrame(
        data=None,
        columns=[f'Segment {x}' for x in range(n_segments)],
        index=[f'Combination {x}' for x in range(n_sample)])

    # Random number generator
    rng = np.random.default_rng()

    lowest_penalty = 500  # just a random high value
    for _ in range(10000):
        subset = rng.choice(combinations, size=n_sample)
        penalty = 0
        for segment in subset.T:
            _, counts = np.unique(segment, return_counts=True)
            penalty += _penalties(counts).sum()
        if penalty < lowest_penalty:
            lowest_penalty = penalty
            best_subset.loc[:, :] = subset

    return best_subset


def select_final_subset(cfg, subsets, prov=None):
    """Select sample with minimal reuse of ensemble segments.

    Final set of eight samples should have with minimal reuse
    of the same ensemble member for the same period.
    From 10.000 randomly selected sets of 8 samples, count
    and penalize re-used segments (1 for 3*reuse, 5 for 4*reuse).
    Choose the set with the lowest penalty.
    """
    n_samples = cfg['n_samples']
    all_scenarios = {}
    for scenario, dataframes in subsets.items():
        # Make a table with the final indices
        LOGGER.info("Selecting %s final samples for scenario %s", n_samples,
                    scenario)
        control = _best_subset(dataframes['control'].combination, n_samples)
        future = _best_subset(dataframes['future'].combination, n_samples)
        table = pd.concat([control, future],
                          axis=1,
                          keys=['control', 'future'])
        all_scenarios[scenario] = table

        # Store the output
        filename = get_diagnostic_filename(f'indices_{scenario}',
                                           cfg,
                                           extension='csv')
        table.to_csv(filename)
        LOGGER.info("Selected recombinations for scenario %s: \n %s", scenario,
                    table)
        LOGGER.info('Output stored as %s', filename)

        # Write provenance information
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(filename, prov)
    return all_scenarios


def _cmip_envelope(datasetlist, variable, target_year):
    """Determine the change in <variable> PDF of each CMIP model.

    Note: using mf_dataset not possible due to different calendars.
    """
    cmip = select_metadata(datasetlist, variable_group=f'{variable}_cmip')
    envelope = []
    ancestors = []
    for data_dict in cmip:
        dataset = xr.open_dataset(data_dict['filename'])[variable]
        control = dataset.sel(time=slice('1981', '2010'))
        future = dataset.sel(time=slice(str(target_year -
                                            15), str(target_year + 15)))

        quantiles = [.05, .1, .25, .5, .75, .90, .95]
        qcontrol = control.groupby('time.season').quantile(quantiles)
        qfuture = future.groupby('time.season').quantile(quantiles)

        if variable == 'tas':
            # absolute diff
            envelope.append(qfuture - qcontrol)
        else:
            # pr; relative diff
            envelope.append((qfuture - qcontrol) / qcontrol * 100)
        ancestors.append(data_dict['filename'])

    cmip = xr.concat(envelope, dim='multimodel')
    provenance = _create_provenance_record(ancestors)

    # Prevent confusion between dimension 'quantile' and method 'quantile'
    return cmip.rename({'quantile': 'percentile'}), provenance


def _recombine(segments, combinations):
    """Recombine segments according to the final combinations."""
    n_segments = len(segments.segment)
    new_climates = []
    for _, indices in combinations.iterrows():
        # Create indexer array
        indexer = xr.DataArray(indices,
                               dims=['segment'],
                               coords={'segment': range(n_segments)})

        # Recombine the segments using the indexer
        resample = segments.sel(ensemble_member=indexer).mean('segment')
        new_climates.append(resample)
    return xr.concat(new_climates, dim='sample')


def _get_climatology(cfg, scenario_name, table):
    """Determine the change in <variable> PDF of each scenario."""
    dataset, _ = _get_data_target_model(cfg)

    future = cfg['scenarios'][scenario_name]['resampling_period']
    segments_control = _segment(dataset, cfg['control_period'])
    segments_future = _segment(dataset, future)

    resampled_control = _recombine(segments_control, table['control'])
    resampled_future = _recombine(segments_future, table['future'])

    quantiles = [.05, .1, .25, .5, .75, .90, .95]
    qcontrol = resampled_control.groupby('time.season').quantile(
        quantiles, dim=['sample', 'time'])
    qfuture = resampled_future.groupby('time.season').quantile(
        quantiles, dim=['sample', 'time'])

    qchange_tas = (qfuture - qcontrol).tas
    qchange_pr = ((qfuture - qcontrol) / qcontrol * 100).pr
    return xr.merge([qchange_tas, qchange_pr])


def make_plots(cfg, scenario_tables):
    """Reproduce figure 5 from the paper."""
    # Note that quantile is applied twice! Once to get the pdf's of seasonal
    # tas/pr and once to get the multimodel pdf of the quantile changes
    metadata = cfg['input_data'].values()

    climates = {}
    for name, info in cfg['scenarios'].items():
        climatology = _get_climatology(cfg, name, table=scenario_tables[name])
        climates[name] = climatology

    for year in [2050, 2085]:
        fig, subplots = plt.subplots(2, 2, figsize=(12, 8))

        for row, variable in zip(subplots, ['pr', 'tas']):
            cmip, prov = _cmip_envelope(metadata, variable, year)

            for axes, season in zip(row, ['DJF', 'JJA']):
                percentiles = cmip.percentile.values
                xlocs = np.arange(len(percentiles))

                # Plot the cmip envelope
                seasondata = cmip.sel(season=season)
                for high, low in [[0.9, 0.1], [0.75, 0.25]]:
                    upper = seasondata.quantile(high, dim='multimodel')
                    lower = seasondata.quantile(low, dim='multimodel')
                    axes.fill_between(xlocs, upper, lower, color='k', alpha=.3)
                    axes.set_title(f'{variable} / {season}')

                # Plot the recombined scenarios
                for name, info in cfg['scenarios'].items():
                    if year == info['scenario_year']:
                        climate = climates[name].sel(season=season)[variable]
                        axes.plot(xlocs, climate, lw=3, label=name)

                axes.set_xticks(xlocs)
                axes.set_xticklabels([f'P{100*x:02.0f}' for x in percentiles])
        subplots[0, 0].set_ylabel('change (%)')
        subplots[1, 0].set_ylabel('change (K)')
        subplots[1, 1].legend()
        filename = get_plot_filename(f'local_validation_{year}', cfg)
        fig.suptitle(f'Year: {year}')
        fig.savefig(filename, bbox_inches='tight', dpi=300)
        LOGGER.info("Envelope figure stored as %s", filename)
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(filename, prov)


def main(cfg):
    """Resample the model of interest."""
    # Step 0: extract segmented subsets and precompute seasonal means
    segment_season_means, provenance = get_segment_season_means(cfg)

    # Step 1: get 1000 combinations
    top1000s = get_all_top1000s(cfg, segment_season_means)

    # Step 2: select samples based on the the percentile bounds
    subsets = get_percentile_subsets(cfg, segment_season_means, top1000s)

    # Step 3: select final set of eight samples
    scenarios = select_final_subset(cfg, subsets, prov=provenance)

    # Step 4: plot the results
    if cfg['write_plots']:
        make_plots(cfg, scenarios)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
