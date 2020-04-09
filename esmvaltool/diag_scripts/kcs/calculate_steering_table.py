"""
"""
import os
from pathlib import Path
import logging
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import iris
try:
    from iris.util import equalise_attributes
except ImportError:   # Iris 2
    from iris.experimental.equalise_cubes import equalise_attributes
import cftime
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvalcore.preprocessor import area_statistics
import kcsutils


REF_PERIOD = (1980, 2009)

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


def num2date(coord, index=None):
    """DUMMY DOCSTRING"""
    if index is None:
        return cftime.num2date(coord.points, str(coord.units), coord.units.calendar,
                               only_use_cftime_datetimes=False)
    return cftime.num2date(coord.points[index], str(coord.units),
                           coord.units.calendar, 
                           only_use_cftime_datetimes=False)


def extract_season(cubes, season):
    """DUMMY DOC-STRING"""
    constraint = iris.Constraint(season=kcsutils.EqualConstraint(season))
    logger.info("Extracting season %s", season)
    for cube in cubes:
        if not cube.coords('season'):
            iris.coord_categorisation.add_season(cube, 'time')
    cubes = list(map(constraint.extract, cubes))
    return cubes


def average_year_cube(cube, season=None):
    """DUMMY DOC-STRING"""
    if season:
        if not cube.coords('season_year'):
            iris.coord_categorisation.add_season_year(cube, 'time')
        mean = cube.aggregated_by('season_year', iris.analysis.MEAN)
    else:
        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')
        mean = cube.aggregated_by('year', iris.analysis.MEAN)
    return mean


def average_year(cubes, season=None):
    """DUMMY DOC-STRING"""
    logger.info("Calculating %s averages", season if season else 'yearly')
    averages = []
    for cube in cubes:
        if season:
            if not cube.coords('season_year'):
                iris.coord_categorisation.add_season_year(cube, 'time')
            average = cube.aggregated_by('season_year', iris.analysis.MEAN)
        else:
            if not cube.coords('year'):
                iris.coord_categorisation.add_year(cube, 'time')
            average = cube.aggregated_by('year', iris.analysis.MEAN)
        averages.append(average)
    return averages


def calc_reference_values(cubes, reference_period, normby='run'):
    """Calculate the reference values for each cube

    If normby is not 'run', the individual reference values are averaged, so that each run
    will later be normalized by the same value.

    """

    constraint = iris.Constraint(year=kcsutils.RangeConstraint(*reference_period))
    values = []
    for cube in cubes:
        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')
        mean = constraint.extract(cube)
        mean = mean.collapsed('time', iris.analysis.MEAN)
        values.append(mean.data)
    values = np.array(values)
    if normby != 'run':
        values = np.full(values.shape, values.mean(), dtype=values.dtype)
    return values


def normalize(cubes, refvalues, relative):
    """Normalize a single cube to its reference value

    If the value is a relative value, i.e., a percentual change, set
    the 'relative' parameter to `True`.

    """

    for cube, refvalue in zip(cubes, refvalues):
        cube.data -= refvalue
        if relative:
            cube.data /= refvalue
            cube.data *= 100
            cube.units = '%'
    return cubes


def calc_steering(dataset, distribution, scenarios, window_rolling_mean=0,
                  rounding=0, timespan=30, maxepoch=2100):
    """Parameters
    ----------
    - dataset: Iris.cube.Cube

        Input data. Historical and future experiment data should be
        concatenated into one cube.  Data should be yearly-averaged
        and normalized. If there are more realizations, these should
        be merged (averaged) into a single cube.

    - distribution: Pandas DataFrame

        CMIP percentile distribution. The index is years, the columns
        are 'mean' and percentiles of interest (e.g., 5, 10, 25, 50,
        75, 90, 95). Calculated with `kcs.tas_change`.

    - scenarios: list of dicts
      The dict should contain a name, year and percentile, e.g.
        scenarios=[{'name': 'G', 'epoch': 2050, 'percentile': 90},
                   {'name': 'L', 'epoch': 2050, 'percentile': 10},]
      The name should be unique.

    """

    if window_rolling_mean > 1:
        distribution = distribution.rolling(window_rolling_mean, center=True).mean().dropna()

    # Limit the dataset to 2085, so we don't try and calculate beyond 2100
    maxyear = distribution.index.max().year
    cube = iris.Constraint(year=lambda year: year <= maxyear).extract(dataset)

    for i, scenario in enumerate(scenarios):
        epoch = datetime(int(scenario['epoch']), 1, 1)
        percentile = scenario['percentile']
        delta_t = distribution.loc[epoch, str(percentile)]
        if rounding:  # nearest multiple of `round`
            rem = delta_t % rounding
            if rounding - rem > rem:
                # Round down
                delta_t -= rem
            else:
                # Round up
                delta_t += rounding - rem

        index = np.argmin(np.abs((cube.data - delta_t).data))
        date = num2date(cube[index].coord('time'))[0]
        # Fix below to use proper years
        period = ((date - timedelta(timespan/2*365.24)).year,
                  (date + timedelta(timespan/2*365.24)).year)
        if period[1] > maxepoch:
            period = maxepoch - timespan, maxepoch
            epoch = maxepoch - timespan//2
            delta = [d.year - epoch for d in num2date(cube.coord('time'))]
            index = np.argmin(np.abs(delta))
        scenarios[i]['cmip_delta_t'] = delta_t
        # Correct for the fact that our previous calculations were all on January 1.
        # We simply equate that to Dec 12 of the previous year, and thus make the
        # end-year of the period inclusive
        scenarios[i]['period'] = period[0], period[1]-1
        model_delta_t = cube[index].data
        scenarios[i]['model_delta_t'] = model_delta_t
        scenarios[i]['factor'] = delta_t / model_delta_t

    return scenarios


def normalize_average_dataset(cubes, season=None, average_years=True, relative=False,
                              reference_period=None):
    """Normalize and average a given iterable of cubes

    The dataset is normalized by, and averaged across, its individual
    ensemble runs. Thus, this works best (only), if the dataset
    belongs to the same model and experiment, and has no other
    ensemble variations other than its realization.

    The dataset should already be concatenated across the historical
    and future experiment, if necessary.

    Each Iris cube inside the dataset should have a scalar realization
    coordinate, with its value given the realization number. If not,
    these are added on the fly, equal to the iteration index of the
    cubes.

    """

    if reference_period is None:
        reference_period = REF_PERIOD

    if season:
        cubes = extract_season(cubes, season)
    if average_years:
        cubes = average_year(cubes, season=season)
    refvalues = calc_reference_values(
        cubes, reference_period=reference_period, normby='run')
    cubes = normalize(cubes, refvalues, relative=relative)

    for i, cube in enumerate(cubes):
        if not cube.coords('realization'):
            coord = iris.coords.AuxCoord(
                i, standard_name='realization', long_name='realization',
                var_name='realization')
            cube.add_aux_coord(coord)

    cubes = iris.cube.CubeList(cubes)
    equalise_attributes(cubes)
    cube2d = cubes.merge_cube()
    mean = cube2d.collapsed('realization', iris.analysis.MEAN)

    return mean


def calc(dataset, percentiles, scenarios, season=None, average_years=True,
         relative=False, reference_period=None,
         timespan=30, window_rolling_mean=0, rounding=None):
    """Calculate the percentile yearly change distribution for the input data

    Also performs extracting of season (optional), averaging of years
    (optional) and normalization to a common reference period (needed
    for a better inter-model comparison), before the percentiles are
    calculated.

    Returns
      2-tuple of

      - Percentiles, as Pandas DataFrame

      - Input dataset, but with possibly extracted seasons and
        averaged years, and normalized data

    """

    if reference_period is None:
        reference_period = REF_PERIOD

    mean = normalize_average_dataset(dataset['cube'], season, average_years,
                                     relative=relative, reference_period=reference_period)

    steering = calc_steering(mean, percentiles, scenarios, timespan=timespan,
                             window_rolling_mean=window_rolling_mean, rounding=rounding)
    return steering


def main(cfg):
    """

    """
    csvfile = Path(cfg['run_dir']).parent / 'script1' / cfg['input']
    percentiles = pd.read_csv(csvfile, index_col=0)
    percentiles.index = pd.to_datetime(percentiles.index)
    models = group_metadata(cfg['input_data'].values(), 'dataset')
    selection = models[cfg['model']]
    selection = group_metadata(selection, 'filename')
    paths = [Path(filename) for filename in selection.keys()]
    dataset = kcsutils.read_data(paths)
    dataset = kcsutils.match(
        dataset, match_by='ensemble', on_no_match='randomrun',
        historical_key='historical')
    # For historical and practical reasons, we concate the matched cubes
    # (Historical: the original code was developed where the model-of-interest
    #  was already concatenated, and in fact, didn't follow the CMIP/CF conventions)
    dataset = kcsutils.concat_cubes(dataset, historical_key='historical')
    columns = ['model', 'experiment', 'realization', 'initialization', 'physics', 'prip', 'var', 'index_match_run']
    print(dataset[columns])

    scenarios = [
        {'name': 'M', 'epoch': 2050, 'percentile': 10},
        {'name': 'W', 'epoch': 2050, 'percentile': 90},
        {'name': 'M', 'epoch': 2085, 'percentile': 10},
        {'name': 'W', 'epoch': 2085, 'percentile': 90},
    ]
    steering = calc(dataset, percentiles, scenarios=cfg['scenarios'],
                    reference_period=cfg['reference_period'],
                    timespan=cfg['timespan_years'],
                    window_rolling_mean=cfg['window_rolling_mean'])
    steering = pd.DataFrame(steering)

    steering.to_csv(cfg['output'], index=False)

    return


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
