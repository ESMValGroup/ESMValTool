"""A script to obtain the global mean temperature rise (change).

The global mean temperature rise is one of the two steering variables,
used in the interpolation of the regional climate model.

The minimum data for this script is one dataset with two different experiments:
historical and rcp (Representative Concentration Pathways).
For example:
- {dataset: ACCESS1-0, exp: historical, ensemble: r1i1p1, ...}
- {dataset: ACCESS1-0, exp: rcp45, ensemble: r1i1p1, ...}
"""
import os
from pathlib import Path
from datetime import datetime
import functools
import multiprocessing
import logging
import warnings
from pprint import pformat
import iris
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvalcore.preprocessor import area_statistics
import kcsutils


PERIOD = (1961, 2099)
MINDATA = {'historical': 20, 'future': 4}

logger = logging.getLogger(os.path.basename(__file__))


def extract_season(cubes, season):
    """DUMMY DOC-STRING"""
    # Use a class instead of a lambda function, so we can pass the
    # constraint to multiprocessing (which doesn't handle lambdas).
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
    if season:
        logger.info("Calculating %s averages", season)
    else:
        logger.info("Calculating yearly averages")
    func = functools.partial(average_year_cube, season=season)
    cubes = list(map(func, cubes))
    return cubes


class ModelReferencePointCalculation:
    """DUMMY DOC-STRING"""

    def __init__(self, dataset, reference_period,
                 historical_key="historical", yearly=True, season=None,
                 normby='run'):
        self.dataset = dataset
        self.historical_key = historical_key
        self.normby = normby
        if yearly:
            self.mindata = MINDATA
        elif season:  # in ['djf', 'mam', 'jja', 'son']:
            # Three months a year
            self.mindata = {key: 3*value for key, value in MINDATA.items()}
        else:
            # Twelve months a year
            self.mindata = {key: 12*value for key, value in MINDATA.items()}
        self.constraint = kcsutils.make_year_constraint_all_calendars(*reference_period)

    def __call__(self, model):
        """DUMMY DOC-STRING"""
        dataset = self.dataset[self.dataset['model'] == model]

        if self.normby == 'model':
            cubes = dataset['cube']
            histcubes = dataset['match_historical_run']

            value = self.calc_refvalue(cubes, histcubes, model)
        elif self.normby == 'experiment':
            value = {}
            for exp, group in dataset.groupby('experiment'):
                cubes = group['cube']
                histcubes = group['match_historical_run']
                value[exp] = self.calc_refvalue(cubes, histcubes, model)
        else:
            value = {}
            for index, row in dataset.iterrows():
                cubes = [row['cube']]
                histcubes = [row['match_historical_run']]
                value[index] = self.calc_refvalue(cubes, histcubes, model)

        return value

    def calc_refvalue(self, cubes, histcubes, model):
        """DUMMY DOC-STRING"""
        avs = {}
        avs['historical'] = self.calc_mean(histcubes, self.mindata['historical'], model)
        avs['future'] = self.calc_mean(cubes, self.mindata['future'], model)
        if not avs['historical'] or not avs['future']:  # Too few data to calculate a decent bias
            logger.warning("%s does not have enough data to compute a reference", model)
            return None
        logger.info("Calculating time-weighted reference value")
        ndata = {}
        mean = {}
        for key, values in avs.items():
            n = len(values)
            # Weighted time average for each section
            ndata[key] = sum(value[1] for value in values) / n
            mean[key] = sum(value[0].data for value in values) / n
        logger.debug("Reference data values: %s , with weights %s", pformat(mean), pformat(ndata))
        value = ((mean['historical'] * ndata['historical'] +
                  mean['future'] * ndata['future']) /
                 (ndata['historical'] + ndata['future']))

        return value

    def calc_mean(self, cubes, mindata, model):
        """DUMMY DOC-STRING"""
        averages = []
        for cube in cubes:
            calendar = cube.coord('time').units.calendar
            excube = self.constraint[calendar].extract(cube)
            if excube is None:
                logger.warning("A cube of %s does not support time range: %s",
                               model, cube.coord('time'))
                continue
            ndata = len(excube.coord('time').points)
            # Not enough of data? Ignore
            if ndata < mindata:
                logger.warning("A cube of %s has only %d data points for its time range",
                               model, ndata)
                continue
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", category=UserWarning,
                    message="Collapsing a non-contiguous coordinate. "
                    "Metadata may not be fully descriptive for 'year'")
                averages.append((excube.collapsed('time', iris.analysis.MEAN), ndata))
        return averages


def calc_reference_values(dataset, reference_period, yearly=False, season=None,
                          historical_key=None,
                          normby='run'):
    """Calculate reference values

    model:  *per model*, so that each realization is
    scaled to the reference period following the average *model* reference
    value

    """

    if not historical_key:
        historical_key = default_config['data']['historical_experiment']
    logger.info("Calculating reference values (period = %s)", reference_period)

    index = dataset['index_match_run']
    hindex = index[index > -1]
    cubes = dataset.loc[hindex, 'cube'].array
    future_data = dataset[index > -1].copy()
    future_data['match_historical_run'] = cubes

    calculation = ModelReferencePointCalculation(
        future_data, reference_period, yearly=yearly, season=season,
        historical_key=historical_key, normby=normby)

    models = dataset['model'].unique()
    reference_values = filter(None, map(calculation, models))

    if normby == 'model':
        ref_values = dataset['model'].map(dict(zip(models, reference_values)))
    elif normby == 'experiment':
        ref_values = []
        for model, values in zip(models, reference_values):
            sel = (dataset['model'] == model)
            experiments = dataset.loc[sel, 'experiment']
            ref_values.append(experiments.map(values))
        ref_values = pd.concat(ref_values)
    else:
        ref_values = pd.concat([pd.Series(values) for values in reference_values])

    return ref_values


def normalize_cube(item, relative=False):
    """Normalize a single cube to its reference value

    If the value is a relative value, i.e., a percentual change, set
    the 'relative' parameter to `True`.

    """

    cube, ref_value = item
    cube.data -= ref_value
    if relative:
        cube.data /= ref_value
        cube.data *= 100
        cube.units = '%'
    return cube


def normalize(dataset, relative=False, normby='run'):
    """Normalize cubes to their reference values

    If the value is a relative value, i.e., a percentual change, set
    the 'relative' parameter to `True`.

    This returns a new, modified dataset.

    """

    dataset['matched_exp'] = ''
    if normby != 'model':
        # Add the (double/triple/etc) historical runs for the future experiments,
        # and get rid of the old historical runs
        # We do this by
        # - selecting the indices & reference values for the matching runs
        # - copy the relevant data from the dataset
        #   We need to copy the Iris cubes, otherwise we'll have
        #   multiple references to the same cube
        # - Set the reference values and the matching future indices
        #   for the newly copied historical data
        # - Concatenate the current and new datasets, and drop all
        #   rows that don't have a reference value (original
        #   historical rows)

        sel = dataset['index_match_run'] > -1
        indices = dataset.loc[sel, 'index_match_run']

        reference_values = dataset.loc[sel, 'reference_value']
        hist_data = dataset.loc[indices, :]
        hist_data['cube'] = hist_data['cube'].apply(lambda cube: cube.copy())
        hist_data['reference_value'] = reference_values.array
        hist_data['index_match_run'] = dataset.loc[sel].index.array
        hist_data['matched_exp'] = dataset.loc[sel, 'experiment'].array
        dataset = pd.concat([dataset, hist_data])
        dataset.dropna(axis=0, subset=['reference_value'], inplace=True)
    logger.info("Normalizing data to the reference period")
    # Tempting, but it's not possible to simply do dataset['cube'] =
    # dataset['cube'] / data['reference_value']
    func = functools.partial(normalize_cube, relative=relative)
    data = zip(dataset['cube'], dataset['reference_value'])
    dataset['cube'] = list(map(func, data))
    return dataset


def calc_percentile_year(dataset, year, average_experiments=False):
    """Calculate the percentile distribution of the cubes for a given year"""

    constraint = iris.Constraint(year=kcsutils.EqualConstraint(year))

    if average_experiments:
        data = []
        for _, group in dataset.groupby(['model', 'experiment', 'matched_exp']):
            cubes = list(filter(None, map(constraint.extract, group['cube'])))
            if cubes:
                data.append(np.mean([cube.data for cube in cubes]))
    else:
        cubes = list(filter(None, map(constraint.extract, dataset['cube'])))
        data = [cube.data for cube in cubes]
    mean = np.mean(data)
    percs = np.percentile(data, [5, 10, 25, 50, 75, 90, 95], overwrite_input=True)
    return dict(zip(['mean', '5', '10', '25', '50', '75', '90', '95'],
                    [mean] + percs.tolist()))


def calc_percentiles(dataset, period=PERIOD, average_experiments=False):
    """DUMMY DOC-STRING"""

    logger.info("Calculating percentiles")

    years = list(range(*period))
    func = functools.partial(calc_percentile_year, dataset, average_experiments=average_experiments)
    percs = list(map(func, years))
    return pd.DataFrame(
        percs, index=pd.DatetimeIndex([datetime(year, 1, 1) for year in years]))


def calc(dataset, reference_period, historical_key=None, season=None, average_years=True,
         relative=False, period=PERIOD, normby='run', average_experiments=False):
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

    if not historical_key:
        historical_key = default_config['data']['attributes']['historical_experiment']

    if season:
        dataset['cube'] = extract_season(dataset['cube'], season)

    if average_years:
        dataset['cube'] = average_year(dataset['cube'], season=season)

    reference_values = calc_reference_values(
        dataset, reference_period, yearly=average_years, season=season,
        historical_key=historical_key,
        normby=normby)
    dataset['reference_value'] = reference_values
    dataset = normalize(dataset, relative=relative, normby=normby)

    percentiles = calc_percentiles(dataset, period=period,
                                   average_experiments=average_experiments)

    return percentiles, dataset


def main(cfg):
    """Return percentiles of temperature at each year.

    Four steps:
    1- obtain average values per each year,
    2- obtain the total average that is called reference_value,
    3- normalize temperature using refernce_value,
    4- calculate percentile for 5-95 perentage per each year.

    """
    # assemble the data dictionary keyed by dataset name
    # this makes use of the handy group_metadata function that
    # orders the data by 'dataset'; the resulting dictionary is
    # keyed on datasets e.g. dict = {'MPI-ESM-LR': [var1, var2...]}
    # where var1, var2 are dicts holding all needed information per variable
    logger.debug("\n\n\nCONFIG:\n %s\n\n\n", cfg)

    my_files_dict = group_metadata(cfg['input_data'].values(), 'filename')
    #pprint(my_files_dict)
    paths = [Path(filename) for filename in my_files_dict.keys()]
    dataset = kcsutils.read_data(paths)
    dataset = kcsutils.match(
        dataset, match_by='ensemble', on_no_match='randomrun',
        historical_key='historical')
    logger.debug("%s", dataset.columns)
    columns = ['model', 'experiment', 'realization', 'initialization', 'physics', 'prip', 'var', 'index_match_run']
    logger.debug("%s", dataset[columns])

    result, _ = calc(dataset,
                     reference_period=(1980, 2009),
                     historical_key="historical",
                     relative=False,
                     period=(1961, 2099), normby='run',
                     average_experiments=False)
    result.to_csv(cfg['output'], index_label="date")


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
