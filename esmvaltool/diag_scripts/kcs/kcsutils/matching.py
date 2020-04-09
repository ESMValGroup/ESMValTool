"""
Match experiment by ensemble ID
Matches historical and future experiments
match_experiments:

What attribute(s) to match the historical and future experiments by
Valid values: 'model', 'ensemble'
The default is by 'ensemble'
match_by: 'ensemble'

What to do with future scenarions that don't have a
matching historical run?
Valid values: 'remove', 'random', 'randomrun', 'exception'
- 'remove': if historical data is missing, remove the
future experiment.
- 'random': pick an entirely random historical run from the model
- 'randomrun': keep i & p the same, find a random
  historical run with a (different) r. If there's still no
  match to be found, an exception is raised.
- 'error' the script raises an exception if no matching
  historical run is found.
The default is 'error': to raise an exception
fix_nonmatching_historical: randomrun
"""

import logging


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


def match(dataset, match_by=None, on_no_match=None, historical_key=None):
    """Match dataset experiments, historical with future experiments

    Anything where the experiment does not match the `historical_key`
    is assumed to be a future experiment (rcp, ssp).

    Parameters
    ----------

    - match_by: 'ensemble' or 'model'

    - on_no_match: 'remove', 'random', 'randomrun', or 'error'

          What to do with a future run that has no matching historical run.

          - 'remove' the dataset

          - 'random': pick an entirely random historical run from the model

          - 'randomrun': keep i # & p the same, find a random
            historical run with a (different) r. If there's still no
            match to be found, an exception is raised.

          - 'error' the script raises an exception if no matching
            historical run is found.

    - historical_key: string

        What value for the 'experiment' attribute indicates a historical run?

    - attributes: dict, or None

        The attribute names for the following info. This information
        is given as a dict, with the keys below. The default values
        are given after each key, and is used when attributes=None.
        - experiment: "experiment"
        - model: "model_id"
        - realization: "realization"
        - initialization: "initialization"
        - physics: "physics_version"

      Note that the dictionary does not have to be complete, depending
      on the choice of `match_by` and `on_no_match`.

    """

    if match_by is None:
        match_by = "ensemble"
    if on_no_match is None:
        on_no_match = "randomrun"
    if historical_key is None:
        historical_key = "history"

    dataset['index_match_run'] = -1

    for model, group in dataset.groupby('model'):
        future_sel = group['experiment'] != historical_key
        hist_sel = group['experiment'] == historical_key
        if not any(hist_sel):
            logger.warning("Model %s has no historical runs", model)
            continue

        for row in group.loc[future_sel, :].itertuples():
            index = row.Index
            experiment = row.experiment
            prip = row.prip
            ensemble = f"r{row.realization}i{row.initialization}p{row.physics}"
            if prip:
                realization, initialization, physics = prip
            else:
                logger.warning("parent RIP not available: matching by ensemble info directly")
                realization = row.realization
                initialization = row.initialization
                physics = row.physics
            if match_by == 'ensemble':
                sel = (hist_sel &
                       (group['realization'] == realization) &
                       (group['initialization'] == initialization) &
                       (group['physics'] == physics))
            else:  # matching by model: any historical run will do
                sel = hist_sel
            if not any(sel):
                if on_no_match == 'error':
                    msg = f"no historical data for {model}"
                    if match_by == 'ensemble':
                        msg += f", {ensemble}"
                    raise ValueError(msg)
                logger.warning("no matching historical run found for %s - %s - %s",
                               model, experiment, ensemble)
                # pylint: disable=consider-using-in
                if on_no_match == 'remove' or on_no_match == 'ignore':
                    continue
                if on_no_match == 'randomrun' and match_by == 'ensemble':
                    logger.info("Finding replacement historical run for %s - %s - %s",
                                model, experiment, ensemble)

                    # Grab from all of the dataset, select by index
                    sel = (hist_sel &
                           (group['initialization'] == initialization) &
                           (group['physics'] == physics))
                    if not any(sel):
                        raise ValueError("No matching random realization")
                elif on_no_match == 'random':
                    sel = hist_sel
                    logger.info("Using replacement historical run for %s - %s", model, experiment)
                else:
                    raise ValueError(f"Parameter 'on_no_match' has invalid value {on_no_match}")
            hrow = group.loc[sel, :].iloc[0, :]
            hindex = group.loc[sel, :].index.array[0]
            # Pick the first (if multiple, e.g. for match_by=='model') historical cubes
            if match_by == 'ensemble':
                logger.debug("Matching %s - %s - %s with %s - r%di%dp%d historical run",
                             model, experiment, ensemble, hrow['model'], hrow['realization'],
                             hrow['initialization'], hrow['physics'])
            else:
                logger.debug("Matching %s - %s with %s - r%di%dp%d historical run",
                             model, experiment, hrow['model'], hrow['realization'],
                             hrow['initialization'], hrow['physics'])

            dataset.loc[index, 'index_match_run'] = hindex

    # Select only the matched dataset
    # Future experiments that were not matched, are removed
    # Historical runs that were not matched, however, are kept
    sel = ((dataset['experiment'] == 'historical') |
           (dataset['index_match_run'] > -1))
    dataset = dataset.loc[sel, :]

    return dataset
