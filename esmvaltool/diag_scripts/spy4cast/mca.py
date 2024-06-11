"""Python example diagnostic."""
import logging
import os
from pathlib import Path
from pprint import pformat


from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    select_metadata,
    sorted_metadata,
)
from spy4cast import spy4cast, Month, Region
from spy4cast.dataset import Dataset
from spy4cast.meteo import Clim

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Average {long_name} between {start_year} and {end_year} "
               "according to {dataset}.".format(**attributes))

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': ['zonal'],
        'authors': [
            'vegas-regidor_javier',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record



def main(cfg):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    variables = group_metadata(input_data, 'variable_group', sort='dataset')
    single_predictand = False
    if len(variables['predictand']) == 1:
        predictand = variables['predictand'][0]
        z = Dataset(name=predictand['filename'], folder=cfg['work_dir']).open(var=predictand['short_name'])
        z_ppcessed = spy4cast.Preprocess(z, period=7, order=4)
        z_ppcessed.save(f'preprocessed_z_{predictand["alias"]}', folder=cfg['work_dir'])
        single_predictand = True #Por si cogemos uno solo, sino nos vamos al alias que tiene toda la info. El Slice lo haces en la recipe, en la receta metes el predictando y el predictor

    groups = group_metadata(input_data, 'alias', sort='alias')
    for dataset in groups:
        logger.info("Processing dataset %s", dataset)
        variables = group_metadata(groups[dataset], 'variable_group')
        if single_predictand:
            if 'predictor' not in variables:
                logger.info(f"Dataset {dataset} does not have predictand: {variables.keys()}")
                continue
        else:
            if 'predictor' not in variables or 'predictand' not in variables:
                logger.info(f"Dataset {dataset} incomplete: variables: {variables.keys()}")
                continue

        predictor = variables['predictor'][0]
        y = Dataset(name=predictor['filename'], folder=cfg['work_dir']).open(var=predictor['short_name'])
        # y.q(Region(-5, 10, -30, -10, Month.JUN, Month.JUL, 1850, 2000), skip=0)
        y_ppcessed = spy4cast.Preprocess(y, period=7, order=4)
        y_ppcessed.save('save_preprocessed_y_', folder=cfg['work_dir'])

        if not single_predictand:
            predictand = variables['predictand'][0]
            z = Dataset(name=predictand['filename'], folder=cfg['work_dir']).open(var=predictand['short_name'])
            # z.slice(Region(35, 45, 345, 360, Month.JUN, Month.JUL, 1850, 2000), skip=0)
            z_ppcessed = spy4cast.Preprocess(z, period=7, order=4)
            z_ppcessed.save(f'preprocessed_z_{predictand["alias"]}', folder=cfg['work_dir'])

        nm = cfg['modes']
        alpha = cfg['alpha']
        mca = spy4cast.MCA(y_ppcessed, z_ppcessed, nm, alpha, sig='monte-carlo', montecarlo_iterations=10)
        #mca = spy4cast.MCA(y_ppcessed, z_ppcessed, nm, alpha)
        logger.info("Saving MCA...")
        mca.save(f'mca_{dataset}', folder=cfg['work_dir']) #Aquí es donde lo guarda
        logger.info("Plotting...")
        mca.plot(save_fig=True, cmap='viridis', folder=cfg['plot_dir'], name=f'mca-{dataset}.png')
        logger.info('Finished!!!')

        if cfg['cross_validation']: #Meter cross si podemos
            cross = spy4cast.Crossvalidation(y_ppcessed, z_ppcessed, nm, alpha)
            cross.save(f'cross_{dataset}', folder=cfg['work_dir'])
            cross.plot(save_fig=True, folder=cfg['plot_dir'], name=f'crossvalidation-{dataset}.png')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
