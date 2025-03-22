"""Python example diagnostic."""
import logging
import os
from pathlib import Path


from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
)
from spy4cast import spy4cast
from spy4cast.dataset import Dataset
from esmvaltool.diag_scripts.shared import ProvenanceLogger


logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
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
        # Por si cogemos uno solo, sino nos vamos al alias que tiene toda la info. El Slice lo haces en la recipe, en la receta metes el predictando y el predictor
        predictand = variables['predictand'][0]
        z = Dataset(name=predictand['filename'], folder=cfg['work_dir']).open(var=predictand['short_name'])
        z_ppcessed = spy4cast.Preprocess(z, period=7, order=4)
        z_ppcessed.save(f'preprocessed_z_{predictand["alias"]}', folder=cfg['work_dir'])
        single_predictand = True 

    groups = group_metadata(input_data, 'alias', sort='alias')
    with ProvenanceLogger(cfg) as provenance_logger:
        
        for dataset in groups:
            logger.info("Processing dataset %s", dataset)
            variables = group_metadata(groups[dataset], 'variable_group')
            if single_predictand:
                if 'predictor' not in variables:
                    logger.info(f"Dataset {dataset} does not have predictor: {variables.keys()}")
                    continue
            else:
                if 'predictor' not in variables or 'predictand' not in variables:
                    logger.info(f"Dataset {dataset} incomplete: variables: {variables.keys()}")
                    continue

            predictor = variables['predictor'][0]
            y = Dataset(name=predictor['filename'], folder=cfg['work_dir']).open(var=predictor['short_name'])
            # y._region = Region(-30, 30, 80, -120, Month.JUN, Month.AUG, 1970, 2000)
            y_ppcessed = spy4cast.Preprocess(y, period=7, order=4)
            y_ppcessed.save('save_preprocessed_y_', folder=cfg['work_dir'])

            if not single_predictand:
                predictand = variables['predictand'][0]
                z = Dataset(name=predictand['filename'], folder=cfg['work_dir']).open(var=predictand['short_name'])
                z_ppcessed = spy4cast.Preprocess(z, period=7, order=4)
                z_ppcessed.save(f'preprocessed_z_{predictand["alias"]}', folder=cfg['work_dir'])

            nm = cfg['modes']
            alpha = cfg['alpha']
            mca = spy4cast.MCA(y_ppcessed, z_ppcessed, nm, alpha, sig='monte-carlo', montecarlo_iterations=10)
            logger.info("Saving MCA...")
            mca.save(f'mca_{dataset}', folder=cfg['work_dir'])
            logger.info("Plotting...")
            mca.plot(save_fig=True, cmap='viridis', folder=cfg['plot_dir'], name=f'mca-{dataset}.png',signs=[True,True,True],width_ratios=[0.5,1,1])
            provenance_logger.log(
                os.path.join(cfg['plot_dir'], f'mca-{dataset}.png'), 
                get_provenance_record(
                    caption=f'MCA analysis for {dataset}', 
                    ancestor_files=[predictand['filename'], predictor['filename']]
                )
            )
            logger.info('Finished!!!')

            if cfg['cross_validation']:
                cross = spy4cast.Crossvalidation(y_ppcessed, z_ppcessed, nm, alpha)
                cross.save(f'cross_{dataset}', folder=cfg['work_dir'])
                cross.plot(save_fig=True, folder=cfg['plot_dir'], name=f'crossvalidation-{dataset}.png')
                provenance_logger.log(
                os.path.join(cfg['plot_dir'], f'crossvalidation-{dataset}.png'), 
                get_provenance_record(
                    caption=f'Cross validated MCA analysis for {dataset}', 
                    ancestor_files=[predictand['filename'], predictor['filename']]
                )
            )


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
