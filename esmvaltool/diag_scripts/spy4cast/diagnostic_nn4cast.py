"""Python example diagnostic."""
import logging
from pathlib import Path
import nn4cast.predefined_classes
import os
import tensorflow as tf



from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
)

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

    hyperparameters = cfg['hyperparameters'].copy()
    hyperparameters['outputs_path'] = cfg['work_dir']

    variables = group_metadata(input_data, 'variable_group', sort='dataset')
   
    single_predictand = False
    if len(variables['predictand']) == 1:
        # Por si cogemos uno solo, sino nos vamos al alias que tiene toda la info. El Slice lo haces en la recipe, en la receta metes el predictando y el predictor
        predictand = variables['predictand'][0]
        hyperparameters['path_y'] = predictand['filename']
        hyperparameters['name_y'] = predictand['short_name']
        predictand = variables['predictand'][0]
        single_predictand = True 

    groups = group_metadata(input_data, 'alias', sort='alias')
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
        hyperparameters['path_x'] = predictor['filename']
        hyperparameters['name_x'] = predictor['short_name']

        if not single_predictand:
            predictand = variables['predictand'][0]
            hyperparameters['path_y'] = predictand['filename']
            hyperparameters['name_y'] = predictand['short_name']
        with tf.device('/CPU'):
            dictionary_preprocess= nn4cast.predefined_classes.Preprocess(dictionary_hyperparams=hyperparameters)
            outputs = nn4cast.predefined_classes.Model_build_and_test(dictionary_hyperparams= hyperparameters, dictionary_preprocess= dictionary_preprocess, cross_validation= False, n_cv_folds=0)
            # outputs= nn4cast.predefined_classes.Model_build_and_test(dictionary_hyperparams= hyperparameters, dictionary_preprocess= dictionary_preprocess, cross_validation= True, n_cv_folds=4, plot_differences=False, importances=True, region_importances=[[50,65],[-25,-10]])
            nn4cast.predefined_classes.Results_plotter(hyperparameters, dictionary_preprocess, rang_x=2.5, rang_y=10, predictions=outputs['predictions'], observations=outputs['observations'], years_to_plot=[1999,2000], plot_with_contours=False, importances=outputs['importances'], region_importances=outputs['region_attributed'])
            # eof_analysis = nn4cast.predefined_classes.PC_analysis(hyperparameters, outputs_cross_validation['predictions'], outputs_cross_validation['observations'], n_modes=4, n_clusters=3, cmap='RdBu_r')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
