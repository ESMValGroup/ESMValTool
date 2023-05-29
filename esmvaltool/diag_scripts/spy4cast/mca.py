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

    # Example of how to loop over variables/datasets in alphabetical order
    # groups = group_metadata(input_data, 'variable_group', sort='dataset')
    # for group_name in groups:
    #     logger.info("Processing variable %s", group_name)
    #     for attributes in groups[group_name]:
    #         logger.info("Processing dataset %s", attributes['dataset'])
    #         input_file = attributes['filename']
    #         ds = Dataset(input_file, cfg['work_dir']).open(attributes['short_name'])
    #         map_clim = Clim(ds, 'map')
    #         map_clim.save('map_climatology', cfg['plot_dir'])
    #         map_clim = Clim.load('map_climatology', cfg['plot_dir'], type='map')
    #         map_clim.plot(
    #             show_plot=False,
    #             save_fig=True,
    #             name=os.path.join(cfg['plot_dir'], f'clim-{attributes["dataset"]}.png'),
    #             cmap='jet',
    #             levels=list(range(15, 36)),
    #             ticks=list(range(15, 36)),
    #         )
    variables = group_metadata(input_data, 'variable_group', sort='dataset')
    single_predictand = False
    if len(variables['predictand']) == 1:
        predictand = variables['predictand'][0]
        z = Dataset(name=predictand['filename'], dir=cfg['work_dir']).open(var=predictand['short_name'])
        z.slice(Region(35, 45, 345, 360, Month.JUN, Month.JUL, 1850, 2000), skip=0)
        z_ppcessed = spy4cast.Preprocess(z)
        z_ppcessed.save(f'preprocessed_z_{predictand["alias"]}', dir=cfg['work_dir'])
        single_predictand = True

    groups = group_metadata(input_data, 'dataset', sort='dataset')
    for dataset in groups:
        logger.info("Processing dataset %s", dataset)
        variables = group_metadata(groups[dataset], 'variable_group')
        if 'predictor' not in variables:
            continue
        predictor = variables['predictor'][0]
        predictand = variables['predictand'][0]

        y = Dataset(name=predictor['filename'], dir=cfg['work_dir']).open(var=predictor['short_name'])
        y.slice(Region(-5, 10, -30, -10, Month.JUN, Month.JUL, 1850, 2000), skip=0)
        y_ppcessed = spy4cast.Preprocess(y)
        y_ppcessed.save('save_preprocessed_y_', dir=cfg['work_dir'])

        if not single_predictand:
            predictand = variables['predictand'][0]
            z = Dataset(name=predictand['filename'], dir=cfg['work_dir']).open(var=predictand['short_name'])
            z.slice(Region(35, 45, 345, 360, Month.JUN, Month.JUL, 1850, 2000), skip=0)
            z_ppcessed = spy4cast.Preprocess(z)
            z_ppcessed.save(f'preprocessed_z_{predictand["alias"]}', dir=cfg['work_dir'])

        nm = 3
        alpha = 0.1
        mca = spy4cast.MCA(y_ppcessed, z_ppcessed, nm, alpha)
        mca.save(f'mca_{dataset}', dir=cfg['work_dir'])
        mca.plot(save_fig=True, cmap='viridis', dir=cfg['plot_dir'], name=f'mca-{dataset}.png')

        if cfg['cross_valdiation']:
            cross = spy4cast.Crossvalidation(y_ppcessed, z_ppcessed, nm, alpha)
            cross.save(f'cross_{dataset}', dir=cfg['work_dir'])
            cross.plot(save_fig=True, dir=cfg['plot_dir'], name=f'crossvalidation-{dataset}.png')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
