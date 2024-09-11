"""Diagnostic script to plot minima and maxima trends.

based on code from Anton Steketee's COSIMA cookbook notebook
https://cosima-recipes.readthedocs.io/en/latest/Examples/Sea_Ice_Area_Concentration_Volume_with_Obs.html
"""

import logging
import os

import iris
import matplotlib.pyplot as plt
from iris import quickplot

from esmvaltool.diag_scripts.shared import run_diagnostic, save_figure

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def _prom_dim_coord(cube, _field, _filename):
    iris.util.promote_aux_coord_to_dim_coord(cube, 'year')


def plot_trends(datagroup, provenance_record, cfg):
    """Create plot for min and max groups"""
    for variable_group, attributes in datagroup.items():
        plt.clf()
        for (fp, sn, dt, offset) in attributes:
            cube = iris.load_cube(fp, sn, _prom_dim_coord)
            if offset:
                cube.coord('year').points = [y + offset for y in
                                             cube.coord('year').points]
            quickplot.plot(cube, label=dt)

        plt.title(f"Trends in Sea-Ice {variable_group.split('_')[1]}ima")
        plt.legend(loc='upper left')
        plt.ylabel('Sea-Ice Area (km2)')

        save_figure(variable_group, provenance_record, cfg, dpi=300)


def main(cfg):
    """Compute sea ice area for each input dataset."""
    provenance_record = {
        'caption': "sea ice trends southern hemisphere",
        'authors': [
            'chun_felicity',
            'steketee_anton'
        ],
        'references': [''],
        'ancestors': list(cfg['input_data'].keys()),
    }
    input_data = cfg['input_data'].values()

    datagroup = {}  # for each variable min and max

    for dataset in input_data:
        # Load the data
        if 'offset_years' in dataset:
            input_file = (dataset['filename'], dataset['short_name'],
                          dataset['dataset'], dataset['offset_years'])
        else:
            input_file = (dataset['filename'], dataset['short_name'],
                          dataset['dataset'], None)
        # key for different models
        logger.info(f"dataset: {dataset['long_name']}")
        if dataset['variable_group'] not in datagroup:
            datagroup[dataset['variable_group']] = []
        datagroup[dataset['variable_group']].append(input_file)

    logger.info(datagroup)
    plot_trends(datagroup, provenance_record, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
