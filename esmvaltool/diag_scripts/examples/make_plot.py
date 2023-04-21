"""Python example diagnostic."""
import logging
from pathlib import Path

import cf_units
import iris
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import run_diagnostic, save_figure

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': "Part of figure 9.3a from IPCC AR6.",
        'authors': [
            'kalverla_peter',
            'andela_bouwe',
        ],
        'references': [
            'acknow_project',  # TODO: replace by correct reference
        ],
        'ancestors': ancestor_files,
    }
    return record


def main(cfg):
    """Compute the time average for each input dataset."""
    for attributes in cfg['input_data'].values():
        logger.info("Processing dataset %s", attributes['alias'])
        input_file = attributes['filename']
        cube = iris.load_cube(input_file)

        time_coord = cube.coord('time')
        time_coord.units = cf_units.Unit(time_coord.units.origin,
                                         calendar='gregorian')
        iris.quickplot.plot(cube, label=attributes['alias'].replace('_', ' '))

    plt.legend()

    ancestor_files = list(cfg['input_data'].keys())
    provenance_record = get_provenance_record(ancestor_files)
    filename = 'IPCC_AR6_figure_9.3a_1850-2100'
    save_figure(filename, provenance_record, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
