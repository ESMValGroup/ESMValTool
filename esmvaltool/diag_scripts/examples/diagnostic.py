"""Python example diagnostic."""
import logging
import os

import iris

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Compute the time average for each input model."""
    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from model %s",
                    attributes['standard_name'], attributes['model'])

        logger.debug("Loading %s", filename)
        cube = iris.load_cube(filename)

        logger.debug("Running example computation")
        cube = cube.collapsed('time', iris.analysis.MEAN)

        name = os.path.splitext(os.path.basename(filename))[0] + '_mean'
        if cfg['write_netcdf']:
            path = os.path.join(
                cfg['work_dir'],
                name + '.nc',
            )
            logger.debug("Saving analysis results to %s", path)
            iris.save(cube, target=path)

        if cfg['write_plots'] and cfg.get('quickplot'):
            path = os.path.join(
                cfg['plot_dir'],
                name + '.' + cfg['output_file_type'],
            )
            logger.debug("Plotting analysis results to %s", path)
            quickplot(cube, filename=path, **cfg['quickplot'])


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
