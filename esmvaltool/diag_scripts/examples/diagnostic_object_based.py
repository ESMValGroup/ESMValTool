"""Python example diagnostic using general object-based functions."""
from pprint import pprint
import logging
import os

import iris

import esmvaltool.diag_scripts.shared as e

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Compute the time average for each input model."""
    models = e.Models(cfg)
    variables = e.Variables(cfg)
    logger.debug("Found models: %s", models)
    logger.debug("Found variables: %s", variables)

    pprint(cfg)

    for path in models:
        logger.info("Processing variable %s from model %s",
                    models.get_standard_name(path),
                    models.get_model(path))

        logger.debug("Loading %s", path)
        cube = iris.load_cube(path)

        logger.debug("Running example computation")
        cube = cube.collapsed(e.TIME, iris.analysis.MEAN)

        name = os.path.splitext(os.path.basename(path))[0] + '_mean'
        if cfg[e.WRITE_NETCDF]:
            filepath = os.path.join(cfg[e.WORK_DIR], name + '.nc')
            logger.debug("Saving analysis results to %s", filepath)
            iris.save(cube, target=filepath)

        if cfg[e.WRITE_PLOTS] and cfg.get('quickplot'):
            filepath = os.path.join(cfg[e.PLOT_DIR],
                                    name + '.' + cfg[e.OUTPUT_FILE_TYPE])
            logger.debug("Plotting analysis results to %s", filepath)
            e.plot.quickplot(cube, filename=filepath, **cfg['quickplot'])


if __name__ == '__main__':

    with e.run_diagnostic() as config:
        main(config)
