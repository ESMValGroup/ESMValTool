"""Python example diagnostic using general object-based functions."""
import logging
import os

import iris

import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Compute the time average for each input dataset."""
    datasets = e.Datasets(cfg)
    variables = e.Variables(cfg)
    logger.debug("Found datasets: %s", datasets)
    logger.debug("Found variables: %s", variables)

    for path in datasets:
        logger.info("Processing variable %s from dataset %s",
                    datasets.get_info(n.STANDARD_NAME, path),
                    datasets.get_info(n.DATASET, path))

        logger.debug("Loading %s", path)
        cube = iris.load_cube(path)

        logger.debug("Running example computation")
        cube = cube.collapsed(n.TIME, iris.analysis.MEAN)

        name = os.path.splitext(os.path.basename(path))[0] + '_mean'
        if cfg[n.WRITE_NETCDF]:
            filepath = os.path.join(cfg[n.WORK_DIR], name + '.nc')
            logger.debug("Saving analysis results to %s", filepath)
            iris.save(cube, target=filepath)

        if cfg[n.WRITE_PLOTS] and cfg.get('quickplot'):
            filepath = os.path.join(cfg[n.PLOT_DIR],
                                    name + '.' + cfg[n.OUTPUT_FILE_TYPE])
            logger.debug("Plotting analysis results to %s", filepath)
            e.plot.quickplot(cube, filename=filepath, **cfg['quickplot'])


if __name__ == '__main__':

    with e.run_diagnostic() as config:
        main(config)
