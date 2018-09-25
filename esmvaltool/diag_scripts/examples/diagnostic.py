"""Python example diagnostic."""
import logging
import os
from pprint import pformat

import iris

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    # Demonstrate use of metadata access convenience functions.
    selection = select_metadata(input_data, short_name='pr', project='CMIP5')
    logger.info("Example of how to select only CMIP5 precipitation data:\n%s",
                pformat(selection))

    selection = sorted_metadata(selection, sort='dataset')
    logger.info("Example of how to sort this selection by dataset:\n%s",
                pformat(selection))

    grouped_input_data = group_metadata(
        input_data, 'standard_name', sort='dataset')
    logger.info("Example of how to group and sort input data by standard_name:"
                "\n%s", pformat(grouped_input_data))

    # Example of how to loop over variables/datasets in alphabetical order
    for standard_name in grouped_input_data:
        logger.info("Processing variable %s", standard_name)
        for attributes in grouped_input_data[standard_name]:
            logger.info("Processing dataset %s", attributes['dataset'])

            filename = attributes['filename']
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
