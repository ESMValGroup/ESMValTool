"""Python example diagnostic."""
import logging
from pathlib import Path
from pprint import pformat

import iris

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
)
from esmvaltool.diag_scripts.shared.plot import quickplot

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
            'andela_bouwe',
            'righi_mattia',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    cube = iris.util.squeeze(cube)
    return cube


def plot_diagnostic(cube, basename, provenance_record, cfg):
    """Create diagnostic data and plot it."""

    # Save the data used for the plot
    save_data(basename, provenance_record, cfg, cube)

    if cfg.get('quickplot'):
        # Create the plot
        quickplot(cube, **cfg['quickplot'])
        # And save the plot
        save_figure(basename, provenance_record, cfg)


def main(cfg):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    # Demonstrate use of metadata access convenience functions.
    selection = select_metadata(input_data, short_name='tas', project='CMIP5')
    logger.info("Example of how to select only CMIP5 temperature data:\n%s",
                pformat(selection))

    selection = sorted_metadata(selection, sort='dataset')
    logger.info("Example of how to sort this selection by dataset:\n%s",
                pformat(selection))

    grouped_input_data = group_metadata(input_data,
                                        'variable_group',
                                        sort='dataset')
    logger.info(
        "Example of how to group and sort input data by variable groups from "
        "the recipe:\n%s", pformat(grouped_input_data))

    # Example of how to loop over variables/datasets in alphabetical order
    groups = group_metadata(input_data, 'variable_group', sort='dataset')
    for group_name in groups:
        logger.info("Processing variable %s", group_name)
        for attributes in groups[group_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            cube = compute_diagnostic(input_file)

            output_basename = Path(input_file).stem
            if group_name != attributes['short_name']:
                output_basename = group_name + '_' + output_basename
            provenance_record = get_provenance_record(
                attributes, ancestor_files=[input_file])
            plot_diagnostic(cube, output_basename, provenance_record, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
