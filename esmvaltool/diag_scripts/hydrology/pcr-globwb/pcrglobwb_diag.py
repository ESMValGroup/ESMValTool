"""PCR-GLOBWB diagnostic."""
import datetime
import logging
import os
from pprint import pformat

import iris

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (ProvenanceLogger,
                                                  get_diagnostic_filename,
                                                  get_plot_filename)

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Average {long_name} between {start_year} and {end_year} "
               "according to {dataset}.".format(**attributes))

    record = {
        'caption': caption,
        'statistics': ['other'],
        'domains': ['global'],
        'authors': ['ande_bo'],
        'references': [' acknow_project'],
    }

    return record


def round_time(cube):
    times = [cell.point for cell in cube.coord('time').cells()]
    rounded_times = [
        datetime.datetime(t.year, t.month, t.day, 0, 0, 0) for t in times
    ]
    cube.coord('time').points = [
        cube.coord('time').units.date2num(cl) for cl in rounded_times
    ]
    cube.coord('time').bounds = None

    return cube


def main(cfg):
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data,
                                        'standard_name',
                                        sort='dataset')
    for standard_name in grouped_input_data:
        logger.info("Processing variable %s", standard_name)
        for attributes in grouped_input_data[standard_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            # Change dimensions order to match PCR-GLOBWB requirement
            cube = iris.load_cube(input_file)
            cube = round_time(cube)
            # Set lat from highest to lowest value
            cube = cube[:, ::-1, ...]
            output_basename = os.path.splitext(
                os.path.basename(input_file))[0] + '_pcrglobwb.nc'
            iris.save(cube, output_basename, fill_value=1E20)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
