"""HYPE diagnostic."""
import logging
from pathlib import Path

import dask.array as da
import iris
import numpy
import pandas

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(attributes):
    """Create a provenance record."""
    ancestor_file = attributes['filename']

    record = {
        'caption': "Forcings for the Hype hydrological model.",
        'domains': ['global'],
        'authors': [
            'pelupessy_inti',
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': [ancestor_file],
    }
    return record


def get_output_stem(attributes):
    """Get output file stem, specific to HYPE."""
    short_to_stem = dict(tas="Tobs",
                         tasmin="TMINobs",
                         tasmax="TMAXobs",
                         pr="Pobs")

    shortname = attributes["short_name"]
    if shortname in short_to_stem:
        stem = Path(short_to_stem[shortname])
    else:
        stem = Path(attributes['filename']).stem + '_hype'

    stem = attributes['alias'] / stem

    return stem


def get_data_times_and_ids(attributes):
    """Get the data table to be written and the times and indices."""
    input_file = attributes['filename']

    cube = iris.load_cube(input_file)

    data = numpy.array(cube.core_data())

    # Round times to integer number of days
    time_coord = cube.coord('time')
    time_coord.points = da.floor(time_coord.core_points())
    time_coord.bounds = None

    times = [x.point.strftime("%Y-%m-%d") for x in time_coord.cells()]
    ids = cube.coord('shape_id').core_points()

    return data, times, ids


def main(cfg):
    """Process data for use as input to the HYPE hydrological model."""
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data,
                                        'long_name',
                                        sort='dataset')

    for long_name in grouped_input_data:
        logger.info("Processing variable %s", long_name)
        for attributes in grouped_input_data[long_name]:
            logger.info("Processing dataset %s", attributes['dataset'])

            output_file = get_diagnostic_filename(get_output_stem(attributes),
                                                  cfg, 'txt')
            Path(output_file).parent.mkdir(exist_ok=True)

            data, times, ids = get_data_times_and_ids(attributes)

            frame = pandas.DataFrame(data, index=times, columns=ids)

            frame.to_csv(output_file,
                         sep=' ',
                         index_label="DATE",
                         float_format='%.3f')

            # Store provenance
            provenance_record = get_provenance_record(attributes)
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
