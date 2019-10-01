"""PCR-GLOBWB diagnostic."""
import logging
from pathlib import Path

import dask.array as da
import iris

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(ancestor_file):
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the PCR-GLOBWB hydrological model.",
        'domains': ['global'],
        'authors': [
            'aerts_jerom',
            'andela_bouwe',
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


def main(cfg):
    """Process data for use as input to the PCR-GLOBWB hydrological model."""
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data,
                                        'standard_name',
                                        sort='dataset')

    for standard_name in grouped_input_data:
        logger.info("Processing variable %s", standard_name)
        for attributes in grouped_input_data[standard_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            cube = iris.load_cube(input_file)

            # Round times to integer number of days
            time_coord = cube.coord('time')
            time_coord.points = da.floor(time_coord.core_points())
            time_coord.bounds = None

            # Set lat from highest to lowest value
            cube = cube[:, ::-1, ...]

            # Save data
            output_file = get_diagnostic_filename(
                Path(input_file).stem + '_pcrglobwb', cfg)
            iris.save(cube, output_file, fill_value=1.e20)

            # Store provenance
            provenance_record = get_provenance_record(input_file)
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
