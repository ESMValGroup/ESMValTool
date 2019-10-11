"""marrmot diagnostic."""
import logging
from pathlib import Path

# import dask.array as da
import iris

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(ancestor_file):
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the marrmot hydrological model.",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            'camphuijsen_jaro',
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
    """Process data for use as input to the marrmot hydrological model """
    input_data = cfg['input_data'].values()
    logger.info(input_data)
    grouped_input_data = group_metadata(input_data,
                                        'standard_name',
                                        sort='dataset')

    # for now just open and save the input/output
    for standard_name in grouped_input_data:
        logger.info("Processing variable %s", standard_name)
        for attributes in grouped_input_data[standard_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            cube = iris.load_cube(input_file)

            # Do stuff
            # - Marrmot is a collection of lumped models, so we need
            #   accumulated values of each variable.
            # - Unit conversion: P = mm/d, PET = mm/d, T = K
            # - Output format:
            # >> A Matlab structure with fields 'precip', 'pet', 'temp',
            # >> and 'delta_t'. 'precip', 'pet', and 'temp' are vectors
            # >> of size 1x[length of time period] each. 'delta_t' is a
            # >> scalar of size 1x1 that specifies the time resolution
            # >> of the 'precip', 'pet' and 'temp' vectors in units [days].
            # >>
            # >> A vector with initial values for each of the model stores,
            # >> of size 1x[number of stores].

            # Save data
            output_file = get_diagnostic_filename(
                Path(input_file).stem + '_marrmot', cfg)
            iris.save(cube, output_file, fill_value=1.e20)

            # Store provenance
            provenance_record = get_provenance_record(input_file)
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
