"""LISFLOOD diagnostic."""
import logging
from pathlib import Path

import iris

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(ancestor_file):
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the LISFLOOD hydrological model.",
        'domains': ['global'],
        'authors': [
            'verh_st',
            'kalverla_peter',
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
    """Process data for use as input to the LISFLOOD hydrological model """
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

            # Compute additional variables as input for lisvap
            # Actual vapour pressure computed with cdo by e=6.11*exp(7.5*(d2m-273.15)/(237.3+(d2m-273.15))) expression
            #
            # windspeed: wspd = sqrt(u10*u10+v10*v10)
            # surface_solar_radiation_downwards: ??
            # era-interim-name: surface_solar_radiation_downwards aka ssrd

            # Unit conversion:
            # precipitation should be in mm
            # temperature should be in degree C
            # solar_radiation_downwards should be in J m**-2 day**-1
            # water vapour pressure should be in hPa
            # wind speed should be in m/s

            # Concatenation
            # lisflood wants one file per variable (not per year)

            # Target output file
            # In the cdo example, the output format was:
            # 2m_temperature_final_1990-1996.nc

            # Save data
            output_file = get_diagnostic_filename(
                Path(input_file).stem + '_lisflood', cfg)
            iris.save(cube, output_file, fill_value=1.e20)

            # Store provenance
            provenance_record = get_provenance_record(input_file)
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
