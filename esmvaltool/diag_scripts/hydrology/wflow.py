"""wflow diagnostic."""
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
        'caption': "Forcings for the wflow hydrological model.",
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
    """Process data for use as input to the wflow hydrological model """
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
            # ...

            # There is a preprocessing record on the jupyter server
            # They do the following steps:
            # 1. Convert to daily statistics with timestamp at midnight
            # 2. Convert 2m temperature to sea level temp using a
            #    lapse rate of 6.5 K/km. This thus also requires an
            #    orography file (or geopotential at the surface).
            # 3. Regrid to a target grid for the desired basin. This
            #    target grid is specified by a DEM called e.g.
            #    wflow_dem_Meuse.nc (which is based on a .map file, ask #    Jerom). Some of these can be found on cartesius:
            #    /lustre1/0/wtrcycle/lorentz-workshop/wflow_sbm/model_input/wflow_sbm_meuse/staticmaps/
            #     This regridding cannot use a preprocessor in the
            #     recipe, but it can use the ESMValCore API regridder,
            #     as it is possible to specify a target cube.
            # 4.  Convert sea level temperature back to the new target
            #     grid elevation.
            #
            # Note: The example preprocessing on the jupyter server uses
            #       evaporation. However, the model needs potential
            #       evaporation (according to Jerom). He has a script to
            #       convert it, but that requires several additional
            #       input variables.

            # Save data
            output_file = get_diagnostic_filename(
                Path(input_file).stem + '_wflow', cfg)
            iris.save(cube, output_file, fill_value=1.e20)

            # Store provenance
            provenance_record = get_provenance_record(input_file)
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
