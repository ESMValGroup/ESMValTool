import logging
import sys
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
import spectra_compute

logger = logging.getLogger(Path(__file__).stem)


def main(cfg):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    logger.info('INPUT DATA', input_data)

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
            logger.info("Input data file: ", input_file)

            output_basename = Path(input_file).stem
            attributes['plot_dir'] = cfg['plot_dir']
            attributes['work_dir'] = cfg['work_dir']

            if attributes['variable_group'] == 'ua':
                for pressure_level in [85000., 20000.]:
                    attributes['cube'] = iris.load_cube(input_file,
                                                        iris.Constraint(air_pressure=pressure_level))

                    if pressure_level == 85000.:
                        attributes['varname'] = 'x_wind_850hPa'
                    elif pressure_level == 20000.:
                        attributes['varname'] = 'x_wind_200hPa'

                    logger.info(attributes)
                    # Call Spectra calculations
                    spectra_compute.WKSpectra(cfg, attributes).wkSpaceTime()
                    spectra_compute.WKSpectra(cfg, attributes).SpectraSeason()

            elif attributes['variable_group'] == 'pr':
                attributes['cube'] = iris.load_cube(input_file)
                attributes['varname'] = 'Precipitation'

                logger.info(attributes)
                # Call Spectra calculations
                spectra_compute.WKSpectra(cfg, attributes).wkSpaceTime()
                spectra_compute.WKSpectra(cfg, attributes).SpectraSeason()

            elif attributes['variable_group'] == 'rlut': # Check rlut
                attributes['cube'] = iris.load_cube(input_file)
                attributes['varname'] = 'toa_outgoing_longwave_flux'

                logger.info(attributes)
                # Call Spectra calculations
                spectra_compute.WKSpectra(cfg, attributes).wkSpaceTime()
                spectra_compute.WKSpectra(cfg, attributes).SpectraSeason()


            if group_name != attributes['short_name']:
                output_basename = group_name + '_' + output_basename
            if "caption" not in attributes:
                attributes['caption'] = input_file

            print('output_basename', output_basename)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
