"""ESMValTool CMORizer for FLUXCOM GPP data.

<We will add some useful info here later>
"""
import logging
from pathlib import Path

import iris

from . import utilities as utils

logger = logging.getLogger(__name__)

def cmorization(in_dir, out_dir, cfg, _):
    """Cmorize the dataset."""

    # Get general information from the config file
    attributes = cfg['attributes']
    variables = cfg['variables']

    for short_name, variable_info in variables.items():
        logger.info("CMORizing variable: %s", short_name)

        # 1a. Find the input data (one file for each year)
        filename_pattern = cfg['filename']
        matches = Path(in_dir).glob(filename_pattern)

        for match in matches:
            # 1b. Load the input data
            input_file = str(match)
            logger.info("found: %s", input_file)
            cube = iris.load_cube(input_file)

            # 2. Apply the necessary fixes
            # 2a. Fix/add coordinate information and metadata
            cube.coord('lat').standard_name = 'latitude'
            cube.coord('lon').standard_name = 'longitude'
            utils.fix_coords(cube)

            # 2b. Fix gpp units
            logger.info("Changing units for gpp from gc/m2/day to kg/m2/s")
            cube.data = cube.core_data() / (1000 * 86400)
            cube.units = 'kg m-2 s-1'

            # 2c. Fix metadata
            cmor_table = cfg['cmor_table']
            cmor_info = cmor_table.get_variable(variable_info['mip'], short_name)
            utils.fix_var_metadata(cube, cmor_info)

            # 2d. Add additional metadata
            utils.set_global_atts(cube, attributes)

            # 3. Save the CMORized data
            all_attributes = {**attributes, **variable_info}
            utils.save_variable(cube=cube, var=short_name, outdir=out_dir, attrs=all_attributes)
