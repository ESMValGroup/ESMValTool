"""ESMValTool CMORizer for NCEP-NCAR-R1 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html

Last access
    20220908

Download and processing instructions
    To facilitate the download, the links to the ftp server are provided.
    Since the filenames are sometimes identical across different
    save the data in two subdirectories in input_dir_path.
    Subdirectory pressure/:
    
    ftp://ftp.cdc.noaa.gov/Projects/Datasets/data.ncep.reanalysis/pressure/
    
    ftp://ftp.cdc.noaa.gov/Projects/Datasets/data.ncep.reanalysis/Dailies/pressure/
        
        
    Subdirectory surface/:
      ftp://ftp.cdc.noaa.gov/Datasets/data.ncep.reanalysis/Monthlies/gaussian_grid
        air.2m.mon.mean.nc
        tcdc.eatm.mon.mean.nc
        pr_wtr.eatm.mon.mean.nc

    #Select the section "Pressure" and "Surface" and download the variables
    #listed below. Since raw data on pressure levels and for surface have the
    #same file and variable name, save the data in two different subdirectories
    #"press" and "surf" in input_dir_path.
    #Specify the time range of the data as YEAR1-YEAR2 below, considering only
    #complete years (Jan to Dec).

Caveats

"""

import logging
import re
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from datetime import datetime, timedelta
from pathlib import Path
from warnings import catch_warnings, filterwarnings
from cf_units import Unit

import iris
import numpy as np
from esmvalcore.cmor.table import CMOR_TABLES
from esmvalcore.preprocessor import daily_statistics, monthly_statistics
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _fix_units(cube, definition):
    """Fix issues with the units."""
    if cube.var_name in {'clt'}:
        # Change units from fraction to percentage
        cube.units = definition.units
        cube.data = cube.core_data() * 100.


def _fix_coordinates(cube, definition):
    # fix flipped latitude
    utils.flip_dim_coord(cube, 'latitude')
    utils.fix_dim_coordnames(cube)
    cube_coord = cube.coord('latitude')
    utils.fix_bounds(cube, cube_coord)
    cube_coord = cube.coord('longitude')
    utils.fix_bounds(cube, cube_coord)

    return cube

def _extract_variable(short_name, var, cfg, raw_filepath, out_dir):
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    # load data
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(raw_filepath)

    utils.set_global_atts(cube, attributes)
    print("first stop")
    
    _fix_units(cube, definition)
    print("third stop")

    # fix time units
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))
    
    cube = _fix_coordinates(cube, definition)
    print("fourth stop")

    utils.save_variable(
        cube,
        short_name,
        out_dir,
        attributes,
        unlimited_dimensions=['time'],
    )


   
def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Run CMORizer for NCEP-NCAR-R1."""
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        var['short_name'] = short_name
        raw_filepath = os.path.join(in_dir, var['file'])
        print(raw_filepath)

        _extract_variable(short_name, var, cfg, raw_filepath, out_dir)
