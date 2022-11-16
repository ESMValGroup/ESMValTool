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
        air.mon.mean.nc
        hgt.mon.mean.nc
        rhum.mon.mean.nc
        shum.mon.mean.nc
        uwnd.mon.mean.nc
        vwnd.mon.mean.nc
        omega.mon.mean.nc    
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/
        uwnd.*?.nc
        vwnd.*.nc

    Subdirectory surface/:
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/
        air.mon.mean.nc
        pr_wtr.mon.mean.nc
        slp.mon.mean.nc
        wspd.mon.mean.nc
        rhum.mon.mean.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/
        air.2m.mon.mean.nc
        prate.sfc.mon.mean.nc
        pevpr.sfc.mon.mean.nc (not included yet)
        tmax.2m.mon.mean.nc
        tmin.2m.mon.mean.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/other_gauss/
        tcdc.eatm.mon.mean.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/
        prate.sft.gauss.*.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/other_gauss/
        ulwrf.ntat.gauss.*.nc
       

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
    definition = cmor_table.get_variable(var['mip'], short_name)
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    #print(attributes)
    #print(cmor_table)
    #print(short_name)
    
    # load data
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(str(raw_filepath), NameConstraint(var_name = raw_var))

    utils.set_global_atts(cube, attributes)
    #print("first stop")
    
    _fix_units(cube, definition)
    #print("third stop")

    #utils.fix_var_metadata(cube, definition)
    #print(cube)
    #print(definition)
    utils.fix_var_metadata(cube, cmor_info)
    
    
    # fix time units
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))
    
    cube = _fix_coordinates(cube, definition)
    #print("fourth stop")

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
        short_name = var['short_name']
        print(short_name)
        raw_filenames = sorted(Path(in_dir).rglob(var['file']))
        #print(raw_filenames)
        for raw_filename in raw_filenames:    
            #filepath = Path(in_dir) / raw_filename
            #print(raw_filename)

            _extract_variable(short_name, var, cfg, raw_filename, out_dir)
