"""ESMValTool CMORizer for cds-satellite-lai-fapar data.

Tier
   Tier 3
Source
   https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-lai-fapar?tab=form
Last access
   20190703

NEEDED TO UPDATE THIS!!!
Download and processing instructions
   - Open in a browser the data source as specified above
   - Put the right ticks:
      - Tick variables LAI and FAPAR
      - Tick satellite SPOT (System Pour l'Observation de la Terre)
      - Tick sensor VGT (Vegetation)
      - Tick horizontal resolution 1km
      - Tick product version V1
      - Tick all available years
      - Tick all available months
      - Tick Nominal day 20
   - Click 'submit form'
   - According to ESMValTool practice, put them in the right rawobsdir folder

Notes
-----
   - This script regrids and cmorizes the above dataset.
   - Request might need to be split into chunks to not exceed download limit

Caveats
   - Fails setting standard name for variable FAPAR

Modification history
   20200512-crezee_bas: adapted to reflect changes in download form by CDS.
   20190703-crezee_bas: written.
"""

import glob
import logging
import os
from copy import deepcopy
from datetime import datetime
from warnings import catch_warnings, filterwarnings

import cf_units
import iris
from esmvalcore.preprocessor import regrid
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

def load_callback(cube, field, filename):
    """
    Callback fucntion for iris.load to remove all attributes from cube
    so they will concatenate into a single cube
    """
    cube.attributes = None

def load_dataset(in_dir, var, cfg, year, month):
    """
    Load the files from an individual month
    """
    filelist = glob.glob(os.path.join(in_dir, var["file"]))
    this_month_year_files = []
    for file in filelist:
       if f"{year}{month:02d}" in file:
           this_month_year_files.append(file)
                
    lai_cube = iris.load(this_month_year_files,
                        NameConstraint(var_name=var["raw"]),
                        callback=load_callback,
            )
    
    return lai_cube.concatenate_cube()

def _cmorize_dataset(cube, var, cfg):

    cmor_table = cfg["cmor_table"]
    definition = cmor_table.get_variable(var["mip"], var["short_name"])
    
    # standard name
    # long name
    cube.var_name = definition.short_name
    if definition.standard_name:
        cube.standard_name = definition.standard_name

    cube.long_name = definition.long_name

    # units
    cube.convert_units(definition.units)

    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # run the cmorization
    # Pass on the workdir to the cfg dictionary
    cfg["work_dir"] = cfg_user.work_dir
    # If it doesn't exist, create it
    if not os.path.isdir(cfg["work_dir"]):
        logger.info(
            "Creating working directory for regridding: %s", cfg["work_dir"]
        )
        os.mkdir(cfg["work_dir"])

    logger.info(f"{cfg=}")
    logger.info(f"{cfg['Parameters']['custom']['regrid_resolution']=}")
    
    for short_name, var in cfg["variables"].items():
        var["short_name"] = short_name
        logger.info("Processing var %s", short_name)

        for year in range(cfg["attributes"]["start_year"],
                         cfg["attributes"]["end_year"]):

            for month in range(1,13):

                logger.info(f"Working with year {year}, month {month}")
        
                # Load orginal data in an indendent function
                lai_cube = load_dataset(in_dir, var, cfg, year, month)

                # Regrdding
                # uses nearest neighbour, skips if resolution = None
                resolution = cfg["Parameters"]["custom"]["regrid_resolution"]
                if resolution == "None":
                    logger.info("No regridding")
                else:
                    logger.info(f"Regridding {cfg["Parameters"]["custom"]["regrid_resolution"]}")
                    lai_cube = regrid(
                        lai_cube, cfg["Parameters"]["custom"]["regrid_resolution"], "nearest"
                        )

                # time bounds
                # This sets time bounds without needing extra loops and checks
                lai_cube.coord('time').guess_bounds()

                # cmorize
                lai_cube = _cmorize_dataset(lai_cube, var, cfg)

                # save cube
                logger.info(f"Saving CMORized cube for variable {lai_cube.var_name}")
                # these should all be the same
                attributes = cfg["attributes"]
                attributes["mip"] = var["mip"]
                utils.save_variable(lai_cube, lai_cube.var_name, out_dir, attributes)