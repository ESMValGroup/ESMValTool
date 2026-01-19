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
import calendar
import numpy as np
import dask.array as da

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

            output = iris.cube.CubeList([])
            for month in range(1,13):

                logger.info(f"Working with year {year}, month {month}")
        
                # Load orginal data in an indendent function
                lai_cube = load_dataset(in_dir, var, cfg, year, month)
                
                # Regrdding
                # uses nearest neighbour, skips if resolution = None
                # This uses a huge amount of resource - be careful
                resolution = cfg["Parameters"]["custom"]["regrid_resolution"]
                if resolution == "None":
                    logger.info("No regridding")
                else:
                    logger.info(f"Regridding {cfg["Parameters"]["custom"]["regrid_resolution"]}")
                    lai_cube = regrid(
                        lai_cube, cfg["Parameters"]["custom"]["regrid_resolution"], "nearest"
                        )

                
                
                # make a daily version with Nan cubes for missing days
                # This will work with 10-day CDS data and 5-day CCI data in updates at a later date
                days_in_month = calendar.monthrange(year, month)[1]
                time_coord = lai_cube.coord('time')
                time_values = time_coord.points
                dts = time_coord.units.num2date(time_values)
                days = [item.day for item in dts]

                # lai cube 0 is the problem, need the zeroth time step form the cuble of 3 time! #####################
                for day in range(1, days_in_month + 1):
                    if day in days:
                        logger.info(f"{day} is in CUBES")
                        #iris.util.new_axis(lai_cube, 'time')
                        new_cube = iris.util.new_axis(lai_cube[days.index(day)], 'time')
                        output.append(new_cube)
                    else:
                        logger.info(f"{day} NOT in CUBES")
                        # nan_cube = _create_nan_cube(lai_cube[0], year, month, day)
                        logger.info(f"{lai_cube[0]=}")
                        nan_cube = create_dask_cube(lai_cube[0], year, month, day)
                        new_cube = iris.util.new_axis(nan_cube, 'time')
                        output.append(new_cube)

                logger.info(f"{output=}")

                output = output.concatenate_cube()
                logger.info(f"{output=}")
                
                # time bounds
                # This sets time bounds without needing extra loops and checks
                output.coord('time').guess_bounds()

                # cmorize
                output = _cmorize_dataset(output, var, cfg)
                #logger.info(f"********{lai_cube=}")
                #print(0/0)
                # save cube
                logger.info(f"Saving CMORized cube for variable {output.var_name}")
                # these should all be the same
                attributes = cfg["attributes"]
                attributes["mip"] = var["mip"]
                utils.save_variable(output, lai_cube.var_name, out_dir, attributes, zlib=True)
                #print(0/0)

# from CCI SNOW CMORISER
def _create_nan_cube(cube, year, month, day):
    """Create cube containing only nan from existing cube."""
    nan_cube = cube.copy()
    nan_cube.data = np.full_like(nan_cube.data, np.nan, dtype=np.float32)

    # Read dataset time unit and calendar from file
    dataset_time_unit = str(nan_cube.coord("time").units)
    dataset_time_calender = nan_cube.coord("time").units.calendar

    # Convert datetime
    newtime = datetime.datetime(year=year, month=month, day=day)
    newtime = cf_units.date2num(
        newtime, dataset_time_unit, dataset_time_calender
    )

    nan_cube.coord("time").points = np.float64(newtime)

    return nan_cube

def create_dask_cube(cube, year, month, day):
    nan_da = da.full(cube.shape, np.nan,
                     chunks='auto', dtype=np.float32)
    
    new_cube = cube.copy()
    new_cube.data = nan_da

    dataset_time_unit = str(new_cube.coord("time").units)
    dataset_time_calender = new_cube.coord("time").units.calendar

    # Convert datetime
    newtime = datetime(year=year, month=month, day=day)
    newtime = cf_units.date2num(
        newtime, dataset_time_unit, dataset_time_calender
    )

    new_cube.coord("time").points = np.float64(newtime)

    return new_cube