"""ESMValTool CMORizer for cds-satellite-lai-fapar data.
Tier
   Tier 3
Source
   https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-lai-fapar
Last access
   20190703
Download and processing instructions
   day 20 of each month
Modification history
   20190703-A_crez_ba: written.
"""


import logging
#from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from datetime import datetime
#from os import cpu_count
from pathlib import Path
from warnings import catch_warnings, filterwarnings
from esmvalcore.preprocessor import regrid
import cf_units

import iris
import numpy as np
import os
import glob

from esmvalcore.cmor.table import CMOR_TABLES

from . import utilities as utils

logger = logging.getLogger(__name__)

#TODO where to set this parameter?
regridres = '0.25x0.25'

def _cmorize_dataset(in_file, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input file '%s'",
                var['short_name'], in_file)
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    with catch_warnings():
        filterwarnings(
            action='ignore',
            message="Ignoring netCDF variable 'tcc' invalid units '(0 - 1)'",
            category=UserWarning,
            module='iris',
        )
        cube = iris.load_cube(
            str(in_file),
            constraint=utils.var_name_constraint(var['raw']),
        )

    # Set correct names
    cube.var_name = definition.short_name
    cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    for coord_name in 'latitude', 'longitude', 'time':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float64') #TODO @bouwe needed?  
   
    # Convert units if required
    cube.convert_units(definition.units)

    # Make latitude increasing
    cube = cube[:, ::-1, ...]

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    logger.info("Saving cube\n%s", cube)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)

    return in_file

def _regrid_cds_satellite_lai_fapar(filepattern,var_info):
    '''regrid each file and write to disk appending 'regrid' in front of filename'''
    filelist = glob.glob(filepattern)
    for infile in filelist:
        infile_head,infile_tail = os.path.split(infile)
        outfile_tail = infile_tail.replace('c3s','c3s_regridded')
        outfile = os.path.join(infile_head,outfile_tail)
        logger.info("Reading: {0}".format(infile))
        lai_cube = iris.load_cube(infile,var_info.standard_name)
        logger.info("Regridding to: {0}".format(regridres))
        lai_cube = regrid(lai_cube,regridres,'nearest')
        logger.info("Saving: {0}".format(outfile))
        iris.save(lai_cube,outfile)

def _concatenate_cds_satellite_lai_fpar_over_time(filepattern,var_info):
    """Concatenates single files over time and returns on single cube"""
    filelist = glob.glob(filepattern)
    cubelist = iris.load(filelist,var_info.standard_name)

    # For saving the identifiers
    identifiers = []
    for n in range(len(cubelist)):
        logger.info("Fixing time bounds for cube {0}".format(n+1))
        time_coverage_start = cubelist[n].attributes.pop('time_coverage_start')
        time_coverage_end = cubelist[n].attributes.pop('time_coverage_end')
        # The identifier is a combination of the parent_identifier (which is not removed) and
        identifier = cubelist[n].attributes.pop('identifier')
        # Keep it in a list
        identifiers.append(identifier)

        # Now put time_coverage_start/end as time_bnds
        # Convert time_coverage_xxxx to datetime
        bnd_a = datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%SZ")
        bnd_b = datetime.strptime(time_coverage_end, "%Y-%m-%dT%H:%M:%SZ")

        # Put in shape for time_bnds
        time_bnds_datetime = [bnd_a,bnd_b]

        # Read dataset time unit and calendar from file
        dataset_time_unit = str(cubelist[n].coord('time').units)
        dataset_time_calender = cubelist[n].coord('time').units.calendar
        # Convert datetime
        time_bnds = cf_units.date2num(time_bnds_datetime, dataset_time_unit,dataset_time_calender)
        # Put them on the file
        cubelist[n].coord('time').bounds = time_bnds
    
    # Now the cubes can be concatenated over the time dimension
    cube = cubelist.concatenate_cube()

    # Add identifiers from each cube to the concatenated cube as a comma separated list
    cube.attributes['identifiers_comma_separated']=','.join(identifiers)
    return cube

def cmorization(in_dir, out_dir, cfg):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']


    # run the cmorization#TODO maybe parallel, see era5 example
    for short_name, var in cfg['variables'].items():
        var['short_name'] = short_name
        # First collect all information
        inpfile = os.path.join(in_dir, var['file'])
        logger.info("CMORizing var %s from file %s", short_name, inpfile)
        var_info = cmor_table.get_variable(var['mip'], short_name)
        raw_info = {'name': var['raw'], 'file': inpfile}
        glob_attrs['mip'] = var['mip']

        # Now start regridding
        # Determine filepattern for all files that need regridding
        regrid_filepattern = os.path.join(in_dir, var['file'])
        #_regrid_cds_satellite_lai_fapar(regrid_filepattern,var_info)

        # Specify the filepattern of the regridded files
        concatenate_filepattern = regrid_filepattern.replace('c3s','c3s_regridded')

        result_cube = _concatenate_cds_satellite_lai_fpar_over_time(concatenate_filepattern,var_info)

        # Write it to disk
        savename = os.path.join(in_dir, var_info.short_name+'_regridded.nc')
        logger.info("saving as: {0}".format(savename))
        iris.save(result_cube,savename)

        with catch_warnings():
            filterwarnings(
                action='ignore',
                message=('WARNING: missing_value not used since it\n'
                         'cannot be safely cast to variable data type'),
                category=UserWarning,
                module='iris',
            )
            print(var_info, raw_info, out_dir, glob_attrs)
            in_file = savename
            _cmorize_dataset(in_file, var, cfg, out_dir)


