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

import iris
import numpy as np
import os
import glob

from esmvalcore.cmor.table import CMOR_TABLES

from . import utilities as utils

logger = logging.getLogger(__name__)

#TODO where to set this parameter?
regridres = '0.25x0.25'

def _extract_variable(in_file, var, cfg, out_dir):
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

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    # Fix coordinates
    cube.coord('latitude').var_name = 'lat'
    cube.coord('longitude').var_name = 'lon'

    for coord_name in 'latitude', 'longitude', 'time':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float64')
        coord.guess_bounds()

    # Convert units if required
    cube.convert_units(definition.units)

    #TODO check, needed yes/no? Make latitude increasing
    #cube = cube[:, ::-1, ...]

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    logger.info("Saving cube\n%s", cube)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)

    return in_file

def _regrid_cds_satellite_lai_fapar(filepattern):
    '''regrid each file and write to disk appending 'regrid' in front of filename'''
    filelist = glob.glob(filepattern)
    for infile in filelist:
        infile_head,infile_tail = os.path.split(infile)
        outfile_tail = infile_tail.replace('c3s','c3s_regridded')
        outfile = os.path.join(infile_head,outfile_tail)
        logger.info("Reading: {0}".format(infile))
        #lai_cube = iris.load_cube(infile,['leaf_area_index'])
        logger.info("Regridding to: {0}".format(regridres))
        #lai_cube = esmvalpp.regrid(lai_cube,regridres,'nearest')
        logger.info("Saving: {0}".format(outfile))
        #iris.save(lai_cube,outfile)

def _concatenate_cds_satellite_lai_fpar_over_time(filepattern):
    filelist = glob.glob(filepattern)
    for n in range(len(cubelist)):
        time_coverage_start = cubelist[n].attributes.pop('time_coverage_start')
        time_coverage_end = cubelist[n].attributes.pop('time_coverage_end')
        #TODO don't loose the identifier, keep included in some metadata
        identifier = cubelist[n].attributes.pop('identifier')
        # Now put time_coverage_start/end as time_bnds
        # Convert time_coverage_xxxx to datetime
        bnd_a = datetime.datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%SZ")
        bnd_b = datetime.datetime.strptime(time_coverage_end, "%Y-%m-%dT%H:%M:%SZ")
    
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

def cmorization(in_dir, out_dir, cfg):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']


    # run the cmorization#TODO maybe parallel, see era5 example
    for var, vals in cfg['variables'].items():
        # First do regridding here
        print(var,vals)

        # Determine filepattern for all files that need regridding
        regrid_filepattern = os.path.join(in_dir, vals['file'])
        _regrid_cds_satellite_lai_fapar(regrid_filepattern)


        import IPython;IPython.embed()
        
        inpfile = os.path.join(in_dir, vals['file'])
        logger.info("CMORizing var %s from file %s", var, inpfile)
        var_info = cmor_table.get_variable(vals['mip'], var)
        raw_info = {'name': vals['raw'], 'file': inpfile}
        glob_attrs['mip'] = vals['mip']
        with catch_warnings():
            filterwarnings(
                action='ignore',
                message=('WARNING: missing_value not used since it\n'
                         'cannot be safely cast to variable data type'),
                category=UserWarning,
                module='iris',
            )
            _extract_variable(var_info, raw_info, out_dir, glob_attrs)


