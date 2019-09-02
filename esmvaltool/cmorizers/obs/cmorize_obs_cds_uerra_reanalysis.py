"""ESMValTool CMORizer for XXXXXXXXXXXXXXXXX

Tier
   Tier 3
Source
   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-uerra-europe-soil-levels?tab=form
Last access
   20190822

Download and processing instructions
   - Open in a browser the data source as specified above
   - Put the right ticks:
      - Tick Origin UERRA-HARMONIE
      - Tick Variable 'Volumetric soil moisture'
      - Tick Soil level 1, 2, 3
      - Tick all available years
      - Tick all available months
      - Tick all available days
   - Click 'submit form'
   - According to ESMValTool practice, put them in the right rawobsdir folder

Notes
-----
   - This script uses the xesmf regridder, which is not standard included in ESMValTool, 
     install it in the esmvaltool environment:  conda install -c conda-forge xesmf


Modification history
   20190821-A_crez_ba: written.
"""

import glob
import logging
import os
from copy import deepcopy
from datetime import datetime
from warnings import catch_warnings, filterwarnings

import cf_units
import iris
import xarray as xr
import xesmf as xe

from .utilities import (constant_metadata, convert_timeunits, fix_coords,
                        fix_var_metadata, save_variable, set_global_atts)
from esmvalcore.preprocessor import regrid
from esmvalcore.preprocessor._regrid import _stock_cube

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _attrs_are_the_same(cubelist):
    # assume they are the same
    attrs_the_same = True
    allattrs = cubelist[0].attributes
    for key in allattrs:
        try:
            unique_attr_vals = {cube.attributes[key] for cube in cubelist}
        # This exception is needed for valid_range, which is an
        # array and therefore not hashable
        except TypeError:
            unique_attr_vals = {tuple(cube.attributes[key])
                                for cube in cubelist}
        if len(unique_attr_vals) > 1:
            attrs_the_same = False
            print("Different values found for {0}-attribute: {1}"
                  .format(key, unique_attr_vals))
    return attrs_the_same

def _cmorize_dataset(in_file, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input file '%s'",
                var['short_name'], in_file)
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']

    cmor_table = cfg['cmor_table']
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    cube = iris.load_cube(
        str(in_file),
        constraint=utils.var_name_constraint(var['raw']))

    #import IPython;IPython.embed()


    # The following lines are essential before applying the common function fix_coords
    # Convert time calendar from proleptic_gregorian to gregorian
    cube.coord('time').units = cf_units.Unit(cube.coord('time').units.origin,'gregorian')

    # Set standard_names for lat and lon
    cube.coord('lat').standard_name = 'latitude'
    cube.coord('lon').standard_name = 'longitude'
    
    cube = fix_coords(cube)

    # Set correct names
    cube.var_name = definition.short_name

    cube.long_name = definition.long_name

    # Convert units if required
    cube.convert_units(definition.units)

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    logger.info("Saving CMORized cube for variable %s", cube.var_name)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)

    return in_file

def _regrid_dataset(in_dir, var, cfg):
    """
    Regridding of original files.

    This function regrids each file and write to disk appending 'regrid'
    in front of filename.
    """
    #TODO move back workdir as cfg['work_dir']
    workdir = cfg['work_dir']
    filelist = glob.glob(os.path.join(in_dir, var['file']))
    for infile in filelist:
        _, infile_tail = os.path.split(infile)
        outfile_tail = infile_tail.replace('uerra', 'uerra_regridded')
        outfile = os.path.join(workdir, outfile_tail)

        targetgrid = _stock_cube(cfg['custom']['regrid'])
        targetgrid_ds = xr.DataArray.from_iris(targetgrid)
        input_ds = xr.open_dataset(infile)
        # Do renaming for consistency of coordinate names
        input_ds = input_ds.rename({'latitude': 'lat','longitude': 'lon'})
        # TODO check selection of soil level
        input_da = input_ds[var['raw']].isel(level=0)
        logger.info("Regridding... ")
        # A workaround to avoid spreading of nan values,
        # related to Github issue
        constantval = 10
        input_da = input_da + constantval
        assert(int((input_da==0.).sum())==0) # Make sure that there 
                                             # are no zero's in the data, 
                                             # since they will be masked out
        #TODO configure such that weightfile is written into workdir
        regridder = xe.Regridder(input_ds,targetgrid_ds,'bilinear',reuse_weights=True)
        da_out = regridder(input_da)
        da_out = da_out.where(da_out!=0.)
        da_out = da_out - constantval
        
        # Save it.
        logger.info("Saving: %s", outfile)
        da_out.to_netcdf(outfile)


def cmorization(in_dir, out_dir, cfg, cfg_user):
    """Cmorization func call."""
    # run the cmorization
    # Pass on the workdir to the cfg dictionary
    #TODO put back cfg_user['work_dir']
    cfg['work_dir'] = '/net/exo/landclim/PROJECTS/C3S/workdir/cmorize_temp/'
    # If it doesn't exist, create it
    if not os.path.isdir(cfg['work_dir']):
        logger.info("Creating working directory for regridding: %s",
                    cfg['work_dir'])
        os.mkdir(cfg['work_dir'])

    for short_name, var in cfg['variables'].items():
        var['short_name'] = short_name
        logger.info("Processing var %s", short_name)

        # Regridding
        logger.info("Start regridding to: %s",
                    cfg['custom']['regrid'])
        #TODO set back
        #_regrid_dataset(in_dir, var, cfg)
        logger.info("Finished regridding")

        #TODO move this block to a function
        for year in range(2008,2009):
            # File concatenation
            logger.info("Concatenating files over time for year {0}".format(year))

            filelist = glob.glob(cfg['work_dir']+'uerra_regridded_{0}??.nc'.format(year))
            cubelist = iris.load(filelist)
            for cube in cubelist:
                cube.remove_coord('time') # Time has strange values, so use 
                                          # forecast_reference_time instead
                cube.coord('forecast_reference_time').rename('time')

            cube_concatenated = cubelist.concatenate()[0]

            in_file = os.path.join(cfg['work_dir'],'concatenated_{0}.nc'.format(year))
            logger.info("Saving as: {0}".format(in_file))
            iris.save(cube_concatenated,in_file)

            # Read in the full dataset here from 'workdir'
            logger.info("CMORizing")
            logger.info("Start CMORization of file %s", in_file)
            _cmorize_dataset(in_file, var, cfg, out_dir)
            logger.info("Finished processing year %s", year)
