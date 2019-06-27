# pylint: disable=invalid-name
"""ESMValTool CMORizer for Duveiller2018 data.

Tier
   Tier 2: other freely-available dataset.

Source
   https://www.nature.com/articles/sdata201814

Last access
   20190430

Download and processing instructions
   Download the dataset albedo_IGBPgen.nc

Modification history
   20190430-A_crez_ba: started with cmorize_obs_Landschuetzer2016.py as an example to follow

"""

import logging
import os
from warnings import catch_warnings, filterwarnings

import iris
from dask import array as da

from .utilities import (set_global_atts, fix_coords, fix_var_metadata,
                        read_cmor_config, save_variable, constant_metadata, convert_timeunits)

logger = logging.getLogger(__name__)

# read in CMOR configuration

#CFG = read_cmor_config('Duveiller2018')



# TODO: maybe not needed?
# pylint: disable=unused-argument
def _fix_fillvalue(cube, field, filename):
    """Create masked array from missing_value."""
    if hasattr(field.cf_data, 'missing_value'):
        # fix for bad missing value definition
        cube.data = da.ma.masked_equal(cube.core_data(),
                                       field.cf_data.missing_value)



#def fix_custom_coord_iTr()

def duveiller2018_callback_function(cube,field,filename):
    cube.coord('Month').rename('time')
    cube.coord('time').units = 'months since 2010-01-01 00:00:00'
    #cube.coord('Vegetation transition code').rename('vegetation_transition_code')

def extract_variable(var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    var = var_info.short_name
    with catch_warnings():
        filterwarnings(
            action='ignore',
            message='Ignoring netCDF variable .* invalid units .*',
            category=UserWarning,
            module='iris',
        )
        cubes = iris.load(raw_info['file'], callback=duveiller2018_callback_function)
    rawvar = raw_info['name']
    print(cubes)
    for cube in cubes:
        if cube.var_name == rawvar:
#            import IPython;IPython.embed()
            for cubecoord in cube.coords():
                if cubecoord.var_name=='iTr':
#                    cubecoord.standard_name = None # CMOR checker raises: iTr: standard_name should be , not None
                    cubecoord.standard_name = ''  # this script raises: ValueError: '' is not a valid standard_name

            fix_var_metadata(cube, var_info)
            fix_coords(cube)
# BAS: this is how it is treated for variables: 
#            if var_info.standard_name == '':
#                cube.standard_name = None
#            else:
#                cube.standard_name = var_info.standard_name
#            _fix_data(cube, var)
            # Rename Month to time coordinate
            set_global_atts(cube, attrs)
            save_variable(
                cube,
                var,
                out_dir,
                attrs,
                local_keys=['positive'],
#                unlimited_dimensions=['time'],
            )

def cmorization(in_dir, out_dir, cfg):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var, vals in cfg['variables'].items():
        inpfile = os.path.join(in_dir, vals['file'])
        logger.info("CMORizing var %s from file %s", var, inpfile)
        var_info = cmor_table.get_variable(vals['mip'], var)
        print("var = ",var)
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
            extract_variable(var_info, raw_info, out_dir, glob_attrs)

