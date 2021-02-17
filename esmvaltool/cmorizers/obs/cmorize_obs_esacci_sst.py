"""ESMValTool CMORizer for ESACCI-SST data.

Tier
   Tier 2: other freely-available dataset.

Source
   http://surftemp.net/regridding/index.html

Last access
   20201214

Download and processing instructions
   Download the following files:
     Go to http://surftemp.net/regridding/index.html
     and request regridded data with the following options:
       Time Resolution: monthly
       Longitude Resolution: 0.5
       Latitude Resolution: 0.5
       Start Date: 1982-01-01
       End Date: 2019-12-31
       Exclude data above sea ice threshold: True
         (Threshold: 100 %)
       Include post-hoc SST bias adjustments: True
       Output Absolute or Anomaly SST: absolute
       Generate Sea Ice Fraction: True
       Error Correlation in Time (Days): 7
       Error Correlation In Space (Degrees): 3.0

Modification history
   20201204-roberts_charles: written.
   20201214-predoi_valeriu: approved.
   20201214-lauer_axel: approved.

"""

import logging
import os

import iris

from esmvalcore.preprocessor import concatenate
from .utilities import (convert_timeunits, fix_coords, fix_var_metadata,
                        save_variable, set_global_atts)

logger = logging.getLogger(__name__)


def extract_variable(var_info, raw_info, attrs, year):
    """Extract to all vars."""
    cubes = iris.load(raw_info['file'])
    rawvar = raw_info['name']

    for cube in cubes:
        if cube.var_name == rawvar:
            fix_var_metadata(cube, var_info)
            convert_timeunits(cube, year)
            fix_coords(cube)
            set_global_atts(cube, attrs)
            return cube


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw'], 'file': vals['file']}
        inpfile = os.path.join(in_dir, cfg['filename'])
        logger.info("CMORizing var %s from file type %s", var, inpfile)
        years = range(1982, 2020)
        months = ["0" + str(mo) for mo in range(1, 10)] + ["10", "11", "12"]
        for year in years:
            monthly_cubes = []
            for month in months:
                raw_info['file'] = inpfile.format(year=year, month=month)
                logger.info("CMORizing var %s from file type %s", var,
                            raw_info['file'])
                cube = extract_variable(var_info, raw_info, glob_attrs, year)
                monthly_cubes.append(cube)
            yearly_cube = concatenate(monthly_cubes)
            save_variable(yearly_cube,
                          var,
                          out_dir,
                          glob_attrs,
                          unlimited_dimensions=['time'])
