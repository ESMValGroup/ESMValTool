"""ESMValTool CMORizer for WOA data.

Tier
   Tier 2: other freely-available dataset.

Source
   WOA18: https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/

Last access
   WOA18: 20210311

Download and processing instructions
   All handled by the script (download only if local data are missing)

Modification history
   20210311-lovato_tomas: consolidate WOA18 and use OBS6
   20200911-bock_lisa: extend to WOA18
   20190328-lovato_tomas: cmorizer revision
   20190131-predoi_valeriu: adapted to v2.
   20190131-demora_lee: written.

"""

import logging
import os
from warnings import catch_warnings, filterwarnings
import requests
import iris
import numpy as np
from cf_units import Unit

from esmvalcore.preprocessor import annual_statistics

from .utilities import (constant_metadata, fix_coords, fix_var_metadata,
                        save_variable, set_global_atts)

logger = logging.getLogger(__name__)


def collect_files(in_dir, var, cfg):
    """Compose input file list and download if missing."""
    file_list = []
    var_dict = cfg['variables'][var]
    in_dir = os.path.join(in_dir, var_dict['name'])
    # create list of monthly data
    for mon in range(1, 13):
        fname = var_dict['file'] + "{:02d}".format(mon) + '_01.nc'
        in_file = os.path.join(in_dir, fname)

        # download if missing
        if not os.path.isfile(in_file):
            if not os.path.isdir(in_dir):
                os.makedirs(in_dir)
            logger.info('Input file %s is missing. Start download ... ', fname)
            url = os.path.join(cfg['attributes']['source'], var_dict['name'],
                               'netcdf', var_dict['file'].split('_')[1],
                               cfg['custom']['resolution'], fname)
            url_file = requests.get(url)
            open(in_file, 'wb').write(url_file.content)

        file_list.append(in_file)

    return file_list


def extract_variable(in_files, out_dir, attrs, raw_info, cmor_table):
    """Extract variables and create OBS dataset."""
    var = raw_info['var']
    var_info = cmor_table.get_variable(raw_info['mip'], var)
    rawvar = raw_info['raw_var']
    with catch_warnings():
        filterwarnings(
            action='ignore',
            message='Ignoring netCDF variable .* invalid units .*',
            category=UserWarning,
            module='iris',
        )
        cubes = iris.load(in_files, rawvar)
    iris.util.equalise_attributes(cubes)
    cube = cubes.concatenate_cube()

    # set reference time
    year = raw_info['reference_year']
    cube.coord('time').points = np.arange(0.5, 12.5)
    cube.coord('time').units = Unit('months since ' + str(year) +
                                    '-01-01 00:00:00',
                                    calendar='gregorian')

    fix_var_metadata(cube, var_info)
    fix_coords(cube)
    set_global_atts(cube, attrs)
    save_variable(cube, var, out_dir, attrs, unlimited_dimensions=['time'])

    # compute derived vars
    if 'srf_var' in raw_info:
        var_info = cmor_table.get_variable(raw_info['mip'],
                                           raw_info['srf_var'])
        logger.info("Extract surface OBS for %s", raw_info['srf_var'])
        level_constraint = iris.Constraint(cube.var_name, depth=0)
        cube = cube.extract(level_constraint)
        fix_var_metadata(cube, var_info)
        save_variable(cube,
                      raw_info['srf_var'],
                      out_dir,
                      attrs,
                      unlimited_dimensions=['time'])

    if raw_info['compute_Oyr']:
        attrs['mip'] = 'Oyr'
        cube = annual_statistics(cube, 'mean')
        save_variable(cube, var, out_dir, attrs, unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        in_files = collect_files(in_dir, var, cfg)
        logger.info("CMORizing var %s from input set %s", var, vals['name'])
        raw_info = cfg['variables'][var]
        raw_info.update({
            'var': var,
            'reference_year': cfg['custom']['reference_year'],
            'compute_Oyr': cfg['custom']['compute_Oyr']
        })
        glob_attrs['mip'] = vals['mip']
        extract_variable(in_files, out_dir, glob_attrs, raw_info, cmor_table)
