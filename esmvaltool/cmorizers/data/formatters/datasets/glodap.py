"""ESMValTool CMORizer for GLODAP data.

Tier
   Tier 2: other freely-available dataset.

Source
   GLODAP: https://www.glodap.info/index.php/mapped-data-product/

Last access
   GLODAP: 20210528

Download and processing instructions
   All handled by the script (download only if local raw data are missing)

Modification history
   20210528-lovato_tomas: written
"""

import logging
import os
import tarfile
from warnings import catch_warnings, filterwarnings

import iris
import requests
from cf_units import Unit

from esmvaltool.cmorizers.data.utilities import (
    constant_metadata,
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube):
        if var in [
                'dissic',
                'talk',
        ]:
            cube /= 1000.  # Convert from umol/kg to mol/m^3
    return cube


def collect_files(
    in_dir,
    var,
    cfg,
):
    """Compose input file list and download if missing."""
    var_dict = cfg['variables'][var]

    fname = '.'.join([var_dict['file'], var_dict['raw_var'], 'nc'])
    in_file = os.path.join(in_dir, fname)

    # check if input file is missing
    if not os.path.isfile(in_file):
        if not os.path.isdir(in_dir):
            os.makedirs(in_dir)

        # check if raw tar file is in place
        tar_file = os.path.basename(cfg['attributes']['source'])
        tar_file = os.path.join(in_dir, tar_file)
        if not os.path.isfile(tar_file):
            logger.info('Input file %s is missing\n', tar_file)
            logger.info('Start download (requested space ~250Mb)... ')
            url_file = requests.get(cfg['attributes']['source'])
            open(tar_file, 'wb').write(url_file.content)

        # get input file from tar archive
        tar_file = tarfile.open(name=tar_file, mode='r')
        tar_base = os.path.basename(cfg['attributes']['source'])[:-7]
        member = tar_file.getmember(os.path.join(tar_base, fname))
        member.name = fname
        tar_file.extract(member, path=in_dir)
        tar_file.close()

    return in_file


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
        cube = iris.load_cube(in_files, rawvar)
        depth = iris.load_cube(in_files, 'Depth')

    # add depth coord
    cube.add_dim_coord(
        iris.coords.DimCoord(depth.data,
                             var_name='depth',
                             units='m',
                             attributes={'positive': 'down'}), 0)
    # add time coord
    year = raw_info['reference_year']
    time = Unit('months since ' + str(year) + '-01-01 00:00:00',
                calendar='gregorian')
    cube = iris.util.new_axis(cube)
    cube.add_dim_coord(
        iris.coords.DimCoord(6.,
                             standard_name='time',
                             units=time,
                             bounds=[0., 12.]), 0)

    fix_var_metadata(cube, var_info)
    fix_coords(cube)
    _fix_data(cube, var)
    set_global_atts(cube, attrs)
    save_variable(cube, var, out_dir, attrs, unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        in_files = collect_files(in_dir, var, cfg)
        logger.info("CMORizing var %s from input set %s", var, vals['file'])
        raw_info = cfg['variables'][var]
        raw_info.update({
            'var': var,
            'reference_year': cfg['custom']['reference_year'],
        })
        glob_attrs['mip'] = vals['mip']
        extract_variable(in_files, out_dir, glob_attrs, raw_info, cmor_table)
