"""ESMValTool CMORizer for CRU data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.04/cruts.2004151855.v4.04/

Last access
    20190516

Download and processing instructions
    Download the following files:
        {raw_name}/cru_ts4.04.1901.2019.{raw_name}.dat.nc.gz
    where {raw_name} is the name of the desired variable(s).

Two files are generated per variable, one with version (e.g. TS4.04),
one with version + _stn1 (e.g. TS4.04_stn1), which is constrained on holding
gridpoint values relying on data from at least one station (i.e. removing
gridpoints solely relying on climatological infilling).

"""

import gzip
import logging
import os
import shutil
import numpy as np
from cf_units import Unit

import iris

from . import utilities as utils

logger = logging.getLogger(__name__)

import IPython
from traitlets.config import get_config
c = get_config()
c.InteractiveShellEmbed.colors = "Linux"
# IPython.embed(config=c)


def _clean(filepath):
    """Remove unzipped input file."""
    if os.path.isfile(filepath):
        os.remove(filepath)
        logger.info("Removed cached file %s", filepath)


def _extract_variable(short_name, var, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    # Fix units
    if 'raw_units' in var:
        cube.units = var['raw_units']
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.convert_units(cmor_info.units)
    # fix time units
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))

    IPython.embed(config=c)

    # Fix coordinates
    utils.fix_coords(cube)
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])

    # build contraints cube on stn < 1
    constraint_var = var.get('constraint', short_name)
    constr_cube = iris.load_cube(filepath,
                                 utils.var_name_constraint(constraint_var))
    utils.fix_coords(constr_cube)
    # fix time units
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))

    cube.data = np.ma.masked_where(constr_cube.data < 1., cube.data)

    # Save variable
    attrs = cfg['attributes_constraint']
    attrs['mip'] = var['mip']

    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def _unzip(short_name, var, raw_filepath, out_dir):
    """Unzip `*.gz` file."""
    raw_var = var.get('raw', short_name)
    zip_path = raw_filepath.format(raw_name=raw_var)
    if not os.path.isfile(zip_path):
        logger.debug("Skipping '%s', file '%s' not found", short_name,
                     zip_path)
        return None
    logger.info("Found input file '%s'", zip_path)
    filename = os.path.basename(zip_path.replace('.gz', ''))
    new_path = os.path.join(out_dir, filename)
    with gzip.open(zip_path, 'rb') as zip_file:
        with open(new_path, 'wb') as new_file:
            shutil.copyfileobj(zip_file, new_file)
    logger.info("Succefully extracted file to %s", new_path)
    return new_path


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        filepath = _unzip(short_name, var, raw_filepath, out_dir)
        if filepath is None:
            continue
        _extract_variable(short_name, var, cfg, filepath, out_dir)
        _clean(filepath)
