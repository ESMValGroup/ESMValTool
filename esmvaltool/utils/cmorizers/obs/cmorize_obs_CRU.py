"""ESMValTool CMORizer for CRU data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.02/cruts.1811131722.v4.02/

Last access
    20190516

Download and processing instructions
    Download the following files:
        {raw_name}/cru_ts4.02.1901.2017.{raw_name}.dat.nc.gz
    where {raw_name} is the name of the desired variable(s).

"""

import gzip
import logging
import os
import shutil

import iris

import esmvaltool.utils.cmorizers.obs.utilities as utils

logger = logging.getLogger(__name__)


def _clean(filepath):
    """Remove unzipped input file."""
    if os.path.isfile(filepath):
        os.remove(filepath)
        logger.info("Removed cached file %s", filepath)


def _extract_variable(raw_var, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))
    utils.fix_var_metadata(cube, cmor_info)
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    if var in ('tas', ):
        utils.add_height2m(cube)
    utils.save_variable(cube,
                        var,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def _unzip(filepath, out_dir):
    """Unzip `*.gz` file."""
    filename = os.path.basename(filepath.replace('.gz', ''))
    new_path = os.path.join(out_dir, filename)
    with gzip.open(filepath, 'rb') as zip_file:
        with open(new_path, 'wb') as new_file:
            shutil.copyfileobj(zip_file, new_file)
    logger.info("Succefully extracted file to %s", new_path)
    return new_path


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    cfg = utils.read_cmor_config('CRU.yml')
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']
    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        raw_var = var_info.get('raw', var)
        zip_file = os.path.join(in_dir, raw_filepath.format(raw_name=raw_var))
        if not os.path.isfile(zip_file):
            logger.debug("Skipping '%s', file '%s' not found", var, zip_file)
            continue
        logger.info("Found input file '%s'", zip_file)
        filepath = _unzip(zip_file, out_dir)
        _extract_variable(raw_var, cmor_info, glob_attrs, filepath, out_dir)
        _clean(filepath)
