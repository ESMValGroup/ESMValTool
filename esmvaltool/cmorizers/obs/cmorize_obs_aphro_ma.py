"""ESMValTool CMORizer for CRU data.

Tier
    Tier 2: other freely-available dataset.

Source
    http://search.diasjp.net/en/dataset/APHRO_MA
    http://aphrodite.st.hirosaki-u.ac.jp/index.html

Last access
    20200225

Download and processing instructions
    Download the following files:
        APHRO_MA_{version}.1951-2007.nc.gz.tar
        version for APHRO_MA: [025deg_V1101, 050deg_V1101]

Refs:
    Yatagai, A., K. Kamiguchi, O. Arakawa, A. Hamada, N. Yasutomi, and
    A. Kitoh, 2012: APHRODITE: Constructing a Long-Term Daily Gridded
    Precipitation Dataset for Asia Based on a Dense Network of Rain Gauges.
    Bull. Amer. Meteor. Soc., 93, 1401–1415
    https://doi.org/10.1175/BAMS-D-11-00122.1

Issues:
    I downloaded the data using the from dias provided python scripts for above
    mentioned tar balls.
"""

import gzip
import logging
import shutil
import os
from warnings import catch_warnings, filterwarnings
from pathlib import Path

import iris

from esmvalcore.preprocessor import monthly_statistics

from . import utilities as utils

logger = logging.getLogger(__name__)


def _clean(filepath):
    """Remove unzipped input file."""
    if os.path.isfile(filepath):
        os.remove(filepath)
        logger.info("Removed cached file %s", filepath)


def _extract_variable(short_name, var, cfg, filepath, out_dir, version):
    """Extract variable."""

    logger.info("CMORizing variable '%s' from input file '%s'",
                short_name, filepath)

    with catch_warnings():
        filterwarnings(
            action='ignore',
            message="Skipping global attribute 'calendar': 'calendar' is .*",
            category=UserWarning,
            module='iris',
        )
        cube = iris.load_cube(
            str(filepath),
            constraint=utils.var_name_constraint(var['raw']),
            )

    # Fix var units
    definition = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.units = var.get('raw_units', short_name)
    cube.convert_units(definition.units)
    utils.fix_var_metadata(cube, definition)

    # fix coordinates
    utils.fix_coords(cube)

    # Fix metadata
    attrs = cfg['attributes'].copy()
    attrs['mip'] = var['mip']
    attrs['version'] = version
    attrs['source'] = attrs['source']
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])

    if 'add_mon' in var.keys():
        if var['add_mon']:
            logger.info("Building monthly means")

            # Calc monthly
            cube = monthly_statistics(cube)
            cube.remove_coord('month_number')
            cube.remove_coord('year')

            # Fix metadata
            attrs['mip'] = 'Amon'

            # Fix coordinates
            utils.fix_coords(cube)

            # Save variable
            utils.save_variable(cube,
                                short_name,
                                out_dir,
                                attrs,
                                unlimited_dimensions=['time'])


def _unzip(short_name, zip_path, out_dir):
    """Unzip `*.gz` file."""
    if not os.path.isfile(zip_path):
        logger.debug("Skipping '%s', file '%s' not found", short_name,
                     zip_path)
        return None
    logger.info("Found input file '%s'", zip_path)
    filename = os.path.basename(zip_path).replace('.gz', '')
    new_path = os.path.join(out_dir, filename)
    with gzip.open(zip_path, 'rb') as zip_file:
        with open(new_path, 'wb') as new_file:
            shutil.copyfileobj(zip_file, new_file)
    logger.info("Succefully extracted file to %s", new_path)
    return new_path


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filename = cfg['filename']

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        for version in cfg['attributes']['version'].values():
            logger.info("CMORizing variable '%s'", short_name)

            filenames = raw_filename.format(version=version)
            for zip_path in sorted(Path(in_dir).glob(filenames)):
                filepath = _unzip(short_name, zip_path, out_dir)
                if filepath is None:
                    continue
                _extract_variable(short_name, var, cfg, filepath, out_dir,
                                  version)
                _clean(filepath)
