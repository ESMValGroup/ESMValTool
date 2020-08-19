"""ESMValTool CMORizer for APHRODITE Monsoon Asia (APHRO-MA) data.

Tier
    Tier 3: restricted dataset.

Source
    http://aphrodite.st.hirosaki-u.ac.jp/download/

Last access
    20200306

Download and processing instructions
    Register at
    http://aphrodite.st.hirosaki-u.ac.jp/download/create/

    Download the following files from
    http://aphrodite.st.hirosaki-u.ac.jp/product/:
        APHRO_V1808_TEMP/APHRO_MA
            025deg_nc/APHRO_MA_TAVE_025deg_V1808.nc.tgz
            050deg_nc/APHRO_MA_TAVE_050deg_V1808.nc.tgz
        APHRO_V1101/APHRO_MA
            025deg_nc/APHRO_MA_025deg_V1101.1951-2007.nc.gz.tar
            050deg_nc/APHRO_MA_050deg_V1101.1951-2007.nc.gz.tar
        APHRO_V1101EX_R1/APHRO_MA
            025deg_nc/APHRO_MA_025deg_V1101_EXR1.nc.tgz
            050deg_nc/APHRO_MA_050deg_V1101_EXR1.nc.tgz

    Please untar / unzip all *.tar *.tgz *.gz files in the same directory
    (no subdirectories!) prior to running the cmorizer!

Issues:
    In input file APHRO_MA_TAVE_050deg_V1808.2015.nc the input variable is
    called ta instead of tave as in the other files.
    Currently resolved using raw_fallback: ta in case of thrown
    iris.exceptions.ConstraintMismatchError

Refs:
    APHRO_V1101 and APHRO_V1101EX_R1
    Yatagai, A., K. Kamiguchi, O. Arakawa, A. Hamada, N. Yasutomi, and
    A. Kitoh, 2012: APHRODITE: Constructing a Long-Term Daily Gridded
    Precipitation Dataset for Asia Based on a Dense Network of Rain Gauges.
    Bull. Amer. Meteor. Soc., 93, 1401–1415
    https://doi.org/10.1175/BAMS-D-11-00122.1

    APHRO_V1808_TEMP
    Yasutomi, N., Hamada, A., Yatagai, A. (2011) Development of a long-term
    daily gridded temperature dataset and its application to rain/snow
    discrimination of daily precipitation,
    Global Environmental Research 15 (2), 165-172
"""

import logging
from warnings import catch_warnings, filterwarnings
from pathlib import Path

import iris

from esmvalcore.preprocessor import monthly_statistics

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, cfg, filepath, out_dir, version):
    """Extract variable."""
    logger.info("CMORizing variable '%s' from input file '%s'", short_name,
                filepath)

    with catch_warnings():
        filterwarnings(
            action='ignore',
            message="Skipping global attribute 'calendar': 'calendar' is .*",
            category=UserWarning,
            module='iris',
        )
        try:
            cube = iris.load_cube(
                str(filepath),
                constraint=utils.var_name_constraint(var['raw']),
            )
        except iris.exceptions.ConstraintMismatchError:
            cube = iris.load_cube(
                str(filepath),
                constraint=utils.var_name_constraint(var['raw_fallback']),
            )

    # Fix var units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.units = var.get('raw_units', short_name)
    cube.convert_units(cmor_info.units)
    utils.fix_var_metadata(cube, cmor_info)

    # fix coordinates
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)
    utils.fix_coords(cube)

    # Fix metadata
    attrs = cfg['attributes'].copy()
    attrs['mip'] = var['mip']
    attrs['version'] = version.replace('_', '-')
    attrs['reference'] = var['reference']
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


def cmorization(in_dir, out_dir, cfg, _, __, ___):
    """Cmorization func call."""
    raw_filename = cfg['filename']

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        for version in var['version'].values():
            logger.info("CMORizing variable '%s'", short_name)
            filenames = raw_filename.format(raw_file_var=var['raw_file_var'],
                                            version=version)
            for filepath in sorted(Path(in_dir).glob(filenames)):
                _extract_variable(short_name, var, cfg, filepath, out_dir,
                                  version)
