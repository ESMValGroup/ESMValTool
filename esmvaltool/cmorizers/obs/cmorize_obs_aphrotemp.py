"""ESMValTool CMORizer for CRU data.

Tier
    Tier 2: other freely-available dataset.

Source
    http://search.diasjp.net/en/dataset/AphroTemp
    http://aphrodite.st.hirosaki-u.ac.jp/index.html

Last access
    20200225

Download and processing instructions
    Download the following files:
        APHRO_MA_TAVE_{version}.1961-2007.nc.tgz
        version for APHRO_MA_TAVE: [025deg_V1204R1, 050deg_V1204R1]

Refs:
    Yasutomi, N., Hamada, A., Yatagai, A. (2011) Development of a long-term
    daily gridded temperature dataset and its application to rain/snow
    discrimination of daily precipitation,
    Global Environmental Research 15 (2), 165-172

Issues:
    Data was downloaded using the from dias provided python scripts for above
    mentioned tar balls (download_APHRO.py). The single files for 2007 in the
    tar were empty. The single files for 2007 have to be downloaded separately.
    Handling of the tgz archives is not implemented.
"""

import logging
import iris

from warnings import catch_warnings, filterwarnings
from pathlib import Path

from esmvalcore.preprocessor import monthly_statistics
from . import utilities as utils

logger = logging.getLogger(__name__)


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
    utils._set_units(cube, var.get('raw_units', short_name))
    cube.convert_units(definition.units)
    utils.fix_var_metadata(cube, definition)

    # fix coordinates
    utils.fix_coords(cube)
    utils.add_scalar_height_coord(cube, 2.)

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
                        unlimited_dimensions=['time']
                        )

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
                                unlimited_dimensions=['time']
                                )


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filename = cfg['filename']

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        for version in cfg['attributes']['version'].values():
            logger.info("CMORizing variable '%s'", short_name)

            filenames = raw_filename.format(version=version)
            for filepath in sorted(Path(in_dir).glob(filenames)):
                _extract_variable(short_name, var, cfg, filepath, out_dir,
                                  version)
