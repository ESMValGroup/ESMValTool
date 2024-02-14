"""ESMValTool CMORizer for ESACCI-SNOW data.

Tier
   Tier 2: other freely-available dataset.

Source
   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/snow/data

Last access
   20240214

Download and processing instructions
   Download the data from:
        scfg/AVHRR-MERGED/v2.0/
        swe/MERGED/v2.0/
      Put all scfg files in a single directory called 'scfg'
      (no subdirectories with years or months).
      Put all swe files in a single directory called 'swe'
      (no subdirectories with years or months).
"""

import glob
import logging
import os
from copy import deepcopy
from datetime import datetime
from dateutil import relativedelta

import iris
from iris.util import equalise_attributes, unify_time_units
from esmvalcore.cmor.table import CMOR_TABLES

from esmvaltool.cmorizers.data.utilities import (
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def _extract_variable(in_files, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input files '%s'",
                var['short_name'], ', '.join(in_files))
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    attributes['raw'] = var['raw']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    # load all input files (1 year) into 1 cube
    cube_list = iris.load_raw(in_files)
    drop_attrs = [
        'source', 'date_created', 'history', 'tracking_id',
        'id', 'time_coverage_start', 'time_coverage_end', 'platform'
    ]

    for cube in cube_list:
        for attr in drop_attrs:
            cube.attributes.pop(attr)
        cube.coord('time').points = cube.coord('time').core_points().astype(
            'float64')
    iris.util.unify_time_units(cube_list)
    cube = cube_list.concatenate_cube()

#    # keep the following raw cube attributes
#    attrs_to_keep = [
#        "institution", "Institution",
#        "institute_id", "VersionID",
#        "experiment_id",
#        "source", "Source",  # overrides empty string default
#        "model_id", "ModelID",
#        "contact", "Contact",
#        "references",
#        "tracking_id",
#        "mip_specs",  # described by "mip" already
#        "source_id", "SourceID",
#        "product", "Product",
#        "frequency", "Frequency",
#        "creation_date",
#        "project_id", "ProjectID",
#        "table_id", "TableID",
#        "title", "Title",
#        "modeling_realm",
#        "doi",
#        "VersionID",  # described by "version" already
#    ]
#
#    attrs_to_keep_exist = [
#        att for att in cube.attributes if att in attrs_to_keep
#    ]
#    for att in attrs_to_keep_exist:
#        attributes[att] = cube.attributes[att]
#
#    utils.set_global_atts(cube, attributes)

    # Set correct names
    cube.var_name = definition.short_name
    # cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name




    # Fix units (if needed)
    # input variable reports m-3 m-3 instead of m3 m-3
    if cube.var_name == "sm":
        cube.units = definition.units
    # Convert units to CMOR units
    cube.convert_units(definition.units)

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    # Roll longitude
    cube.coord('longitude').points = cube.coord('longitude').points + 180.
    nlon = len(cube.coord('longitude').points)
    cube.data = da.roll(cube.core_data(), int(nlon / 2), axis=-1)

    # Fix coordinates
    cube = _fix_coordinates(cube, definition)

    cube.coord('latitude').attributes = None
    cube.coord('longitude').attributes = None

    cube = _fix_time_monthly(cube)

    logger.debug("Saving cube\n%s", cube)
    logger.debug("Setting time dimension to UNLIMITED while saving!")
    utils.save_variable(cube, cube.var_name,
                        out_dir, attributes,
                        unlimited_dimensions=['time'])
    logger.info("Finished CMORizing %s", ', '.join(in_files))

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize ESACCI-SNOW dataset."""
    glob_attrs = cfg['attributes']

    logger.info("Starting cmorization for tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    logger.info('Cmorizing ESACCI-SNOW version %s', glob_attrs['version'])

    if start_date is None:
        start_date = datetime(1985,1,1) #1979, 1, 1)
    if end_date is None:
        end_date = datetime(1985,1,1) #2019, 12, 31)

    for short_name, var in cfg['variables'].items():
        if 'short_name' not in var:
            var['short_name'] = short_name
        loop_date = start_date
        while loop_date <= end_date:
            filepattern = os.path.join(in_dir, var['file'].format(year=loop_date.year))
            in_files = glob.glob(filepattern)
            if not in_files:
                logger.info(f'{loop_date.year}: no data not found for variable {short_name}')
            else:
                _extract_variable(in_files, var, cfg, out_dir)

            loop_date += relativedelta.relativedelta(years=1)
