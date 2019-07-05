"""ESMValTool CMORizer for cds-satellite-lai-fapar data.

Tier
   Tier 3
Source
   https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-lai-fapar?tab=form
Last access
   20190703

Download and processing instructions
   - Open in a browser the data source as specified above
   - Put the right ticks:
      - Tick variables FAPAR and LAI
      - Tick SPOT
      - Tick V1
      - Tick all available years
      - Tick all available months
      - Tick Nominal day 20
   - Click 'submit form'
   - According to ESMValTool practice, put them in the right rawobsdir folder

Notes
-----
   - This script regrids and cmorizes the above dataset.

Caveats
   - Handling of fAPAR not correct yet (standard_name can not be set)
   - This dataset has custom time bounds (see dataset documentation
     available from source)
   - Currently files are written to rawobsdir and not cleaned-up

Modification history
   20190703-A_crez_ba: written.
"""

import glob
import logging
import os
from copy import deepcopy
from datetime import datetime
from warnings import catch_warnings, filterwarnings

import cf_units
import iris

from esmvalcore.cmor.table import CMOR_TABLES
from esmvalcore.preprocessor import regrid
from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _cmorize_dataset(in_file, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input file '%s'",
                var['short_name'], in_file)
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']

    cmor_table = cfg['cmor_table']
    definition = cmor_table.get_variable(var['mip'], var['short_name'])
    if not definition:
        cmor_table = CMOR_TABLES['custom']
        definition = cmor_table.get_variable(var['mip'], var['short_name'])

    with catch_warnings():
        filterwarnings(
            action='ignore',
            message="Ignoring netCDF variable 'tcc' invalid units '(0 - 1)'",
            category=UserWarning,
            module='iris',
        )
        cube = iris.load_cube(
            str(in_file),
            constraint=utils.var_name_constraint(var['raw']),
        )

    # Set correct names
    cube.var_name = definition.short_name
    try:
        cube.standard_name = definition.standard_name
    except ValueError:
        logger.warning(
            "Could not set standard name %s for variable (var_name=) %s",
            definition.standard_name, cube.var_name)
    cube.long_name = definition.long_name

    # Convert units if required
    cube.convert_units(definition.units)

    # Make latitude increasing
    cube = cube[:, ::-1, ...]

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    logger.info("Saving cube\n%s", cube)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)

    return in_file


def _regrid_dataset(in_dir, var, cfg):
    """
    Regridding of original files.

    This function regrids each file and write to disk appending 'regrid'
    in front of filename.
    """
    filelist = glob.glob(os.path.join(in_dir, var['file']))
    for infile in filelist:
        infile_head, infile_tail = os.path.split(infile)
        outfile_tail = infile_tail.replace('c3s', 'c3s_regridded')
        outfile = os.path.join(infile_head, outfile_tail)
        lai_cube = iris.load_cube(infile,
                                  constraint=utils.var_name_constraint(
                                      var['raw']))
        lai_cube = regrid(lai_cube, cfg['custom']['regrid_resolution'],
                          'nearest')
        logger.info("Saving: %s", outfile)
        iris.save(lai_cube, outfile)


def _concatenate_dataset_over_time(in_dir, var):
    """Concatenate single files over time and returns on single cube."""
    # This is a complicated expression, but necessary to keep local
    # variables below the limit, otherwise prospector complains.
    cubelist = iris.load(glob.glob(os.path.join(
        in_dir, var['file'].replace('c3s', 'c3s_regridded'))))

    # For saving the identifiers
    identifiers = []
    # The purpose of the following loop is to remove any attributes
    # that differ between cubes (otherwise concatenation over time fails).
    # In addition, care is taken of the time coordinate, by adding the
    # time_coverage attributes as time_bnds to the time coordinate.
    for n_cube, _ in enumerate(cubelist):
        time_coverage_start = cubelist[n_cube].\
            attributes.pop('time_coverage_start')
        time_coverage_end = cubelist[n_cube].\
            attributes.pop('time_coverage_end')
        # Remove identifier from attrs, add to list for later usage
        identifiers.append(cubelist[n_cube].attributes.pop('identifier'))

        # Now put time_coverage_start/end as time_bnds
        # Convert time_coverage_xxxx to datetime
        bnd_a = datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%SZ")
        bnd_b = datetime.strptime(time_coverage_end, "%Y-%m-%dT%H:%M:%SZ")

        # Put in shape for time_bnds
        time_bnds_datetime = [bnd_a, bnd_b]

        # Read dataset time unit and calendar from file
        dataset_time_unit = str(cubelist[n_cube].coord('time').units)
        dataset_time_calender = cubelist[n_cube].coord('time').units.calendar
        # Convert datetime
        time_bnds = cf_units.date2num(time_bnds_datetime, dataset_time_unit,
                                      dataset_time_calender)
        # Put them on the file
        cubelist[n_cube].coord('time').bounds = time_bnds

    # Now the cubes can be concatenated over the time dimension
    cube = cubelist.concatenate_cube()

    # Add identifiers from each cube to the concatenated cube
    # as a comma separated list
    cube.attributes['identifiers_comma_separated'] = ','.join(identifiers)
    return cube


def cmorization(in_dir, out_dir, cfg):
    """Cmorization func call."""
    # run the cmorization
    for short_name, var in cfg['variables'].items():
        var['short_name'] = short_name
        # Now call every function with the following pattern (in_dir, var, cfg)
        logger.info("Start regridding to: %s",
                    cfg['custom']['regrid_resolution'])
        _regrid_dataset(in_dir, var, cfg)
        logger.info("Finished regridding")

        # First collect all information
        inpfile = os.path.join(in_dir, var['file'])
        logger.info("CMORizing var %s from file %s", short_name, inpfile)

        logger.info("Start file concatenation over time")
        result_cube = _concatenate_dataset_over_time(
            in_dir, var)
        # Write it to disk
        savename = os.path.join(in_dir, var['short_name'] + '_regridded.nc')
        logger.info("saving as: %s", savename)
        iris.save(result_cube, savename)
        logger.info("Finished file concatenation over time")

        with catch_warnings():
            filterwarnings(
                action='ignore',
                message=('WARNING: missing_value not used since it\n'
                         'cannot be safely cast to variable data type'),
                category=UserWarning,
                module='iris',
            )
            in_file = savename
            logger.info("Start CMORization of file")
            _cmorize_dataset(in_file, var, cfg, out_dir)
        logger.info("Finished regridding and CMORizing %s", in_file)
