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
      - Tick variables LAI and FAPAR
      - Tick satellite SPOT (System Pour l'Observation de la Terre)
      - Tick sensor VGT (Vegetation)
      - Tick horizontal resolution 1km
      - Tick product version V1
      - Tick all available years
      - Tick all available months
      - Tick Nominal day 20
   - Click 'submit form'
   - According to ESMValTool practice, put them in the right rawobsdir folder

Notes
-----
   - This script regrids and cmorizes the above dataset.
   - Request might need to be split into chunks to not exceed download limit

Caveats
   - Fails setting standard name for variable FAPAR

Modification history
   20200512-crezee_bas: adapted to reflect changes in download form by CDS.
   20190703-crezee_bas: written.
"""

import glob
import logging
import os
from copy import deepcopy
from datetime import datetime
from warnings import catch_warnings, filterwarnings

import cf_units
import iris

from esmvalcore.preprocessor import regrid
from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _attrs_are_the_same(cubelist):
    # assume they are the same
    attrs_the_same = True
    allattrs = cubelist[0].attributes
    for key in allattrs:
        try:
            unique_attr_vals = {cube.attributes[key] for cube in cubelist}
        # This exception is needed for valid_range, which is an
        # array and therefore not hashable
        except TypeError:
            unique_attr_vals = {tuple(cube.attributes[key])
                                for cube in cubelist}
        if len(unique_attr_vals) > 1:
            attrs_the_same = False
            print("Different values found for {0}-attribute: {1}"
                  .format(key, unique_attr_vals))
    return attrs_the_same


def _cmorize_dataset(in_file, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input file '%s'",
                var['short_name'], in_file)
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']

    cmor_table = cfg['cmor_table']
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    cube = iris.load_cube(
        str(in_file),
        constraint=utils.var_name_constraint(var['raw']))

    # Set correct names
    cube.var_name = definition.short_name
    if definition.standard_name:
        cube.standard_name = definition.standard_name

    cube.long_name = definition.long_name

    # Convert units if required
    cube.convert_units(definition.units)

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    logger.info("Saving CMORized cube for variable %s", cube.var_name)
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
        _, infile_tail = os.path.split(infile)
        outfile_tail = infile_tail.replace('c3s', 'c3s_regridded')
        outfile = os.path.join(cfg['work_dir'], outfile_tail)
        with catch_warnings():
            filterwarnings(
                action='ignore',
                # Full message:
                # UserWarning: Skipping global attribute 'long_name':
                #              'long_name' is not a permitted attribute
                message="Skipping global attribute 'long_name'",
                category=UserWarning,
                module='iris',
            )
            lai_cube = iris.load_cube(infile,
                                      constraint=utils.var_name_constraint(
                                          var['raw']))
        lai_cube = regrid(lai_cube, cfg['custom']['regrid_resolution'],
                          'nearest')
        logger.info("Saving: %s", outfile)

        iris.save(lai_cube, outfile)


def _set_time_bnds(in_dir, var):
    """Set time_bnds by using attribute and returns a cubelist."""
    # This is a complicated expression, but necessary to keep local
    # variables below the limit, otherwise prospector complains.
    cubelist = iris.load(glob.glob(os.path.join(
        in_dir, var['file'].replace('c3s', 'c3s_regridded'))))

    # The purpose of the following loop is to remove any attributes
    # that differ between cubes (otherwise concatenation over time fails).
    # In addition, care is taken of the time coordinate, by adding the
    # time_coverage attributes as time_bnds to the time coordinate.
    for n_cube, _ in enumerate(cubelist):
        time_coverage_start = cubelist[n_cube].\
            attributes.pop('time_coverage_start')
        time_coverage_end = cubelist[n_cube].\
            attributes.pop('time_coverage_end')

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

    return cubelist


def cmorization(in_dir, out_dir, cfg, cfg_user):
    """Cmorization func call."""
    # run the cmorization
    # Pass on the workdir to the cfg dictionary
    cfg['work_dir'] = cfg_user['work_dir']
    # If it doesn't exist, create it
    if not os.path.isdir(cfg['work_dir']):
        logger.info("Creating working directory for regridding: %s",
                    cfg['work_dir'])
        os.mkdir(cfg['work_dir'])

    for short_name, var in cfg['variables'].items():
        var['short_name'] = short_name
        logger.info("Processing var %s", short_name)

        # Regridding
        logger.info("Start regridding to: %s",
                    cfg['custom']['regrid_resolution'])
        _regrid_dataset(in_dir, var, cfg)
        logger.info("Finished regridding")

        # File concatenation
        logger.info("Start setting time_bnds")
        cubelist = _set_time_bnds(cfg['work_dir'], var)

        # Loop over two different platform names
        for platformname in ['SPOT-4', 'SPOT-5']:
            # Now split the cubelist on the different platform
            logger.info("Start processing part of dataset: %s", platformname)
            cubelist_platform = cubelist.extract(iris.AttributeConstraint(
                platform=platformname))
            for n_cube, _ in enumerate(cubelist_platform):
                cubelist_platform[n_cube].attributes.pop('identifier')
            if cubelist_platform:
                assert _attrs_are_the_same(cubelist_platform)
                cube = cubelist_platform.concatenate_cube()
            else:
                logger.warning("No files found for platform %s \
                               (check input data)", platformname)
                continue
            savename = os.path.join(cfg['work_dir'],
                                    var['short_name'] + platformname + '.nc')
            logger.info("Saving as: %s", savename)
            iris.save(cube, savename)
            logger.info("Finished file concatenation over time")
            in_file = savename
            logger.info("Start CMORization of file %s", in_file)
            _cmorize_dataset(in_file, var, cfg, out_dir)
            logger.info("Finished regridding and CMORizing %s", in_file)
