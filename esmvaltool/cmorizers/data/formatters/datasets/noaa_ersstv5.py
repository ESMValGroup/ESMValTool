"""ESMValTool CMORizer for NOAA ERSST data, version 5.

   This is the CMORizer script for the NOAA Extended Reconstructed Sea Surface
   Temperature (ERSST) data of version 5.

Tier
    Tier 2: open dataset.

Source
    https://doi.org/10.7289/V5T72FNM

Last access
    20200520

Download and processing instructions
    The data is available via an ftp-server provided by NOAA. The script
    below can be used to automate the process.

.. code-block:: bash

    #!/bin/env bash

    url=ftp://ftp.ncdc.noaa.gov/pub/data/cmb/ersst/v5/netcdf
    yr=1854
    while [ $yr -le 2019 ]
    do
        nm=1
        while [ $nm -le 12 ]
        do
            if [ $nm -lt 10 ]
            then
                nm=0$nm
            fi
            wget $url/ersst.v5.$yr$nm.nc
            nm=$(expr $nm + 1)
        done
        yr=$(expr $yr + 1)
    done

"""

import logging
import os, re

import iris
import cf_units

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

def _get_filepaths(in_dir, basename): 
    """Find correct name of file (extend basename with timestamp).
        return 2 lists as files differ from 2008"""
    # load all files in folder get sst into one cube, each have 1 time
    regex = re.compile(basename) 
    return_files = []
    return_files_gr08 = []
    for file in os.listdir(in_dir):

        if regex.match(file):
            yr = file.split('.')[2][:4]  # ersst.v5.$yr$nm.nc
            if int(yr)<2008:
                return_files.append(os.path.join(in_dir, file))
            else:
                return_files_gr08.append(os.path.join(in_dir, file))

    return return_files, return_files_gr08

def _fix_time_coord(cube, field, filename):
    """Set time points to central day of month."""
    time_coord = cube.coord('time')
    _unit = time_coord.units  # [0] slice for time_coord?
    new_time = [d.replace(day=15) for d in _unit.num2date(time_coord.points)]
    time_coord.points = _unit.date2num(new_time)
    time_coord.units = cf_units.Unit(time_coord.units.origin, calendar='standard')
    time_coord.long_name = 'Time'

def modify_cal(cube, field, filename):
    tcoord = cube.coord('time')
    tcoord.units = cf_units.Unit(tcoord.units.origin, calendar='standard')


def _extract_variable(raw_var, cmor_info, attrs, filepaths, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    # cube = iris.load_cube(filepath, raw_var)

    cubels = iris.load(filepaths, raw_var, _fix_time_coord)
    removed_attr = iris.util.equalise_attributes(cubels)
    iris.util.unify_time_units(cubels)
    cube = cubels.concatenate_cube()
    cube = iris.util.squeeze(cube)

    ## regrid - time metadata for before 2008, 

    utils.fix_var_metadata(cube, cmor_info)
    utils.fix_coords(cube)

    utils.set_global_atts(cube, attrs)
    utils.save_variable(cube,
                        var,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    filepaths = _get_filepaths(in_dir,cfg['filename'])
    # filepath = os.path.join(in_dir, cfg['filename'])

    if len(filepaths[0]) > 0 or len(filepaths[1]) > 0:
        totalfiles = len(filepaths[0]) + len(filepaths[1])
        logger.info("%d files before 2008", len(filepaths[0]))
        logger.info("Found %d input files in '%s'", totalfiles, in_dir)
    else:
        logger.info("No files found, basename: %s", cfg['filename'])

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        raw_var = var_info.get('raw', var)
        _extract_variable(raw_var, cmor_info, glob_attrs, filepaths[0], out_dir)
        _extract_variable(raw_var, cmor_info, glob_attrs, filepaths[1], out_dir)
