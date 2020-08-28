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
import os

import iris

from . import utilities as utils

logger = logging.getLogger(__name__)


def _fix_time_coord(cube):
    """Set time points to central day of month."""
    time_coord = cube.coord('time')
    _unit = time_coord[0].units
    new_time = [d.replace(day=15) for d in _unit.num2date(time_coord.points)]
    time_coord.points = _unit.date2num(new_time)


def _extract_variable(raw_var, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))
    cube = iris.util.squeeze(cube)
    _fix_time_coord(cube)
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)
    utils.save_variable(cube,
                        var,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']
    filepath = os.path.join(in_dir, cfg['filename'])
    logger.info("Found input file '%s'", filepath)

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        raw_var = var_info.get('raw', var)
        _extract_variable(raw_var, cmor_info, glob_attrs, filepath, out_dir)
