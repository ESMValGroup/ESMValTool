"""
ESMValTool CMORizer for BerkeleyEarth data.

Tier
    Tier 2: other freely-available dataset.

Source
    http://berkeleyearth.org/data/
    Monthly Land + Ocean
    Average Temperature with Air Temperatures at Sea Ice
    (Recommended; 1850 – Recent)
    1º x 1º Latitude-Longitude Grid (~400 MB)

Last access
    20200225

Download and processing instructions
    Download the following file:
    http://berkeleyearth.lbl.gov/auto/Global/Gridded/Land_and_Ocean_LatLong1.nc

"""

import logging
import os
import re
from warnings import catch_warnings, filterwarnings

import numpy as np
import cf_units
import iris
from iris import coord_categorisation

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def reinit_broken_time(cube_anom, cube_clim, climstart, climend):
    """
    Fix broken time.

    The time coordinates are a big mess (given as floats in years A.D.)
    best to reinitialize them from scratch
    """
    logger.info("Reinitializing broken time coordinate")
    time_raw = cube_anom.coord('time')

    n_years, n_add_mon = len(time_raw.points) // 12, len(time_raw.points) % 12
    start_year = int(time_raw.points[0])
    n_days = (n_years + n_add_mon / 12) * 365.25 + 50  # have some extra length
    climcenter = (climend - climstart) // 2

    times = iris.coords.DimCoord(np.arange(int(n_days), dtype=float),
                                 var_name='time',
                                 standard_name='time',
                                 long_name='time',
                                 units=cf_units.Unit(
                                     'days since {}-01-01 00:00:00'.format(
                                         start_year),
                                     calendar=cf_units.CALENDAR_STANDARD))

    # init a dummy cube to enable coord_categorisation
    dummycube = iris.cube.Cube(np.zeros(int(n_days), np.int),
                               dim_coords_and_dims=[(times, 0)])
    coord_categorisation.add_year(dummycube, 'time', name='year')
    coord_categorisation.add_month_number(dummycube, 'time', name='month')

    # build timecoord for the anomaly cube
    dummycube = dummycube.aggregated_by(['year', 'month'], iris.analysis.MEAN)
    dummycube = dummycube[:(n_years * 12 + n_add_mon)]
    timecoord_anom = dummycube.coord('time')

    # build timecoord for the climatology cube
    dummycube_clim = dummycube.extract(iris.Constraint(
        year=lambda cell: cell == climstart + climcenter))
    timecoord_clim = dummycube_clim.coord('time')

    # change to the new time coordinates
    cube_anom.remove_coord('time')
    cube_anom.add_dim_coord(timecoord_anom, 0)
    cube_clim.add_dim_coord(timecoord_clim, 0)

    # convert time units to standard
    utils.convert_timeunits(cube_anom, 1950)
    utils.convert_timeunits(cube_clim, 1950)

    return (cube_anom, cube_clim)


def calc_abs_temperature(cube_anom, cube_clim, short_name):
    """Derive absolute tas values."""
    logger.info("Deriving absolute temperature fields")

    # prepare cubes
    for cube in [cube_anom, cube_clim]:
        cube.attributes.pop('valid_max')
        cube.attributes.pop('valid_min')

    # declare attributes
    fill_value = cube_clim.data.fill_value
    dtype = cube_anom.dtype
    var_name = short_name
    units = cube_clim.units

    # init data array
    crds_n_dims = [(cor.copy(), i) for i, cor in enumerate(cube_anom.coords())]
    shape = [x[0].shape[0] for x in crds_n_dims]
    array = np.ma.ones(shape, dtype=dtype) * fill_value
    array.mask = True
    array.fill_value = fill_value

    # calculate abs fields
    for i in range(0, cube_anom.coord('time').shape[0]):
        with catch_warnings():
            filterwarnings(
                action='ignore',
                message='.* not used since it\ncannot be safely cast to'
                        ' variable data type *',
                category=UserWarning,
                module='iris',
            )
            array[i] = cube_clim[i % 12].data + cube_anom[i].data

    # build absolute tas cube
    cube_abs = iris.cube.Cube(array, var_name=var_name,
                              units=units,
                              dim_coords_and_dims=crds_n_dims)

    return cube_abs


def _extr_var_n_calc_abs_tas(short_name, var, cfg, filepath, out_dir):
    """Extract variable."""
    # load tas anomaly, climatology and sftlf
    with catch_warnings():
        filterwarnings(
            action='ignore',
            message='Ignoring netCDF variable .* invalid units .*',
            category=UserWarning,
            module='iris',
        )
        filterwarnings(
            action='ignore',
            message='.* not used since it\ncannot be safely cast to variable'
                    ' data type *',
            category=UserWarning,
            module='iris',
        )
        cubes = iris.load(filepath)

    # tas anomaly
    raw_var = var.get('raw', short_name)
    cube_anom = cubes.extract(utils.var_name_constraint(raw_var))[0]

    # tas climatology
    raw_var_clim = var.get('rawclim', short_name)
    cube_clim = cubes.extract(utils.var_name_constraint(raw_var_clim))[0]
    # information on time for the climatology are only present in the long_name
    climstart, climend = [int(x) for x in re.findall(r"\d{4}",
                                                     cube_clim.long_name)]

    # redo the broken time coordinate
    cube_anom, cube_clim = reinit_broken_time(cube_anom, cube_clim,
                                              climstart, climend)

    # derive absolute tas values
    cube_abs = calc_abs_temperature(cube_anom, cube_clim, short_name)

    # fix coordinates
    logger.info("Fixing coordinates")
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    short_names = [short_name, var['short_anom']]
    for s_name, cube in zip(short_names, [cube_abs, cube_anom]):
        cmor_info = cfg['cmor_table'].get_variable(var['mip'], s_name)

        utils.fix_coords(cube)
        if 'height2m' in cmor_info.dimensions:
            utils.add_height2m(cube)

        cube.units = var['raw_units']
        if s_name != 'tasa':
            cube.convert_units(cmor_info.units)

        utils.fix_var_metadata(cube, cmor_info)

    # save temperature data
    logger.info("Saving temperature data")
    comments = {'tas': "Temperature time-series calculated from the anomaly "
                       "time-series by adding the temperature climatology "
                       "for {}-{}".format(climstart, climend),
                'tasa': "Temperature anomaly with respect to the period"
                        " {}-{}".format(climstart, climend)}

    for s_name, cube in zip(short_names, [cube_abs, cube_anom]):
        attrs['comment'] = comments[s_name]
        utils.set_global_atts(cube, attrs)
        utils.save_variable(cube,
                            s_name,
                            out_dir,
                            attrs,
                            unlimited_dimensions=['time'])

    # sftlf
    # extract sftlf
    raw_var_sftlf = var.get('rawsftlf', short_name)
    cube_sftlf = cubes.extract(utils.var_name_constraint(raw_var_sftlf))[0]

    # fix coordinates
    utils.fix_coords(cube_sftlf)

    # cmorize sftlf units
    cmor_info_sftlf = cfg['cmor_table'].get_variable(var['rawsftlf_mip'],
                                                     var['rawsftlf_varname'])
    attrs_sftlf = cfg['attributes']
    attrs_sftlf['mip'] = var['rawsftlf_mip']
    if 'rawsftlf_units' in var:
        if 'rawsftlf_units' in var:
            cube_sftlf.units = var['rawsftlf_units']
        cube_sftlf.convert_units(cmor_info_sftlf.units)

    # fix metadata and save
    logger.info("Saving sftlf")
    utils.fix_var_metadata(cube_sftlf, cmor_info_sftlf)
    utils.set_global_atts(cube_sftlf, attrs_sftlf)
    utils.save_variable(cube_sftlf,
                        var['rawsftlf_varname'],
                        out_dir,
                        attrs_sftlf)


def cmorization(in_dir, out_dir, cfg, _, __, ___):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        _extr_var_n_calc_abs_tas(short_name, var, cfg, raw_filepath, out_dir)
