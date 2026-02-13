"""ESMValTool CMORizer for Sea Ice Concentration CDR.

Tier
   Tier 3: restricted dataset.

Source
   https://nsidc.org/data/g02202/versions/4

Last access
   20231213

Download and processing instructions
    Download data from:
    https://noaadata.apps.nsidc.org/NOAA/G02202_V4/south/monthly
    lat and lon from:
    https://noaadata.apps.nsidc.org/NOAA/G02202_V4/ancillary/
    area file:
    ftp://sidads.colorado.edu/DATASETS/seaice/polar-stereo/tools/
    pss25area_v3.dat

    https://nsidc.org/sites/default/files/g02202-v004-userguide_1_1.pdf

"""

import logging
import os
import re

import iris
import numpy as np
from cf_units import Unit
from esmvalcore.cmor._fixes.common import OceanFixGrid
from esmvalcore.cmor.fixes import get_time_bounds
from iris.coords import AuxCoord

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_filepaths(in_dir, basename, yyyy):
    """Find correct name of file (extend basename with timestamp)."""
    f_name = basename.format(year=yyyy)
    regex = re.compile(f_name)
    return_files = []
    for files in os.listdir(in_dir):
        if regex.match(files):
            return_files.append(os.path.join(in_dir, files))

    return return_files


def _fix_time_coord(cube, _field, _filename):
    """Set time points to central day of month."""
    time_coord = cube.coord("time")
    new_unit = Unit("days since 1850-01-01 00:00:00", calendar="standard")
    time_coord.convert_units(new_unit)
    old_time = new_unit.num2date(time_coord.points)
    new_time = [d.replace(day=15) for d in old_time]
    time_coord.points = new_unit.date2num(new_time)


def _prom_dim_coord(cube, _field, _filename):
    iris.util.promote_aux_coord_to_dim_coord(cube, "time")


def _create_coord(cubes, var_name, standard_name):
    cube = cubes.extract_cube(standard_name)
    coord = AuxCoord(
        cube.data,
        standard_name=standard_name,
        long_name=cube.long_name,
        var_name=var_name,
        units="degrees",
    )
    return coord


def _extract_variable(raw_var, cmor_info, attrs, filepath, out_dir, latlon):
    """Extract variable from all files."""
    var = cmor_info.short_name
    cubes = iris.load(filepath, raw_var, _prom_dim_coord)
    iris.util.equalise_attributes(cubes)

    cube = cubes.concatenate_cube()
    iris.util.promote_aux_coord_to_dim_coord(cube, "projection_y_coordinate")
    iris.util.promote_aux_coord_to_dim_coord(cube, "projection_x_coordinate")

    cube.add_aux_coord(latlon[0], (1, 2))
    cube.add_aux_coord(latlon[1], (1, 2))

    # add coord typesi
    area_type = AuxCoord(
        [1.0],
        standard_name="area_type",
        var_name="type",
        long_name="Sea Ice area type",
    )
    cube.add_aux_coord(area_type)

    cube.units = "%"
    cube.data[cube.data > 100] = np.nan
    cube = cube * 100

    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)
    # latlon are multidimensional, create bounds
    siconc = OceanFixGrid(cmor_info)
    cube = siconc.fix_metadata(cubes=[cube])[0]
    # time bounds
    cube.coord("time").bounds = get_time_bounds(
        cube.coord("time"),
        cmor_info.frequency,
    )

    utils.save_variable(
        cube,
        var,
        out_dir,
        attrs,
        unlimited_dimensions=["time"],
    )

    return cube


def _create_areacello(cfg, in_dir, sample_cube, glob_attrs, out_dir):
    if not cfg["custom"].get("create_areacello", False):
        return
    var_info = cfg["cmor_table"].get_variable("Ofx", "areacello")
    glob_attrs["mip"] = "Ofx"
    lat_coord = sample_cube.coord("latitude")

    area_file = os.path.join(in_dir, cfg["custom"]["area_file"])
    with open(area_file, "rb") as datfile:
        areasdmnd = np.fromfile(datfile, dtype=np.int32).reshape(
            lat_coord.shape,
        )

    # Divide by 1000 to get km2 then multiply by 1e6 to m2 ...*1000
    ardata = areasdmnd * 1000

    cube = iris.cube.Cube(
        ardata,
        standard_name=var_info.standard_name,
        long_name=var_info.long_name,
        var_name=var_info.short_name,
        units="m2",
        # time is index 0, add cell index dim
        dim_coords_and_dims=[
            (sample_cube.coords()[1], 0),
            (sample_cube.coords()[2], 1),
        ],
    )
    cube.add_aux_coord(lat_coord, (0, 1))
    cube.add_aux_coord(sample_cube.coord("longitude"), (0, 1))
    utils.fix_var_metadata(cube, var_info)
    utils.set_global_atts(cube, glob_attrs)
    utils.save_variable(
        cube,
        var_info.short_name,
        out_dir,
        glob_attrs,
        zlib=True,
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg["attributes"]
    cmor_table = cfg["cmor_table"]

    # get aux nc file
    cubesaux = iris.load(os.path.join(in_dir, "G02202-cdr-ancillary-sh.nc"))
    lat_coord = _create_coord(cubesaux, "lat", "latitude")
    lon_coord = _create_coord(cubesaux, "lon", "longitude")

    year = 1978
    # split by year..
    sample_cube = None
    for year in range(1979, 2022, 1):
        filepaths = _get_filepaths(in_dir, cfg["filename"], year)

        if len(filepaths) > 0:
            logger.info(
                "Year %d: Found %d files in '%s'",
                year,
                len(filepaths),
                in_dir,
            )

            for var, var_info in cfg["variables"].items():
                logger.info("CMORizing variable '%s'", var)
                glob_attrs["mip"] = var_info["mip"]
                cmor_info = cmor_table.get_variable(var_info["mip"], var)
                raw_var = var_info.get("raw", var)
                sample_cube = _extract_variable(
                    raw_var,
                    cmor_info,
                    glob_attrs,
                    filepaths,
                    out_dir,
                    [lat_coord, lon_coord],
                )

        else:
            logger.info(
                "No files found year: %d basename: %s",
                year,
                cfg["filename"],
            )

    if sample_cube is not None:
        _create_areacello(cfg, in_dir, sample_cube, glob_attrs, out_dir)
