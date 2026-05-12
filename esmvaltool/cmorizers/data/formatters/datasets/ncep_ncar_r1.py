"""ESMValTool CMORizer for NCEP-NCAR-R1 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html

Last access
    20221116

Download and processing instructions
    To facilitate the download, the links to the ftp server are provided.
    Since the filenames are sometimes identical across different
    save the data in two subdirectories in input_dir_path.
    Subdirectory pressure/:
      ftp://ftp.cdc.noaa.gov/Projects/Datasets/data.ncep.reanalysis/pressure/
        air.mon.mean.nc
        hgt.mon.mean.nc
        rhum.mon.mean.nc
        shum.mon.mean.nc
        uwnd.mon.mean.nc
        vwnd.mon.mean.nc
        omega.mon.mean.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/
        uwnd.*.nc
        vwnd.*.nc

    Subdirectory surface/:
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/
        air.mon.mean.nc
        pr_wtr.mon.mean.nc
        slp.mon.mean.nc
        wspd.mon.mean.nc
        rhum.mon.mean.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/
        air.2m.mon.mean.nc
        prate.sfc.mon.mean.nc
        tmax.2m.mon.mean.nc
        tmin.2m.mon.mean.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/other_gauss/
        tcdc.eatm.mon.mean.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/
        prate.sft.gauss.*.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/other_gauss/
        ulwrf.ntat.gauss.*.nc


    #Select the section "Pressure" and "Surface" and download the variables
    #listed below. Since raw data on pressure levels and for surface have the
    #same file and variable name, save the data in two different subdirectories
    #"press" and "surf" in input_dir_path.

Caveats

"""

import logging
import re
from copy import deepcopy
from pathlib import Path
from warnings import catch_warnings, filterwarnings

import iris
from cf_units import Unit
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _fix_units(cube, definition):
    """
    Fix issues with the units.

    Exception is `pr` since the units in the
    raw file are not recognized correctly.
    """
    if cube.var_name != "pr":
        cube.convert_units(definition.units)


def _fix_coordinates(cube, definition):
    cube = utils.fix_coords(cube)

    if "height2m" in definition.dimensions:
        utils.add_height2m(cube)
    if "height10m" in definition.dimensions:
        utils.add_scalar_height_coord(cube, height=10.0)

    for coord_def in definition.coordinates.values():
        axis = coord_def.axis
        coord = cube.coord(axis=axis)
        if axis == "Z":
            coord.convert_units(coord_def.units)
        coord.standard_name = coord_def.standard_name
        coord.var_name = coord_def.out_name
        coord.long_name = coord_def.long_name
        coord.points = coord.core_points().astype("float64")
        if coord.var_name == "plev":
            coord.attributes["positive"] = "down"

    return cube


def _extract_variable(short_name, var, cfg, raw_filepath, out_dir):
    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip"]
    definition = cfg["cmor_table"].get_variable(var["mip"], short_name)
    if definition.positive != "":
        attributes["positive"] = definition.positive

    # load data
    raw_var = var.get("raw", short_name)
    with catch_warnings():
        filterwarnings(
            "ignore",
            message="Ignoring netCDF variable .* invalid units .*",
            category=UserWarning,
            module="iris",
        )
        cube = iris.load_cube(
            str(raw_filepath),
            NameConstraint(var_name=raw_var),
        )

    utils.set_global_atts(cube, attributes)

    # Set correct names
    cube.var_name = definition.short_name
    if definition.standard_name:
        cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    _fix_units(cube, definition)

    utils.fix_var_metadata(cube, definition)

    # fix time units
    cube.coord("time").convert_units(
        Unit("days since 1950-1-1 00:00:00", calendar="gregorian"),
    )

    cube = _fix_coordinates(cube, definition)

    if var.get("make_negative"):
        cube.data = -1 * cube.data

    utils.save_variable(
        cube,
        short_name,
        out_dir,
        attributes,
        unlimited_dimensions=["time"],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Run CMORizer for NCEP-NCAR-R1."""
    # Run the cmorization
    for short_name, var in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", short_name)
        short_name = var["short_name"]
        raw_filenames = Path(in_dir).rglob("*.nc")
        filenames = []
        for raw_filename in raw_filenames:
            if re.search(var["file"], str(raw_filename)) is not None:
                filenames.append(raw_filename)

        for filename in sorted(filenames):
            _extract_variable(short_name, var, cfg, filename, out_dir)
