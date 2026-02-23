"""
ESMValTool CMORizer for JRA-55 data.

Tier
    Tier 2: other freely-available dataset.

Source
    Research Data Archive (RDA):
    https://rda.ucar.edu/datasets/ds628.1/

Last access
    20230322

Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/jra_55.py
"""

import copy
import logging
import os

import iris
import xarray as xr
from cf_units import Unit

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _load_jra55_grib(filenames, var):
    """Load data from GRIB file and return list of cubes."""
    leveltype = var.get("typeOfLevel")
    cubelist = []
    if leveltype is not None:
        dataset = xr.open_mfdataset(
            filenames,
            engine="cfgrib",
            filter_by_keys={"typeOfLevel": leveltype},
        )
    else:
        dataset = xr.open_mfdataset(filenames, engine="cfgrib")
    varnames = list(dataset.data_vars)
    for varname in varnames:
        da_tmp = dataset[varname]
        # conversion to Iris cubes requires a valid standard_name
        da_tmp.attrs["standard_name"] = var["standard_name"]
        cube = da_tmp.to_iris()
        # remove auxiliary coordinate 'time'
        cube.remove_coord("time")
        # rename coordinate from 'forecast_reference_time' to 'time
        timecoord = cube.dim_coords[0]
        timecoord.rename("time")
        # convert unit string to cf_unit object
        # (calendar (calendar=coord.units.calendar) must be irgnored
        # or conversion fails
        timecoord.units = Unit(timecoord.units)
        # add forecast period to time coordinate to get the actual time
        # for which the data are valid
        forecast = cube.coord("forecast_period")  # forecast period in hours
        timecoord.points = timecoord.points + forecast.points * 3600
        # remove unneeded scalar variables to prevent warnings
        auxcoordnames = [
            "step",
            "entireAtmosphere",
            "number",
            "isobaricLayer",
            "surface",
            "nominalTop",
            "heightAboveGround",
        ]
        for aux_coord in cube.coords(dim_coords=False):
            if aux_coord.var_name in auxcoordnames:
                cube.remove_coord(aux_coord)
        cubelist.append(cube)

    return cubelist


def _extract_variable(short_name, var, in_files, cfg, out_dir):
    """Extract variable."""
    # load data (returns a list of cubes)
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)
    var["standard_name"] = cmor_info.standard_name
    cubes = _load_jra55_grib(in_files, var)

    # apply operators (if any)
    if len(cubes) > 1:
        if var.get("operator", "") == "sum":
            # Multiple variables case using sum operation
            cube = None
            for in_cube in cubes:
                if cube is None:
                    cube = in_cube
                else:
                    cube += in_cube
        elif var.get("operator", "") == "diff":
            # two variables case using diff operation
            if len(cubes) != 2:
                errmsg = (
                    f"operator diff selected for variable {short_name} "
                    f"expects exactly two input variables and two input "
                    f"files"
                )
                raise ValueError(errmsg)
            cube = cubes[0] - cubes[1]
        else:
            oper = var.get("operator")
            raise ValueError(
                f"multiple input files found for variable {short_name} "
                f"with unknown operator {oper}",
            )
    else:
        cube = cubes[0]

    # Fix metadata
    attrs = copy.deepcopy(cfg["attributes"])
    attrs["mip"] = var["mip"]
    utils.fix_var_metadata(cube, cmor_info)

    if cube.var_name in [
        "hfls",
        "hfss",
        "rlus",
        "rlut",
        "rlutcs",
        "rsus",
        "rsuscs",
        "rsut",
        "rsutcs",
    ]:
        attrs["positive"] = "up"

    if cube.var_name in [
        "rlds",
        "rldscs",
        "rsds",
        "rsdscs",
        "rsdt",
        "rtmt",
        "tauu",
        "tauv",
    ]:
        attrs["positive"] = "down"

    # fix longitudes and z-coordinate (if present)
    for coord in cube.dim_coords:
        coord_type = iris.util.guess_coord_axis(coord)
        if coord_type == "X":
            # -> shift longitude coordinate by one grid box
            # to match obs4mips/CREATE-IP grid
            coord.points = coord.points + 360 / len(coord.points)
        if coord_type == "Z":
            coord.standard_name = "air_pressure"
            coord.long_name = "pressure"
            coord.var_name = "plev"
            coord.attributes["positive"] = "down"
            if coord.units == "hPa":
                coord.convert_units("Pa")
            utils.flip_dim_coord(cube, coord.standard_name)

    utils.fix_dim_coordnames(cube)
    utils.fix_coords(cube)
    if "height2m" in cmor_info.dimensions:
        utils.add_height2m(cube)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(
        cube,
        short_name,
        out_dir,
        attrs,
        unlimited_dimensions=["time"],
        local_keys=["positive"],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    if start_date is None:
        start_date = 1958
    else:
        start_date = start_date.year
    if end_date is None:
        end_date = 2023
    else:
        end_date = end_date.year
    for short_name, var in cfg["variables"].items():
        short_name = var["short_name"]
        filename = []
        for year in range(start_date, end_date + 1):
            if "file" in var:
                filename.append(
                    os.path.join(in_dir, var["file"].format(year=year)),
                )
            elif "files" in var:
                for file in var["files"]:
                    filename.append(
                        os.path.join(in_dir, file.format(year=year)),
                    )
            else:
                raise ValueError(
                    f"No input file(s) specified for variable {short_name}.",
                )

        logger.info(
            "CMORizing variable '%s' from file '%s'",
            short_name,
            filename,
        )
        _extract_variable(short_name, var, filename, cfg, out_dir)
