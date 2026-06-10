"""ESMValTool CMORizer for ESACCI-BIOMASS above-ground biomass (agb) data.

Tier
    Tier 2: other freely-available dataset.

Source
    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps

Last access
    20250716

Download and processing instructions
    Download 10 km file:
      v6.0/netcd/ESACCI-BIOMASS-L4-AGB-MERGED-10000m-fv6.0.nc
    Put file in Tier2/ESACCI-BIOMASS
    Or use automatic download script:
      esmvaltool data download ESACCI-BIOMASS
"""

import datetime
import glob
import logging
import os
from copy import deepcopy

import iris
import numpy as np
from dask import array as da
from esmvalcore.preprocessor import extract_time

from esmvaltool.cmorizers.data.utilities import (
    flip_dim_coord,
    save_variable,
    set_global_atts,
)

# Configure logging
logger = logging.getLogger(__name__)


def _extract_variable(in_files, var, cfg, out_dir):
    logger.info(
        "CMORizing variable '%s' from input files '%s'",
        var["short_name"],
        ", ".join(in_files),
    )
    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip"]
    attributes["raw"] = var["raw"]
    attributes["frequency"] = var["frequency"]
    definition = cfg["cmor_table"].get_variable(var["mip"], var["short_name"])

    # load all input files (1 year) into 1 cube
    cube_list = iris.load(in_files, var["raw"])

    drop_attrs = ["valid_max", "valid_min"]

    for cube in cube_list:
        # set correct names
        cube.var_name = definition.short_name
        cube.standard_name = definition.standard_name
        cube.long_name = definition.long_name

        for attr in drop_attrs:
            if attr in cube.attributes:
                cube.attributes.pop(attr)

        cube.coord("time").points = (
            cube.coord("time").core_points().astype("float64")
        )

        # fix units
        cube.convert_units(definition.units)

        # set global attributes
        set_global_atts(cube, attributes)

        # roll longitude (-180...180 --> 0...360)
        cube.coord("longitude").points = cube.coord("longitude").points + 180.0
        nlon = len(cube.coord("longitude").points)
        cube.data = da.roll(cube.core_data(), int(nlon / 2), axis=2)

        # remove rouding errors introduced by da.roll
        loncoord = cube.coord("longitude")
        latcoord = cube.coord("latitude")
        loncoord.points = np.round(loncoord.core_points(), 3)
        latcoord.points = np.round(latcoord.core_points(), 3)

        # flip latitudes
        flip_dim_coord(cube, "latitude")

        # fix coordinates
        cube = _fix_coordinates(cube, definition)
        cube.coord("latitude").attributes = None
        cube.coord("longitude").attributes = None

        # save each year to a separate output file
        timecoord = cube.coord("time")
        for time in timecoord.units.num2date(timecoord.points):
            # extract current year
            outcube = extract_time(cube, time.year, 1, 1, time.year, 12, 31)
            # adjust time bounds to (year-01-01 00:00, year+1-01-01 00:00)
            out_timecoord = outcube.coord("time")
            start_date = datetime.datetime(time.year, 1, 1)
            end_date = datetime.datetime(time.year + 1, 1, 1)
            out_timecoord.bounds = np.array(
                [
                    out_timecoord.units.date2num(start_date),
                    out_timecoord.units.date2num(end_date),
                ],
            )
            # write output to file
            logger.debug("Saving cube\n%s", outcube)
            logger.debug("Setting time dimension to UNLIMITED while saving!")
            save_variable(
                outcube,
                var["short_name"],
                out_dir,
                attributes,
                unlimited_dimensions=["time"],
            )

    logger.info("Finished CMORizing %s", ", ".join(in_files))


def _fix_coordinates(cube, definition):
    """Fix coordinates."""
    axis2def = {"T": "time", "X": "longitude", "Y": "latitude"}
    axes = ["T", "X", "Y"]

    for axis in axes:
        coord_def = definition.coordinates.get(axis2def[axis])
        if coord_def:
            coord = cube.coord(axis=axis)
            if axis == "T":
                coord.convert_units("days since 1850-1-1 00:00:00.0")
                coord.points = coord.core_points().astype("float64")
                if coord.bounds is not None:
                    coord.bounds = None

            if len(coord.points) > 1:
                if coord.bounds is not None:
                    coord.bounds = None
                coord.guess_bounds()
            coord.standard_name = coord_def.standard_name
            coord.var_name = coord_def.out_name
            coord.long_name = coord_def.long_name

    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    glob_attrs = cfg["attributes"]

    logger.info(
        "Starting cmorization for tier%s OBS files: %s",
        glob_attrs["tier"],
        glob_attrs["dataset_id"],
    )
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    logger.info("CMORizing ESACCI-BIOMASS version %s", glob_attrs["version"])

    for short_name, var in cfg["variables"].items():
        filepattern = os.path.join(in_dir, var["filename"])
        in_files = glob.glob(filepattern)
        if "short_name" not in var:
            var["short_name"] = short_name
        if not in_files:
            msg = f"no data not found for variable {short_name}"
            raise ValueError(msg)
        _extract_variable(in_files, var, cfg, out_dir)
