# pylint: disable=unused-argument
# pylint: disable=too-many-arguments
# pylint: disable=too-many-function-args
# pylint: disable=R0917
# pylint: disable=E1121
# flake8: noqa
"""ESMValTool CMORizer for TROPFLUX data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://incois.gov.in/tropflux

Last access
    20250115

Download and processing instructions
    Register and sign-in here:

    https://incois.gov.in/tropflux/DataHome.jsp

    Data is available per year from 1979 to 2018
    at daily and monthly frequency.

    Products availables are:

    - Air Temperature at 2m
    - Latent Heat Flux
    - Long Wave Radiation
    - Meridional Wind Stress
    - Net Surface Heat Flux
    - Sea Surface Temperature
    - Sensible Heat Flux
    - Shortwave Radiation
    - Specific Humidity at 2m
    - Wind Speed at 10m
    - Wind Stress magnitude
    - Zonal Wind Stress

Caveats
"""

import logging
import re
from copy import deepcopy
from pathlib import Path

import cf_units
import iris
from iris import NameConstraint
from iris.util import equalise_attributes

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

# scitools-iris.readthedocs.io/en/stable/generated/api/iris.html#iris.Future
try:
    iris.FUTURE.date_microseconds = True
    iris.FUTURE.save_split_attrs = True
except AttributeError as e:
    # Handle cases where FUTURE or the attributes don't exist
    logger.warning("AttributeError: %s", e)
except (TypeError, ValueError) as e:
    # Handle specific errors if these might occur
    logger.warning("TypeError or ValueError: %s", e)
except BaseException as e:
    # Fallback for rare or unknown issues, but avoid catching Exception
    logger.warning("An unexpected error occurred: %s", e)


def _fix_coordinates(cube, definition):
    cube = utils.fix_coords(cube)

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


def _extract_variable(short_name, var, cfg, raw_filepaths, out_dir):
    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip"]
    definition = cfg["cmor_table"].get_variable(var["mip"], short_name)

    if definition.positive != "":
        attributes["positive"] = definition.positive

    # load data
    raw_var = var.get("raw", short_name)
    cubes = iris.load(raw_filepaths, NameConstraint(var_name=raw_var))
    equalise_attributes(cubes)
    for cube in cubes:
        cube.rename(cubes[0].name())
    cube = cubes.concatenate_cube()
    utils.set_global_atts(cube, attributes)

    # Set correct names
    cube.var_name = definition.short_name
    if definition.standard_name:
        cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    utils.fix_var_metadata(cube, definition)

    # fix time units
    cube.coord("time").convert_units(
        cf_units.Unit("days since 1950-1-1 00:00:00", calendar="gregorian")
    )

    cube = _fix_coordinates(cube, definition)

    if raw_var in ["taux", "lhf"]:
        cube.data = -1 * cube.data

    utils.save_variable(
        cube,
        short_name,
        out_dir,
        attributes,
        unlimited_dimensions=["time"],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Run CMORizer for TROPFLUX."""
    for short_name, var in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", short_name)
        short_name = var["short_name"]
        raw_filenames = Path(in_dir).rglob("*.nc")
        filenames = []
        for raw_filename in raw_filenames:
            if re.search(var["file"], str(raw_filename)):
                filenames.append(raw_filename)

        _extract_variable(short_name, var, cfg, filenames, out_dir)
