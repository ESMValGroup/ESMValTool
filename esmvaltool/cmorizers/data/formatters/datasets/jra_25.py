"""
ESMValTool CMORizer for JRA-25 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://data.rda.ucar.edu/ds625.1
Last access
    2023-06-14

Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/jra_25.py
"""

import copy
import glob
import logging
import os
from datetime import datetime

import iris
from dateutil import relativedelta
from esmvalcore.preprocessor import extract_levels
from iris import NameConstraint
from iris.util import equalise_attributes

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, in_files, cfg, out_dir):
    """Extract variable."""
    # load data
    raw_var = var.get("raw", short_name)
    rawcubes = iris.load(in_files, NameConstraint(var_name=raw_var))

    equalise_attributes(rawcubes)

    # check if data are on hybrid levels

    hybrid_levels = False

    for coord in rawcubes[0].coords():
        if (
            coord.standard_name
            == "atmosphere_hybrid_sigma_pressure_coordinate"
        ):
            hybrid_levels = True
            break

    # If (3-dim) data are on hybrid levels then regrid data to pressure levels.
    # The pressure levels are taken from the JRA-25 files that provide JRA-25
    # data on pressure levels (e.g. anl_p.yyyymm.nc).

    if hybrid_levels:
        pcubes = iris.cube.CubeList([])

        for cube in rawcubes:
            pcube = extract_levels(
                cube,
                [
                    100000,
                    92500,
                    85000,
                    70000,
                    60000,
                    50000,
                    40000,
                    30000,
                    25000,
                    20000,
                    15000,
                    10000,
                    7000,
                    5000,
                    3000,
                    2000,
                    1000,
                    700,
                    500,
                    300,
                    200,
                    100,
                    40,
                ],
                "linear",
                coordinate="air_pressure",
            )
            # remove auxiliary coordinate 'Surface_pressure'
            pcube.remove_coord("Surface_pressure")
            # rename dimension air_pressure to plev
            coord = pcube.coord("air_pressure")
            coord.rename("plev")
            coord.var_name = "plev"

            pcubes.append(pcube)

        cube = pcubes.concatenate_cube()
    else:
        cube = rawcubes.concatenate_cube()

    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)

    try:
        cube.convert_units(cmor_info.units)
    except Exception as ex:  # pylint: disable=broad-except
        logger.warning(
            "Warning: could not convert units from %s to %s (%r)",
            cube.units,
            cmor_info.units,
            ex,
        )

    # Fix coordinates
    utils.fix_coords(cube)
    for coord in cube.coords():
        if coord.var_name == "plev":
            coord.standard_name = "air_pressure"
            coord.long_name = "pressure"
            coord.attributes["positive"] = "up"

    # Fix metadata
    attrs = copy.deepcopy(cfg["attributes"])
    attrs["mip"] = var["mip"]
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Add height2m if needed
    if "height2m" in cmor_info.dimensions:
        utils.add_height2m(cube)

    # Save variable
    utils.save_variable(
        cube,
        short_name,
        out_dir,
        attrs,
        unlimited_dimensions=["time"],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    if start_date is None:
        start_date = datetime(1979, 1, 1)
    if end_date is None:
        end_date = datetime(2007, 12, 31)
    loop_date = start_date

    while loop_date <= end_date:
        year = loop_date.year

        for short_name, var in cfg["variables"].items():
            if "short_name" not in var:
                var["short_name"] = short_name

            # Now get list of files
            filepattern = os.path.join(
                in_dir + "/" + var["file"].format(year=year),
            )
            print("*** " + filepattern)
            in_files = glob.glob(filepattern)
            if not in_files:
                logger.warning("Warning: no data found for %d", year)
                continue
            _extract_variable(short_name, var, in_files, cfg, out_dir)

        loop_date += relativedelta.relativedelta(years=1)
