"""ESMValTool CMORizer for ESACCI-SOILMOISTURE data.

Tier
    Tier 2: other freely-available dataset.

Source
    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/soil_moisture/data/

Last access
    20240626

Download and processing instructions
    Download the data from:
      daily_files/COMBINED/v08.1/
      ancillary/v08.1/
    Put all files under a single directory (no subdirectories with years).
    in ${RAWOBS}/Tier2/ESACCI-SOILMOISTURE

"""

import glob
import logging
import os
from datetime import datetime

import iris
from cf_units import Unit
from esmvalcore.preprocessor import concatenate, monthly_statistics

from ...utilities import (
    fix_bounds,
    fix_dim_coordnames,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def fix_coords(cube):
    """Fix coordinates to CMOR standards.

    Fixes coordinates eg time to have correct units, bounds etc;
    longitude to be CMOR-compliant 0-360deg; fixes some attributes
    and bounds - the user can avert bounds fixing by using supplied
    arguments; if bounds are None they will be fixed regardless.

    Parameters
    ----------
    cube: iris.cube.Cube
        data cube with coordinates to be fixed.


    Returns
    -------
    cube: iris.cube.Cube
        data cube with fixed coordinates.
    """
    # First fix any completely missing coord var names
    fix_dim_coordnames(cube)

    # Convert longitude from -180...180 to 0...360
    cube = cube.intersection(longitude=(0.0, 360.0))

    # Fix individual coords
    for cube_coord in cube.coords():
        # Fix time
        if cube_coord.var_name == "time":
            logger.info("Fixing time...")
            cube.coord("time").convert_units(
                Unit(
                    "days since 1970-01-01T00:00:00+00:00",
                    calendar="proleptic_gregorian",
                )
            )

        # Fix latitude
        if cube_coord.var_name == "lat":
            logger.info("Fixing latitude...")
            cube = iris.util.reverse(cube, cube_coord)

        # Fix bounds of all coordinates
        fix_bounds(cube, cube_coord)

    return cube


def extract_variable(raw_info):
    """Extract variables."""
    rawvar = raw_info["name"]
    constraint = iris.Constraint(name=rawvar)
    if rawvar == "sm_uncertainty":
        sm_cube = iris.load_cube(
            raw_info["file"], iris.NameConstraint(var_name="sm")
        )
        ancillary_var = sm_cube.ancillary_variable(
            "Volumetric Soil Moisture Uncertainty"
        )
        cube = sm_cube.copy(ancillary_var.core_data())
    else:
        cube = iris.load_cube(raw_info["file"], constraint)

    # Remove dysfunctional ancillary data without standard names
    for ancillary_variable in cube.ancillary_variables():
        cube.remove_ancillary_variable(ancillary_variable)

    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    glob_attrs = cfg["attributes"]
    if not start_date:
        start_date = datetime(1978, 1, 1)
    if not end_date:
        end_date = datetime(2022, 12, 31)

    # run the cmorization
    for var_name, vals in cfg["variables"].items():
        all_data_cubes = []
        if not isinstance(vals, dict):  # Ensure vals is a dictionary
            raise ValueError(
                f"Invalid format for variable {var_name}: {type(vals)}"
            )
        var_info = cfg["cmor_table"].get_variable(vals["mip"], var_name)
        glob_attrs["mip"] = vals["mip"]
        raw_info = {"name": vals["raw"]}
        inpfile_pattern = os.path.join(in_dir, vals["filename"])
        logger.info(
            "CMORizing var %s from file type %s", var_name, inpfile_pattern
        )

        for year in range(start_date.year, end_date.year + 1):
            year_inpfile_pattern = inpfile_pattern.format(year=year)
            inpfiles = sorted(glob.glob(year_inpfile_pattern))
            for inpfile in inpfiles:
                raw_info["file"] = inpfile
                cube = extract_variable(raw_info)
                all_data_cubes.append(cube)
        final_cube = concatenate(all_data_cubes)
        fix_var_metadata(final_cube, var_info)
        final_cube = fix_coords(final_cube)
        set_global_atts(final_cube, glob_attrs)

        save_variable(
            final_cube,
            var_name,
            out_dir,
            glob_attrs,
            unlimited_dimensions=["time"],
        )

        # For sm, also save monthly means
        if var_name == "sm":
            monthly_mean_cube = monthly_statistics(final_cube, "mean")
            glob_attrs["mip"] = "Lmon"
            monthly_mean_cube.attributes.update(glob_attrs)
            save_variable(
                monthly_mean_cube,
                var_name,
                out_dir,
                glob_attrs,
                unlimited_dimensions=["time"],
            )
