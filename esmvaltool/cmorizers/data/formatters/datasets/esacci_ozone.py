"""ESMValTool CMORizer for ESACCI-OZONE data.

Tier
    Tier 2: other freely-available dataset.

Source

Last access
    20250107

Download and processing instructions

    GTO-ECV total column (variable toz)
    Select the following from the CDS:
    https://cds.climate.copernicus.eu/datasets/satellite-ozone-v1?tab=download

    Processing Level = "Level 3"
    Variable = "Atm. mole content of ozone"
    Vertical aggregation = "Total column"
    Sensor = "MERGED-UV"
    Year = select all (1995-2023)
    Month = select all (1-12)
    Version = "v2000"

    SAGE-CCI-OMPS (variable o3)
    Select the following options from the same link:
    https://cds.climate.copernicus.eu/datasets/satellite-ozone-v1?tab=download

    Processing Level = "Level 3"
    Variable = "Mole concentration of ozone in air"
    Vertical aggregation = " Vertical profiles from limb sensors"
    Sensor = " CMZM (Monthly zonal mean merged concentration product from limb
            sensors ACE, GOMOS, MIPAS, OMPS, OSIRIS, SAGE-2 and SCIAMACHY)"
    Year = select all (1984-2022)
    Month = select all (1-12)
    Version = "v0008"

    Put all files under a single directory (no subdirectories with years).
    in Tier2/ESACCI-OZONE

"""

import logging
from datetime import datetime
from pathlib import Path

import iris
import iris.util
from cf_units import Unit
from esmvalcore.cmor._fixes.native_datasets import NativeDatasetFix
from esmvalcore.preprocessor import concatenate

from ...utilities import (
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def _convert_units(cubes, short_name, var):
    """Perform variable-specific calculations."""
    cube = cubes.extract_cube(var["raw"])
    if short_name == "o3":  # Ozone mole fraction
        gas_constant = 8.31446261815324  # Ideal gas constant (J mol-1 K-1)
        t_cube = cubes.extract_cube("air_temperature")
        p_cube = cubes.extract_cube("air_pressure")
        p_cube.convert_units("Pa")
        air_mol_concentration = p_cube / (gas_constant * t_cube)  # mol m-3
        cube = cube / air_mol_concentration
        cube.units = "mol mol-1"

    elif short_name == "toz":  # Total ozone column (m)
        # Convert from mol m-2 to m
        # -------------------------
        # 1e-5 m (gas @ T = 273 K and p = 101325 Pa) ~ 2.69e20 molecules m-2
        # (see https://en.wikipedia.org/wiki/Dobson_unit)
        # 1 m (O3 @ standard conditions) = 1e5 * 2.69e20 molecules m-2
        #                              = 2.69e25 moelcules m-2 / N_A mol m-2
        #                              = 2.69e25 molecules m-2 /
        #                                    (6.022e23 molecules mol-1) mol m-2
        #                              = 44.67 mol m-2
        cube = cube / 44.67
        cube.units = "m"
    return cube


def _extract_variable(short_name, var, cfg, filename, year, month):
    """Extract variable, add time coordinate, and scalar longitude."""
    mip = var["mip"]
    cmor_info = cfg["cmor_table"].get_variable(mip, short_name)

    cubes = iris.load(filename)

    cube = _convert_units(cubes, short_name, var)

    logger.info("Checking CMOR info for %s: %s", short_name, cmor_info)
    if cmor_info is None:
        raise ValueError(f"CMOR info for {short_name} in MIP {mip} not found!")

    day = 15  # Mid-month
    time_units = Unit("days since 1950-01-01")
    time_points = time_units.date2num(datetime(year, month, day))
    logger.info(
        "Filename: '%s', Extracted Year:'%d', Month: '%d'",
        filename,
        year,
        month,
    )

    # Add time coordinate to cube.
    time_coord = iris.coords.DimCoord(
        time_points,
        var_name="time",
        standard_name="time",
        long_name="time",
        units=time_units,
    )
    time_coord.guess_bounds(monthly=True)
    cube.add_aux_coord(time_coord, ())
    cube = iris.util.new_axis(cube, time_coord)

    # Add longitude coordinate to cube only for o3.
    if short_name == "o3":
        lon_coord = iris.coords.DimCoord(
            [180.0],
            bounds=[[0.0, 360.0]],
            var_name="lon",
            standard_name="longitude",
            long_name="longitude",
            units="degrees_east",
        )
        cube.add_aux_coord(lon_coord, ())
        cube = iris.util.new_axis(cube, lon_coord)
        cube.transpose([1, 3, 2, 0])
        NativeDatasetFix.fix_alt16_metadata(cube)
    fix_var_metadata(cube, cmor_info)
    cube = fix_coords(cube)
    set_global_atts(cube, cfg["attributes"])
    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization process with dataset-specific time ranges."""
    glob_attrs = cfg["attributes"]

    for var_name, var in cfg["variables"].items():
        # Define dataset-specific time ranges
        if var_name == "toz":  # GTO-ECV
            dataset_start = datetime(1995, 7, 1)
            dataset_end = datetime(2023, 4, 30)
        elif var_name == "o3":  # SAGE-CCI-OMPS
            dataset_start = datetime(1984, 10, 1)
            dataset_end = datetime(2022, 12, 31)
        else:
            raise ValueError(f"Unknown dataset for variable {var_name}")

        # Adjust start and end dates if not provided
        start_date_x = start_date or dataset_start
        end_date_x = end_date or dataset_end

        # Ensure the requested date range falls within the dataset limits
        start_date_x = max(start_date_x, dataset_start)
        end_date_x = min(end_date_x, dataset_end)

        all_data_cubes = []
        glob_attrs["mip"] = var["mip"]
        for year in range(start_date_x.year, end_date_x.year + 1):
            for month in range(1, 13):
                # Skip months outside the dataset range
                current_date = datetime(year, month, 1)
                if current_date < dataset_start or current_date > dataset_end:
                    continue

                date_str = f"{year}{month:02}"  # YYYYMM format
                filename = Path(in_dir) / f"{date_str}-{var['filename']}"
                if not filename.is_file():
                    logger.warning(
                        "No file found for %s in %s-%02d",
                        var_name,
                        year,
                        month,
                    )
                    continue

                logger.info(
                    "CMORizing variable '%s' from file '%s'",
                    var_name,
                    filename,
                )
                cube = _extract_variable(
                    var_name,
                    var,
                    cfg,
                    filename,
                    year,
                    month,
                )
                all_data_cubes.append(cube)

        if not all_data_cubes:
            raise ValueError(
                f"No valid data found for {var_name} within the selected time"
                " range",
            )

        final_cube = concatenate(all_data_cubes)
        save_variable(
            final_cube,
            var_name,
            out_dir,
            glob_attrs,
            unlimited_dimensions=["time"],
        )
