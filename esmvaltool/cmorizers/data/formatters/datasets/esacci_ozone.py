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

    MEGRIDOP (variable o3)
    Processing Level = "Level 3"
    Variable = "At, model content of ozone"
    Vertical aggregation = "Vertical profiles from limb sensors"
    Sensor = "CLLG"
    Year = select all (2001-2024)
    Month = select all (1-12)
    Version = "v0005"

    Put all files under a single directory (no subdirectories with years).
    in ${RAWOBS}/Tier2/ESACCI-OZONE

    ---------------------------------------------------------------------------

    IASI (variable o3)
    Download from BIRA WebDAV server: https://webdav.aeronomie.be
    Path: /guest/o3_cci/webdata/Nadir_Profiles/L3/IASI_MG_FORLI/
    Username: o3_cci_public
    No password (leave empty)
    Download each year (yyyy) into separate folders named "IASI_yyyy"
"""

import glob
import logging
import os
from datetime import datetime
from pathlib import Path
import dask.array as da
import numpy as np

import iris
import iris.util
import iris.experimental.stratify
from cf_units import Unit
from esmvalcore.cmor._fixes.native_datasets import NativeDatasetFix
from esmvalcore.preprocessor import concatenate, monthly_statistics

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
        if var["var_name"] == "o3_iasi":  # IASI merged profiles
            air_mol_concentration = cubes.extract_cube(
                "air_partial_column_profile")
            cube = cube / air_mol_concentration
            cube.units = "mol mol-1"
        else: # SAGE-CCI-OMPS or MEGRIDOP profiles
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


def _extract_variable(in_files, var, cfg, out_dir, year, month):
    """Extract variable, add time coordinate, and scalar longitude."""
    short_name = var["output"]
    mip = var["mip"]
    cmor_info = cfg["cmor_table"].get_variable(mip, short_name)

    # add time coordinate(s)
    # try to extract day from filename(s)
    new_list = iris.cube.CubeList()
    for fname in in_files:
        print(fname)
        cubes = iris.load(fname)
        i = fname.find(f"{year}{month:02}")
        strday = fname[i+6:i+8]
        try:
            day = int(strday)
        except ValueError:
            day = 15

        time_units = Unit("days since 1950-01-01")
        time_points = time_units.date2num(datetime(year, month, day))

        # Add time coordinate to cube.
        time_coord = iris.coords.DimCoord(
            time_points,
            var_name="time",
            standard_name="time",
            long_name="time",
            units=time_units,
        )
        for cube in cubes:
            cube.add_aux_coord(time_coord, ())
            cube = iris.util.new_axis(cube, time_coord)
            new_list.append(cube)

    logger.info("Checking CMOR info for %s: %s", short_name, cmor_info)
    if cmor_info is None:
        raise ValueError(f"CMOR info for {short_name} in MIP {mip} not found!")

    cubes = new_list.concatenate()
    cube = _convert_units(cubes, short_name, var)
    # fix nan's
    cube.data = da.ma.fix_invalid(cube.core_data())
    # calculate monthly means (IASI)
    cube = monthly_statistics(cube, operator="mean")

    # Add longitude coordinate to cube only for o3_sage_omps.
    if var["var_name"] == "o3_sage_omps":
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
    elif var["var_name"] == "o3_megridop":
        NativeDatasetFix.fix_alt16_metadata(cube)
    # add latitude, longitude and nlev coordiantes to cube for o3_iasi
    elif var["var_name"] == "o3_iasi":
        # add named coordinates
        # longitude
        lon_cube = cubes.extract_cube("longitude")
        lon_coord = iris.coords.DimCoord(
            lon_cube.core_data()[0,:],
            var_name="lon",
            standard_name="longitude",
            long_name="longitude",
            units="degrees_east",
        )
        cube.add_aux_coord(lon_coord, 2)
        iris.util.promote_aux_coord_to_dim_coord(cube, "longitude")
        # latitude
        lat_cube = cubes.extract_cube("latitude")
        lat_coord = iris.coords.DimCoord(
            lat_cube.core_data()[0,:],
            var_name="lat",
            standard_name="latitude",
            long_name="latitude",
            units="degrees_north",
        )
        cube.add_aux_coord(lat_coord, 3)
        iris.util.promote_aux_coord_to_dim_coord(cube, "latitude")
        # level
#        lev_coord = iris.coords.DimCoord(
#            np.arange(0, 1, 1/cube.shape[1]),
#            var_name="lev",
#            standard_name="atmosphere_hybrid_sigma_pressure_coordinate",
#            long_name="hybrid sigma pressure coordinate",
#            units="1",
#            attributes={'positive': 'down'},
#        )
        lev_coord = iris.coords.DimCoord(
            np.arange(1, 0, -1/cube.shape[1]),
            var_name="lev",
            standard_name="atmosphere_hybrid_sigma_pressure_coordinate",
            long_name="hybrid sigma pressure coordinate",
            units="1",
            attributes={'positive': 'up'},
        )
        cube.add_dim_coord(lev_coord, 1)
###        cube.add_aux_coord(lev_coord, 1)
###        iris.util.promote_aux_coord_to_dim_coord(cube, "atmosphere_hybrid_sigma_pressure_coordinate")
###        cube.add_dim_coord(lev_coord, 1)
        # reorder dimensions to [time, lev, lat, lon]
        cube.transpose([0,1,3,2])
        # air pressure (aux coordinate)
        # calculate full levels from half levels
        ap_half = cubes.extract_cube("atmosphere_pressure_grid")
        ap_half.data = np.ma.masked_invalid(ap_half.data)
        ap_half = monthly_statistics(ap_half, operator="mean")
###        ap_half.data = da.ma.masked_less(ap_half.core_data(), 1)
###        ap_half = ap_half.collapsed('time', iris.analysis.MEAN)
        # reoder air pressure (aux coordinate) before adding to cube
        # (otherwise, data of aux coordiante remains in original order)
        ap_half.transpose([0,1,3,2])
        ap_full = cube.copy()
        for k in range(0, ap_full.shape[1]):
            ap_full.core_data()[:,k,:,:] = 0.5 * (ap_half.core_data()[:,k,:,:] + ap_half.core_data()[:,k+1,:,:])
#            da.ma.masked_invalid(ap_full.core_data())
#        ap_full = np.ma.masked_where(ap_full.data <= 1, ap_full.data, copy=False)
#        ap_full.fill_value = 1e20
        ap_means = ap_full.core_data().mean(axis=[2, 3], keepdims=True)
        ap_full.data = da.where(da.ma.getmaskarray(ap_full.core_data()), ap_means, ap_full.core_data())
#        ap_full.data = da.array(ap_full.core_data())
        ap_coord = iris.coords.AuxCoord(
            ap_full.core_data(),
            var_name="plev",
            standard_name="air_pressure",
            long_name="pressure",
            units="Pa",
#            attributes={'_FillValue': '1e20'},
        )
        cube.add_aux_coord(ap_coord, [0,1,2,3])

    fix_var_metadata(cube, cmor_info)
###    cube = fix_coords(cube)
    set_global_atts(cube, cfg["attributes"])
    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization process with dataset-specific time ranges."""
    glob_attrs = cfg["attributes"]
    if "version" in glob_attrs:
        glob_version = glob_attrs["version"]
    else:
        glob_version = ""

    for var_name, var in cfg["variables"].items():
        # Define dataset-specific time ranges
        if var_name == "toz_gto_ecv":  # GTO-ECV
            dataset_start = datetime(1995, 7, 1)
            dataset_end = datetime(2023, 4, 30)
        elif var_name == "o3_sage_omps":  # SAGE-CCI-OMPS
            dataset_start = datetime(1984, 10, 1)
            dataset_end = datetime(2022, 12, 31)
        elif var_name == "o3_megridop":  # MEGRIDOP
            dataset_start = datetime(2001, 11, 1)
            dataset_end = datetime(2023, 12, 31)
        elif var_name == "o3_iasi":  # IASI
            dataset_start = datetime(2008, 1, 1)
            dataset_end = datetime(2008, 12, 31)  #(2023, 12, 31)
        else:
            raise ValueError(f"Unknown dataset for variable {var_name}")

        var["var_name"] = var_name

        # Adjust start and end dates if not provided
        start_date_x = start_date or dataset_start
        end_date_x = end_date or dataset_end

        # Ensure the requested date range falls within the dataset limits
        start_date_x = max(start_date_x, dataset_start)
        end_date_x = min(end_date_x, dataset_end)

        all_data_cubes = []
        glob_attrs["mip"] = var["mip"]
        if "version" in var:
            glob_attrs["version"] = var["version"]
        else:
            glob_attrs["version"] = glob_version
        output_var = var["output"]
        for year in range(start_date_x.year, end_date_x.year + 1):
            if var_name == "o3_iasi":
                subfolder = f"IASI_{year}"
            else:
                subfolder = ""
            for month in range(1, 13):
                # Skip months outside the dataset range
                current_date = datetime(year, month, 1)
                if current_date < dataset_start or current_date > dataset_end:
                    continue

                date_str = f"{year}{month:02}"  # YYYYMM format
                monstr = f"{month:02}"
                
                filepattern = os.path.join(
                    in_dir, subfolder, var['filename'].format(year=year, month=monstr)
                )
                in_files = glob.glob(filepattern)
                if not in_files:
                    logger.info(f'{var_name}: no data not found for '
                                f'{year}-{monstr}')
                    continue
                else:
                    cube = _extract_variable(in_files, var, cfg, out_dir, year, month)

                logger.info(
                    "CMORizing variable '%s' from file(s) '%s'",
                    output_var,
                    in_files,
                )
                all_data_cubes.append(cube)

        if not all_data_cubes:
            raise ValueError(
                f"No valid data found for {var_name} within the selected time"
                " range"
            )

        final_cube = concatenate(all_data_cubes)
        save_variable(
            final_cube,
            output_var,
            out_dir,
            glob_attrs,
            unlimited_dimensions=["time"],
        )
