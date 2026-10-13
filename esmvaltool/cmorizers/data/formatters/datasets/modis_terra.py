"""ESMValTool CMORizer for MODIS terra data.

Tier
    Tier 3

Source
    https://ladsweb.modaps.eosdis.nasa.gov/search/order

Requirements
  python hdf reader, https://pypi.org/project/pyhdf/

Modification history
    20260216-Fruttarol_Noah: Written based on the original MODIS CMORization NCL script created by evaldsson_martin, which can be found in at ESMValTool/esmvaltool/cmorizers/data/formatters/datasets/_modis_terra.ncl


"""

import logging
import os
import re
from datetime import datetime, timedelta

import iris
import numpy as np
from dask import array as da
from esmvalcore.cmor.table import CMOR_TABLES
from pyhdf.SD import SD, SDC

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def collect_files(in_dir: str, cfg) -> list:
    """Compose input file list and download if missing."""

    attrs = cfg["attributes"]
    # rootpath / Tier{tier}/{dataset}/{version}/{frequency}/{short_name}
    # in_dir  = rootpath / Tier{tier}/{dataset}
    in_dir = os.path.join(in_dir)
    fname_pattern = (
        f"{attrs['version']}\\.A.*\\.hdf$"  # e.g. "MOD08_M3.A*.hdf"
    )

    files_a = next(os.walk(in_dir), (None, None, []))[
        2
    ]  # get list of files in in_dir
    files = []
    for filename in files_a:
        if re.search(fname_pattern, filename):
            files.append(filename)

    if len(files) == 0:
        raise Exception("No files found", in_dir)

    return files


def get_year_from_filepath(filepath: str, cfg) -> int:
    """get year from first four characters of last section of basename"""

    attrs = cfg["attributes"]
    basename = os.path.splitext(os.path.basename(filepath))[0]
    offset = len(
        f"{attrs['version']}.A"
    )  # length of version string + dot + 'A' character
    year = int(basename[offset : offset + 4])
    return year


def get_month_from_filepath(filepath: str, cfg) -> int:
    """get month from timerange section of basename"""

    attrs = cfg["attributes"]
    basename = os.path.splitext(os.path.basename(filepath))[0]
    offset = len(
        f"{attrs['version']}.A"
    )  # length of version string + dot + 'A' character
    doy = int(basename[offset + 4 : offset + 7])  # day of year
    year = int(basename[offset : offset + 4])
    # Convert day-of-year to month
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    month = date.month
    logger.debug(
        f"Extracted month {month}, doy {doy}, year {year} from filepath {filepath}"
    )
    return month


def group_files_by_year(filepaths: list, cfg) -> dict:
    """group filepaths by year using get_year_from_filepath function"""
    years_D = dict()
    for filename in filepaths:
        year = get_year_from_filepath(filename, cfg)

        if year not in years_D:
            years_D[year] = [filename]
        else:
            years_D[year].append(filename)

    return years_D


def read_hdf(
    in_dir: str,
    filepath: str,
    year: int,
    month: int,
    raw_name: str,
    short_name: str,
    **extras,
):
    """Read HDF file and build iris cube with auxiliary data for special variables."""
    f = SD(os.path.join(in_dir, filepath), SDC.READ)

    data_obj = f.select(raw_name)
    data = data_obj.get().astype(np.float32)

    # mask fill values
    if hasattr(data_obj, "fill_value"):
        data[data == data_obj.fill_value] = np.nan

    # mask values outside valid range
    if hasattr(data_obj, "valid_range"):
        valid_range = data_obj.valid_range
        logger.debug(f"Masking {short_name} data outside range {valid_range}")
        data[(data < valid_range[0]) | (data > valid_range[1])] = np.nan

    # apply scale factor and add offset if they exist
    if hasattr(data_obj, "scale_factor"):
        scale_factor = data_obj.scale_factor
        logger.debug(f"Scaling {short_name} data by {scale_factor}")
        data *= scale_factor
    if hasattr(data_obj, "add_offset"):
        add_offset = data_obj.add_offset
        logger.debug(f"Offsetting {short_name} data by {add_offset}")
        data += add_offset

    # for certain variables, we need to read additional datasets to compute the final variable

    # liquid water path in-cloud, need to use cirrus fraction and cloud fraction to estimate liquid fraction, then multiply by in-cloud lwp to get grid-box average lwp
    if short_name in ["clwvi"]:
        cfirFm_obj = f.select("Cirrus_Fraction_Infrared_FMean")
        cfin = cfirFm_obj.get().astype(np.float32)
        # unpack with scale_factor and add_offset
        cif = (cfin * cfirFm_obj.scale_factor) + cfirFm_obj.add_offset
        cif[cif > 0.999] = np.nan

        cfmm_obj = f.select("Cloud_Fraction_Mean_Mean")
        cfmm = cfmm_obj.get().astype(np.float32)
        ctot = (cfmm * cfmm_obj.scale_factor) + cfmm_obj.add_offset

        # liquid fraction estimated assuming random overlap:
        # ctot = 1 - (1 - cif) * (1 - lif)
        # --> lif = 1 - (1 - ctot) / (1 - cif)
        lif = 1.0 - (1.0 - ctot) / (1.0 - cif)
        lif[lif < 0] = 0

        # ice water path
        cwpimm_obj = f.select("Cloud_Water_Path_Ice_Mean_Mean")
        cwpimm = cwpimm_obj.get().astype(np.float32)
        cwpimm_fv = cwpimm_obj._FillValue
        cwpimm[cwpimm == cwpimm_fv] = np.nan  # mask all fill values

        # convert iwp in-cloud to gridbox avg using masked cirrus fraction
        iwp = cwpimm * cif

        # liquid water path
        # convert lwp in-cloud value to grid-box avg
        lwp = data * lif

        data = lwp + iwp

    # these variables are all in-cloud values, so we need to multiply by cirrus fraction to get grid-box average
    elif short_name in ["clivi", "iwpStderr", "lwpStderr"]:
        # cirrus fraction (0-1)
        cfswirFm_obj = f.select("Cirrus_Fraction_Infrared_FMean")
        cfswirFm = cfswirFm_obj.get().astype(np.float32)
        # unpack with scale_factor and add_offset
        cf = (cfswirFm * cfswirFm_obj.scale_factor) + cfswirFm_obj.add_offset
        cf[cf > 0.999] = np.nan

        # lwpStderr is only for liquid water path, so we want to use the liquid fraction to convert to grid-box average
        if short_name in ["lwpStderr"]:
            cfmm_obj = f.select("Cloud_Fraction_Mean_Mean")
            cfmm = cfmm_obj.get().astype(np.float32)
            # unpack with scale_factor and add_offset
            ctot = (cfmm * cfmm_obj.scale_factor) + cfmm_obj.add_offset

            cif = cf
            cif[cf > 0.999] = np.nan

            cf = 1.0 - (1.0 - ctot) / (1.0 - cif)
            cf[cf < 0] = 0

        data *= cf  # "grid-box average" lwp/iwp

    # mask all null values (nan, inf)
    data = np.ma.masked_invalid(data)

    # Handle units - convert invalid units to '1'
    units_str = data_obj.attributes().get("units", "1")
    if units_str in ["none", "None", ""]:
        units_str = "1"

    # Create time coordinate for the start of the month, with bounds covering the whole month
    time_point = datetime(
        year=year, month=month, day=15
    )  # use 15th as representative time for the month
    time_bounds_lower = datetime(year=year, month=month, day=1)
    time_bounds_upper = datetime(
        year=year + (month == 12), month=month + 1 - (month == 12) * 12, day=1
    ) - timedelta(
        days=1
    )  # end of month is the day before the first day of the next month

    # Convert time to days since 1850-01-01
    time_units = datetime(1850, 1, 1)
    delta_1850 = (time_point - time_units).days
    # After creating time delta, add bounds:
    delta_bounds_lower = (time_bounds_lower - time_units).days
    delta_bounds_upper = (time_bounds_upper - time_units).days

    time_coord = iris.coords.DimCoord(
        points=[delta_1850],
        standard_name="time",
        units="days since 1850-01-01 00:00:00",
        bounds=[[delta_bounds_lower, delta_bounds_upper]],
    )
    data = data[np.newaxis, :, :]  # Add time dimension at start

    lats = f.select("YDim").get().astype(np.float64)
    lat_coord = iris.coords.DimCoord(
        lats,
        standard_name="latitude",
        units="degrees_north",
        bounds=[
            (  # calculate bounds as midpoint between adjacent latitudes, assuming regular grid
                lats[i] - 0.5 * np.abs(lats[1] - lats[0]),
                lats[i] + 0.5 * np.abs(lats[1] - lats[0]),
            )
            for i in range(len(lats))
        ],
    )

    lons = f.select("XDim").get().astype(np.float64)
    lon_coord = iris.coords.DimCoord(
        lons,
        standard_name="longitude",
        units="degrees_east",
        bounds=[
            (  # calculate bounds as midpoint between adjacent longitudes, assuming regular grid
                lons[i] - 0.5 * np.abs(lons[1] - lons[0]),
                lons[i] + 0.5 * np.abs(lons[1] - lons[0]),
            )
            for i in range(len(lons))
        ],
    )

    if short_name == "od550aer":
        # add auxiliary coordinate for wavelength (550 nm for this variable, which is AOD at 550 nm)
        wavelength_coord = iris.coords.AuxCoord(
            [550],
            standard_name="radiation_wavelength",
            var_name="wavelength",
            units="nm",
        )
        aux_coords_and_dims = [
            (
                wavelength_coord,
                None,
            ),  # wavelength is a scalar auxiliary coordinate
        ]
    else:
        aux_coords_and_dims = (
            None  # no auxiliary coordinates for other variables
        )

    cube = iris.cube.Cube(
        data,
        long_name=data_obj.attributes().get("long_name", raw_name),
        var_name=short_name,
        units=units_str,
        dim_coords_and_dims=[
            (time_coord, 0),
            (lat_coord, 1),
            (lon_coord, 2),
        ],
        aux_coords_and_dims=aux_coords_and_dims,
    )

    return cube


def convert(
    cube: iris.cube.Cube,
    cmor_info,
    attrs,
) -> iris.cube.Cube:
    """Convert cube to CMOR standards based on data type and cmor_info."""

    if attrs.get("reference") is None:
        attrs["reference"] = cmor_info.attributes["reference"]

    cube.convert_units(cmor_info.units)

    if cube.coord("longitude").points.min() < 0:
        # convert from [-180, 180] to [0, 360]
        cube.coord("longitude").points = cube.coord("longitude").points + 180.0
        # roll the data as part of the longitude conversion to maintain correct order (assuming regular grid and that longitude is the last dimension)
        nlon = len(cube.coord("longitude").points)
        dim_lon = 2
        cube.data = da.roll(cube.core_data(), int(nlon / 2), axis=dim_lon)

    if np.diff(cube.coord("latitude").points)[0] < 0:
        # convert [90,-90] to [-90,90]
        cube.coord("latitude").points = cube.coord("latitude").points[::-1]
        # flip the data
        cube.data = cube.data[:, ::-1, :]  # latitude is axis=1

    utils.set_global_atts(cube, attrs)
    utils.fix_var_metadata(cube, cmor_info)
    utils.fix_dim_coordnames(cube)

    utils.fix_coords(cube)

    return cube


def _extract_variable(
    filepaths: list,
    var_d,
    attrs,
    cmor_table,
    cfg,
    in_dir: str,
    out_dir: str,
    year: int,
) -> None:
    """Extract variable data from input files, convert to CMOR standards, and save."""

    # add mip to attributes for use in cmor table lookup and global attributes
    attrs["mip"] = var_d["mip"]

    # get the cmor table for the variable
    cmor_info = cmor_table.get_variable(
        var_d.get("mip"), var_d.get("short_name")
    )

    # if this variable has a reference specified in the config, use that, otherwise use the one from the cmor table
    if var_d.get("reference"):
        attrs["reference"] = var_d.get("reference")

    if not cmor_table:
        raise Exception(
            f'project_id "{var_d.get("project_id", "None")}" is invalid, no CMOR table found. Valid options: {", ".join(list(CMOR_TABLES.keys()))}'
        )

    # check if cmor info was found for the variable, if not log an error and skip
    if cmor_info is None:
        logger.error(
            f"CMOR info not found for variable {var_d.get('short_name')} in mip {var_d.get('mip')}"
        )
        return
    logger.debug(
        f"CMOR info for variable {var_d.get('short_name')}: {cmor_info}"
    )

    logger.info(
        f"CMORizing variable {var_d.get('short_name')} from input set {var_d.get('raw_name')}"
    )

    cubes = iris.cube.CubeList()
    for filepath in filepaths:
        month = get_month_from_filepath(filepath, cfg)
        logger.debug(
            f"File: {filepath} -> var={var_d.get('short_name')} year={year}, month={month}"
        )
        cube = read_hdf(in_dir, filepath, year, month, **var_d)

        logger.debug(
            f"Read cube for variable {var_d.get('short_name')}: {cube}, with attributes: {cube.attributes}"
        )
        cubes.append(cube)

    if len(cubes) > 0:
        cube = cubes.concatenate()[
            0
        ]  # concatenate along time dimension and take the resulting cube
    else:
        cube = cubes[0]

    cube = convert(cube, cmor_info, attrs)

    # add dataset_id to global attributes for provenance
    cube_attrs = cube.attributes.globals
    cube_attrs["dataset_id"] = attrs["dataset_id"]
    logger.debug(f"Saving variable {cube} with attributes: {cube_attrs}")

    utils.save_variable(
        cube,
        var_d.get("short_name"),
        out_dir,
        cube_attrs,
        unlimited_dimensions=["time"],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Run CMORizer for MODIS."""

    cmor_table = cfg["cmor_table"]
    glob_attrs = cfg["attributes"]

    # get relevant files
    filepaths = collect_files(in_dir, cfg)

    # group by year
    years_d = group_files_by_year(filepaths, cfg)
    years = sorted(list(years_d.keys()))

    # iterate over the variables to be CMORized
    for var, var_d in cfg["variables"].items():
        logger.info(
            "CMORizing var %s from input set %s", var, var_d["raw_name"]
        )
        var_d["short_name"] = var
        logger.debug("\n" + format(var_d))

        # cmorize by year, merging all months at the end
        for year in years:
            _extract_variable(
                years_d[year],
                var_d,
                glob_attrs,
                cmor_table,
                cfg,
                in_dir,
                out_dir,
                year,
            )
