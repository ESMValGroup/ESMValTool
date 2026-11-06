"""ESMValTool CMORizer for MISR data.

Tier
    Tier 3

Source
    ftp://l5ftl01.larc.nasa.gov/MISR/MIL3MAEN.004/
"""

import copy
import datetime as dt
import logging
import os
from pathlib import Path

import numpy as np
import xarray as xr
from dask import array as da
from esmvalcore.cmor.table import CMOR_TABLES

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

band550 = {"name": "green_558nm", "lambda": 558}


def _extract_variable(short_name, var, cfg, in_dir, out_dir):
    attrs = copy.deepcopy(cfg["attributes"])
    attrs["mip"] = var["mip"]
    ver = attrs["version"]
    files = attrs["files"]
    raw_var = var.get("raw_name", short_name)

    cmor_table = CMOR_TABLES[attrs["project_id"]]
    cmor_info = cmor_table.get_variable(var["mip"], short_name)

    """Extract variable."""
    # load data
    logger.debug(
        "Loading data from file(s) '%s' in directory '%s' with version '%s'",
        files,
        in_dir,
        ver,
    )
    for filepath in Path(os.path.join(in_dir)).glob(files):
        xrds = xr.open_dataset(filepath, group="Aerosol_Parameter_Average")
        xrvar = xrds.sel(Band=band550["name"], Optical_Depth_Range="all")[
            raw_var
        ]

        # change order of latitude and longitude coordinates
        # xrvar = xrvar.transpose()

        # Add additional coordinates before converting to an iris cube, as this is easier with xarray

        # Time not present in source data, needs to be added manually
        # Determine time from filename:
        fileparts = str(filepath).split("_")
        year = int(fileparts[-3])
        monthstr = fileparts[-4]
        month = [
            "JAN",
            "FEB",
            "MAR",
            "APR",
            "MAY",
            "JUN",
            "JUL",
            "AUG",
            "SEP",
            "OCT",
            "NOV",
            "DEC",
        ].index(monthstr) + 1
        days_since_1850 = dt.date(year, month, 15) - dt.date(1850, 1, 1)
        lb_since_1850 = dt.date(year, month, 1) - dt.date(1850, 1, 1)
        if month == 12:
            ub_since_1850 = dt.date(year + 1, 1, 1) - dt.date(1850, 1, 1)
        else:
            ub_since_1850 = dt.date(year, month + 1, 1) - dt.date(1850, 1, 1)

        xrvar = xrvar.assign_coords(time=days_since_1850.days)
        xrvar = xrvar.expand_dims("time", axis=2)
        xrvar["time"].attrs["units"] = "days since 1850-01-01 00:00:00"

        if short_name in ["od550aer", "abs550aer"]:
            xrvar = xrvar.assign_coords(radiation_wavelength=band550["lambda"])
            xrvar["radiation_wavelength"].attrs["units"] = "nm"

        cube = xrvar.to_iris()

        # Fix metadata
        cube.coord("Geodetic Latitude").rename("latitude")
        cube.coord("Geodetic Longitude").rename("longitude")

        # add time bounds
        cube.coord("time").bounds = np.array(
            [ub_since_1850.days, lb_since_1850.days]
        )

        utils.fix_var_metadata(cube, cmor_info)
        utils.set_global_atts(cube, attrs)

        utils.fix_dim_coordnames(cube)

        # When Dask tries to roll this cube, it fails because it can't chunk this properly
        # So here we replicate the part of fix_coords that does that, except with numpy.roll
        # instead of dask.roll.
        # However, it seems that those last two lines had no effect on the results...
        # cube.data has shape (360, 720, 1) i.e. (lat, lon, time)
        # so the original code tried to roll the cube on the time axis
        # I could hard-code the lon axis (like they did), but instead I try to autodetect it.
        cube_coord = cube.coord("longitude")
        logger.debug("Fixing longitude...")
        if cube_coord.ndim == 1:
            if cube_coord.points[0] < 0.0 and cube_coord.points[-1] < 181.0:
                cube_coord.points = cube_coord.points + 180.0
                cube.attributes["geospatial_lon_min"] = 0.0
                cube.attributes["geospatial_lon_max"] = 360.0
                nlon = len(cube_coord.points)
                axis = cube.coords().index(cube_coord)
                shift = nlon // 2
                cube.data = da.roll(cube.core_data(), shift, axis=axis)

        if np.diff(cube.coord("latitude").points)[0] < 0:
            # convert [90,-90] to [-90,90]
            cube.coord("latitude").points = cube.coord("latitude").points[::-1]
            # flip the data
            cube.data = cube.data[:, ::-1, :]  # latitude is axis=1

        lat_bounds = []
        for lat in cube.coord("latitude").points:
            lat_bounds.append([lat - 0.25, lat + 0.25])
        lat_bounds = np.array(lat_bounds)
        cube.coord("latitude").bounds = lat_bounds.reshape(-1, 2)

        lon_bounds = []
        for lon in cube.coord("longitude").points:
            lon_bounds.append([lon - 0.25, lon + 0.25])
        lon_bounds = np.array(lon_bounds)
        cube.coord("longitude").bounds = lon_bounds.reshape(-1, 2)

        utils.fix_coords(cube)

        # fix the wavelength coordinate information.
        if short_name in ["od550aer", "abs550aer"]:
            cube.coord("radiation_wavelength").var_name = "wavelength"
            cube.coord("wavelength").standard_name = "radiation_wavelength"

        utils.set_global_atts(cube, attrs)

        # Save variable
        logger.debug(f"Saving Cube: {cube}, in directory: {out_dir}")
        utils.save_variable(
            cube, short_name, out_dir, attrs, unlimited_dimensions=["time"]
        )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Run CMORizer for MISR."""
    cfg.pop("cmor_table")

    for short_name, var in cfg["variables"].items():
        logger.info(f"CMORizing variable '{short_name}'")
        _extract_variable(short_name, var, cfg, in_dir, out_dir)
