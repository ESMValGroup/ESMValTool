"""ESMValTool CMORizer for Eppley-VGPM-MODIS data from Oregon State University.

Tier

Source
   http://orca.science.oregonstate.edu/data/1x2/monthly/eppley.r2018.m.chl.m.sst/hdf

Last access
   20190515

Download and processing instructions
   Download and unpack all the *.tar files under a single directory
   (no subdirectories with years) in ${RAWOBS}/Tier2/Eppley-VGPM-MODIS

Modification history
   20190515-lovato_tomas: written.
"""

import glob
import logging
import os

import iris
import numpy as np
import pandas as pd
import xarray as xr

from esmvaltool.cmorizers.data.utilities import (
    constant_metadata,
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube):
        if var == "intpp":
            cube /= 1000.0 * 12.01 * 86400.0
    return cube


def extract_variable(var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    var = var_info.short_name
    cubes = iris.load(raw_info["file"])
    rawvar = raw_info["name"]

    for cube in cubes:
        if cube.var_name == rawvar:
            fix_var_metadata(cube, var_info)
            cube = fix_coords(cube)
            _fix_data(cube, var)
            set_global_atts(cube, attrs)
            save_variable(
                cube,
                var,
                out_dir,
                attrs,
                local_keys=["coordinates"],
                unlimited_dimensions=["time"],
            )


def merge_data(in_dir, out_dir, raw_info):
    """Merge all data into a single file."""
    data_array = []
    var = raw_info["name"]
    filelist = sorted(glob.glob(in_dir + "/" + raw_info["file"] + "*.hdf"))
    for filename in filelist:
        dataset = (
            xr.open_rasterio(filename)
            .rename({"y": "lat", "x": "lon"})
            .squeeze()
            .drop("band")
        )
        # create coordinates
        dataset = dataset.assign_coords(
            time=pd.to_datetime(filename[-11:-4], format="%Y%j")
        )
        dataset = dataset.expand_dims(dim="time", axis=0)
        spacing = 90.0 / dataset.lat.size
        dataset = dataset.assign_coords(
            lat=np.linspace(-90.0 + spacing, 90.0 - spacing, dataset.lat.size)
        )
        dataset.lat.attrs = {"long_name": "Latitude", "units": "degrees_north"}
        dataset = dataset.assign_coords(
            lon=np.linspace(
                -180.0 + spacing, 180.0 - spacing, dataset.lon.size
            )
        )
        dataset.lon.attrs = {"long_name": "Longitude", "units": "degrees_east"}
        # get current file data
        data_array.append(dataset)
    damerge = xr.concat(data_array, dim="time")

    # need data flip to match coordinates
    damerge.data = np.fliplr(damerge.data)

    # save to file
    dataset = damerge.to_dataset(name=var)
    thekeys = {
        "lat": {"_FillValue": False},
        "lon": {"_FillValue": False},
        "time": {"calendar": "gregorian"},
        var: {"_FillValue": -9999.0},
    }
    filename = os.path.join(out_dir, raw_info["file"] + "_merged.nc")
    dataset.to_netcdf(filename, encoding=thekeys, unlimited_dims="time")

    logger.info("Merged data written to: %s", filename)

    return filename


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg["cmor_table"]
    glob_attrs = cfg["attributes"]

    # run the cmorization
    for var, vals in cfg["variables"].items():
        var_info = cmor_table.get_variable(vals["mip"], var)
        glob_attrs["mip"] = vals["mip"]
        raw_info = {"name": vals["raw"], "file": vals["file"]}

        # merge data
        inpfile = merge_data(in_dir, out_dir, raw_info)

        logger.info("CMORizing var %s from file %s", var, inpfile)
        raw_info["file"] = inpfile
        extract_variable(var_info, raw_info, out_dir, glob_attrs)

    # Remove temporary input file
    os.remove(inpfile)
