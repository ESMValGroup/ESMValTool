"""ESMValTool CMORizer for GLWD data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://figshare.com/articles/dataset/Global_Lakes_and_Wetlands_Database_GLWD_version_2_0/28519994

Last access
    20250701

Download and processing instructions
    Download the file GLWD_v2_0_combined_classes_tif.zip

"""

import logging
import shutil
import zipfile
from datetime import datetime
from pathlib import Path

import iris
import numpy as np
from cf_units import Unit
from iris.coords import AuxCoord, CellMethod, DimCoord
from osgeo import gdal

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


TIME_UNITS = Unit("days since 1950-01-01 00:00:00", calendar="standard")


def _create_time_coord():
    """Create time coordinate."""
    # Time bounds of the climatology are set to 1984-2020 following the
    # corresponding publication: https://doi.org/10.5194/essd-17-2277-2025
    time_points = TIME_UNITS.date2num([datetime(2002, 7, 2)])
    time_bounds = [datetime(1984, 1, 1), datetime(2020, 12, 31)]
    time_bounds = TIME_UNITS.date2num(time_bounds)
    # Add new time coordinate to cube
    return DimCoord(
        time_points,
        bounds=time_bounds,
        standard_name="time",
        long_name="time",
        var_name="time",
        units=TIME_UNITS,
        climatological=True,
    )


def _get_bounds(points, resolution):
    """Compute bounds following points and resolution."""
    lower = points - resolution / 2
    upper = points + resolution / 2
    return np.stack([lower, upper], axis=1)


def _create_lat_lon_coords(n_lat, n_lon):
    """Create latitude/longitude coordinates."""
    # The product is covering the area from 180° West to 180° East
    # and from 56° South to 84° North with a resolution of 15 arc-second.
    lat_start = -56
    lon_start = -180
    res = 1 / 240  # 15 / 3600
    # Create coordinate points
    lat_points = lat_start + res * (0.5 + np.arange(n_lat))
    lon_points = lon_start + res * (0.5 + np.arange(n_lon))
    # Define bounds for lat/lon coordinates
    lat_bounds = _get_bounds(lat_points, res)
    lon_bounds = _get_bounds(lon_points, res)
    # Define coordinates
    latitude = DimCoord(
        lat_points,
        bounds=lat_bounds,
        standard_name="latitude",
        var_name="lat",
        long_name="Longitude",
        units="degrees",
    )
    longitude = DimCoord(
        lon_points,
        bounds=lon_bounds,
        standard_name="longitude",
        var_name="lon",
        long_name="Longitude",
        units="degrees",
    )
    return latitude, longitude


def _create_typewetla_coord():
    """Create wetland type coordinate."""
    typewetla = AuxCoord(
        "wetland",
        var_name="typewetla",
        standard_name="area_type",
        long_name="Wetland",
        units=Unit("no unit"),
    )
    return typewetla


def _extract_variable(var, var_info, cmor_info, attrs, filedir, out_dir, cfg):
    """Extract variable."""
    logger.info("Loading input files...")

    # Load data of wetland area
    ds = gdal.Open(Path(filedir) / cfg["area_file"])
    array = ds.ReadAsArray()
    n_lat, n_lon = array.shape

    # Get ocean/fill_value mask from main class array
    # Classes in [0=dry-land, 1,..., 33], fill_value(ocean) = 255
    dl = gdal.Open(Path(filedir) / cfg["main_class_file"])
    main_class = dl.ReadAsArray()
    mask = main_class == 255

    logger.info("Fixing data and creating coordinates...")

    # Fix data:
    #   - mask oceans (fill_value = 255) + set value to 0
    #   - flip latitude axis
    array = np.where(array <= 100, array, 0)
    array = np.ma.array(array, mask=mask)
    array = np.flip(array, axis=0)

    # Time coordinate
    time_coord = _create_time_coord()

    # Latitude and longitude coordinates
    latitude_coord, longitude_coord = _create_lat_lon_coords(n_lat, n_lon)

    # Type wetland coordinate
    typewetla_coord = _create_typewetla_coord()

    # Cube data
    logger.info("Setting up the cube for variable %s", var)
    cube = iris.cube.Cube(
        array,
        standard_name=cmor_info.standard_name,
        units=cmor_info.units,
        dim_coords_and_dims=[(latitude_coord, 0), (longitude_coord, 1)],
    )
    cube.add_aux_coord(time_coord, ())
    cube.add_aux_coord(typewetla_coord, ())

    # Add coordinate time axis of size 1
    cube = iris.util.new_axis(cube, "time")

    # Fix cell methods
    cube.add_cell_method(CellMethod("mean within years", coords=time_coord))
    cube.add_cell_method(CellMethod("mean over years", coords=time_coord))

    # Fix coords
    cube = utils.fix_coords(cube)

    # Fix var metadata
    utils.fix_var_metadata(cube, cmor_info)

    # Fix global metadata
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube, var, out_dir, attrs)


def _unzip(filepath, out_dir):
    """Unzip `*.zip` file."""
    extracted_dir = Path(out_dir) / "tmp_extracted_files"
    logger.info("Starting extraction of %s to %s", filepath, extracted_dir)
    with zipfile.ZipFile(filepath, "r") as zip_ref:
        zip_ref.extractall(extracted_dir)
    logger.info("Succefully extracted file to %s", extracted_dir)
    return extracted_dir


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg["cmor_table"]
    glob_attrs = cfg["attributes"]

    # Run the cmorization
    for var, var_info in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs["mip"] = var_info["mip"]
        cmor_info = cmor_table.get_variable(var_info["mip"], var)

        # Extract file from ZIP archive
        zip_file = Path(in_dir) / cfg["archive_filename"]
        if not Path(zip_file).is_file():
            logger.debug("Skipping '%s', file '%s' not found", var, zip_file)
            continue
        logger.info("Found input file '%s'", zip_file)
        filedir = _unzip(zip_file, out_dir)

        # Extract and save variable file
        _extract_variable(
            var,
            var_info,
            cmor_info,
            glob_attrs,
            filedir,
            out_dir,
            cfg,
        )

        # Remove extracted directory
        shutil.rmtree(Path(filedir))
        logger.info("Removed cached input directory %s", filedir)
