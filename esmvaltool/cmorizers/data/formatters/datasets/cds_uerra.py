"""ESMValTool CMORizer for CDS-UERRA version UERRA-HARMONIE.

Tier
   Tier 3
Source
   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-uerra-europe-soil-levels
Last access
   20191104

Download and processing instructions
   - Open in a browser the data source as specified above
   - Put the right ticks:
      - Tick Origin UERRA-HARMONIE
      - Tick Variable 'Volumetric soil moisture'
      - Tick Soil level 1, 2, 3
      - Tick all available years
      - Tick all available months
      - Tick all available days
      - Tick all available timesteps
   - Click 'submit form'
   - According to ESMValTool practice, put them in the right rawobsdir folder

Notes
-----
   - It might be needed to split up the request into smaller chunks

Modification history
   20190821-A_crezee_bas: written.
"""

import glob
import logging
import os
import re
from copy import deepcopy

import cf_units
import iris
import numpy as np
import xarray as xr
import xesmf as xe
from esmvalcore.preprocessor._regrid import _global_stock_cube
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

# Regular expression to parse a "MxN" cell-specification.
_CELL_SPEC = re.compile(
    r"""\A
        \s*(?P<dlon>\d+(\.\d+)?)\s*
        x
        \s*(?P<dlat>\d+(\.\d+)?)\s*
        \Z
     """,
    re.IGNORECASE | re.VERBOSE,
)

# Default fill-value.
_MDI = 1e20

# Stock cube - global grid extents (degrees).
_LAT_MIN = -90.0
_LAT_MAX = 90.0
_LAT_RANGE = _LAT_MAX - _LAT_MIN
_LON_MIN = 0.0
_LON_MAX = 360.0
_LON_RANGE = _LON_MAX - _LON_MIN


def _generate_cube_from_dimcoords(latdata, londata, circular: bool = False):
    """Generate cube from lat/lon points.

    Parameters
    ----------
    latdata : np.ndarray
        List of latitudes.
    londata : np.ndarray
        List of longitudes.
    circular : bool
        Wrap longitudes around the full great circle. Bounds will not be
        generated for circular coordinates.

    Returns
    -------
    :class:`~iris.cube.Cube`
    """
    lats = iris.coords.DimCoord(
        latdata,
        standard_name="latitude",
        units="degrees_north",
        var_name="lat",
        circular=circular,
    )

    lons = iris.coords.DimCoord(
        londata,
        standard_name="longitude",
        units="degrees_east",
        var_name="lon",
        circular=circular,
    )

    if not circular:
        # cannot guess bounds for wrapped coordinates
        lats.guess_bounds()
        lons.guess_bounds()

    # Construct the resultant stock cube, with dummy data.
    shape = (latdata.size, londata.size)
    dummy = np.empty(shape, dtype=np.dtype("int8"))
    coords_spec = [(lats, 0), (lons, 1)]
    cube = iris.cube.Cube(dummy, dim_coords_and_dims=coords_spec)

    return cube


def parse_cell_spec(spec):
    """Parse an MxN cell specification string.

    Parameters
    ----------
    spec: str
        ``MxN`` degree cell-specification for the global grid.

    Returns
    -------
    tuple
        tuple of (float, float) of parsed (lon, lat)

    Raises
    ------
    ValueError
        if the MxN cell specification is malformed.
    ValueError
        invalid longitude and latitude delta in cell specification.
    """
    cell_match = _CELL_SPEC.match(spec)
    if cell_match is None:
        emsg = "Invalid MxN cell specification for grid, got {!r}."
        raise ValueError(emsg.format(spec))

    cell_group = cell_match.groupdict()
    dlon = float(cell_group["dlon"])
    dlat = float(cell_group["dlat"])

    if (np.trunc(_LON_RANGE / dlon) * dlon) != _LON_RANGE:
        emsg = (
            "Invalid longitude delta in MxN cell specification "
            "for grid, got {!r}."
        )
        raise ValueError(emsg.format(dlon))

    if (np.trunc(_LAT_RANGE / dlat) * dlat) != _LAT_RANGE:
        emsg = (
            "Invalid latitude delta in MxN cell specification "
            "for grid, got {!r}."
        )
        raise ValueError(emsg.format(dlat))

    return dlon, dlat


def _cmorize_dataset(in_file, var, cfg, out_dir):
    logger.info(
        "CMORizing variable '%s' from input file '%s'",
        var["short_name"],
        in_file,
    )
    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip"]

    cmor_table = cfg["cmor_table"]
    definition = cmor_table.get_variable(var["mip"], var["short_name"])

    cube = iris.load_cube(
        str(in_file), constraint=NameConstraint(var_name=var["raw"])
    )

    # Time has strange values, so use forecast_reference_time instead
    cube.remove_coord("time")
    cube.coord("forecast_reference_time").rename("time")

    # The following lines are essential before applying
    # the common function fix_coords
    # Convert time calendar from proleptic_gregorian to gregorian
    cube.coord("time").units = cf_units.Unit(
        cube.coord("time").units.origin, "gregorian"
    )

    # Set standard_names for lat and lon
    cube.coord("lat").standard_name = "latitude"
    cube.coord("lon").standard_name = "longitude"

    cube = utils.fix_coords(cube)

    # The above command does not return bounds for longitude
    # so explicitly get them here.
    if not cube.coord("longitude").has_bounds():
        cube.coord("longitude").guess_bounds()

    # Set correct names
    cube.var_name = definition.short_name

    cube.long_name = definition.long_name

    # Convert units if required
    cube.units = "1"
    cube.convert_units(definition.units)

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    logger.info("Saving CMORized cube for variable %s", cube.var_name)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)

    return in_file


def _regrid_dataset(in_dir, var, cfg):
    """Regridding of original files.

    This function regrids each file and write to disk appending 'regrid'
    in front of filename.
    """
    # Match any year here
    filepattern = var["file"].format(year="????", month="??")
    filelist = glob.glob(os.path.join(in_dir, filepattern))
    regridder = None
    for infile in filelist:
        _, infile_tail = os.path.split(infile)
        outfile = os.path.join(cfg["work_dir"], infile_tail)
        targetgrid_ds = xr.DataArray.from_iris(
            _global_stock_cube(cfg["custom"]["regrid"])
        )
        input_ds = xr.open_dataset(infile)
        # Do renaming for consistency of coordinate names
        input_ds = input_ds.rename({"latitude": "lat", "longitude": "lon"})
        # Select uppermoist soil level (index 0)
        input_da = input_ds[var["raw"]].isel(soilLayer=0)
        logger.info("Regridding... ")
        # A workaround to avoid spreading of nan values,
        # related to Github issue
        constantval = 10
        input_da = input_da + constantval
        assert int((input_da == 0.0).sum()) == 0  # Make sure that there
        # are no zero's in the data,
        # since they will be masked out
        shapes = (
            input_ds["lon"].shape,
            input_ds["lat"].shape,
            targetgrid_ds["lon"].shape,
            targetgrid_ds["lat"].shape,
        )
        if regridder is None:
            regridder = xe.Regridder(input_ds, targetgrid_ds, "bilinear")
            ref_shapes = shapes
        else:
            assert shapes == ref_shapes
        da_out = regridder(input_da)
        da_out = da_out.where(da_out != 0.0)
        da_out = da_out - constantval

        # Save it.
        logger.info("Saving: %s", outfile)
        da_out.to_netcdf(outfile)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # run the cmorization
    # Pass on the workdir to the cfg dictionary
    cfg["work_dir"] = cfg_user.work_dir
    # If it doesn't exist, create it
    if not os.path.isdir(cfg["work_dir"]):
        logger.info(
            "Creating working directory for regridding: %s", cfg["work_dir"]
        )
        os.mkdir(cfg["work_dir"])

    for short_name, var in cfg["variables"].items():
        var["short_name"] = short_name
        logger.info("Processing var %s", short_name)

        # Regridding
        logger.info("Start regridding to: %s", cfg["custom"]["regrid"])
        _regrid_dataset(in_dir, var, cfg)
        logger.info("Finished regridding")
        logger.info("Start CMORizing")

        for year in range(1961, 2019):
            # File concatenation
            for month in range(1, 13):
                in_file = os.path.join(
                    cfg["work_dir"],
                    var["file"].format(year=year, month=f"{month:02}"),
                )
                if os.path.isfile(in_file):
                    # Read in the full dataset here from 'workdir'
                    logger.info("Start CMORization of file %s", in_file)
                    _cmorize_dataset(in_file, var, cfg, out_dir)
                    logger.info("Finished processing %s-%s", year, month)
                else:
                    logger.info("No files found for %s-%s", year, month)
                    continue

        logger.info("Finished CMORIZATION")
