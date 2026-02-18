"""ESMValTool CMORizer for ESACCI-PERMAFROST data.

Tier
   Tier 2: other freely-available dataset.

Source
   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/permafrost/data

Last access
   20240227

Download and processing instructions
   Download the data from:
     active_layer_thickness/L4/area4/pp/v03.0
     ground_temperature/L4/area4/pp/v03.0
     permafrost_extent/L4/area4/pp/v03.0
   Put all files in a single directory.
"""

import logging
import os
import os.path
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import iris
import numpy as np
from cdo import Cdo
from dateutil import relativedelta
from esmvalcore.cmor.table import CMOR_TABLES
from netCDF4 import Dataset

from esmvaltool.cmorizers.data.utilities import save_variable, set_global_atts

logger = logging.getLogger(__name__)


def _fix_coordinates(cube: iris.cube.Cube, definition):
    """Fix coordinates."""
    axis2def = {"T": "time", "X": "longitude", "Y": "latitude", "Z": "sdepth"}
    axes = ["T", "X", "Y", "Z"]

    for axis in axes:
        coord_def = definition.coordinates.get(axis2def[axis])
        if coord_def:
            coord = cube.coord(axis=axis)
            if axis == "T":
                coord.convert_units("days since 1850-1-1 00:00:00.0")
            coord.standard_name = coord_def.standard_name
            coord.var_name = coord_def.out_name
            coord.long_name = coord_def.long_name
            coord.points = coord.core_points().astype("float64")
            if len(coord.points) > 1:
                if coord.bounds is not None:
                    coord.bounds = None
                coord.guess_bounds()

    return cube


def _regrid_infile(infile, outfile, weightsfile):
    """Regrid infile to 0.5 deg x 0.5 deg grid using cdo."""
    cdo = Cdo()
    # ESACCI-PERMAFROST v3.0 dimensions of raw input data
    xsize = 14762
    ysize = 10353
    totalsize = xsize * ysize

    # Description of the ESACCI-PERMAFROST v3.0 polar stereographic
    # grid for CDO.
    # All data below are included as global attibutes and as attributes of
    # variable "polar_stereographic" in the ESACCI-PERMAFROST netCDF files
    # downloaded. A ready-made grid definition file for CDO can also be found
    # in the CDO forum here: https://code.mpimet.mpg.de/boards/2/topics/12769

    esagrid = (
        f"gridtype  = projection\n"
        f"gridsize  = {totalsize}\n"
        f"xsize     = {xsize}\n"
        f"ysize     = {ysize}\n"
        f"xname     = x\n"
        f'xlongname = "x coordinate of projection"\n'
        f'xunits    = "m"\n'
        f"yname     = y\n"
        f'ylongname = "y coordinate of projection"\n'
        f'yunits    = "m"\n'
        f"xfirst    = -6111475.22239475\n"
        f"xinc      = 926.625433138333\n"
        f"yfirst    = 4114895.09469662\n"
        f"yinc      = -926.625433138333\n"
        f"grid_mapping = polar_stereographic\n"
        f"grid_mapping_name = polar_stereographic\n"
        f"straight_vertical_longitude_from_pole = 0.\n"
        f"false_easting = 0.\n"
        f"false_northing = 0.\n"
        f"latitude_of_projection_origin = 90.\n"
        f"standard_parallel = 71.\n"
        f"longitude_of_prime_meridian = 0.\n"
        f"semi_major_axis = 6378137.\n"
        f"inverse_flattening = 298.257223563\n"
        f'proj_params = "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0'
        f' +k=1" +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"\n'
    )

    esagrid_file = "./esacci_grid.txt"

    # write grid description to ASCII file
    with open(esagrid_file, "w", encoding="ascii") as file:
        file.write(esagrid)
        file.close()

    # define dimensions of target grid (regular lat-lon grid)
    target_dimx = 720  # delta_lon = 0.5 deg
    target_dimy = 360  # delta_lat = 0.5 deg
    target_grid = f"r{target_dimx}x{target_dimy}"

    # check if suitable weights file already exists
    # (e.g. from previous call to _regrid_file)

    weightsfile_ok = False

    if os.path.isfile(weightsfile):
        weights = Dataset(weightsfile, "r")
        # make sure dimensions of source and target grids match
        # expected values
        src = weights.variables["src_grid_dims"]
        dst = weights.variables["dst_grid_dims"]
        if (
            xsize == src[0]
            and ysize == src[1]
            and target_dimx == dst[0]
            and target_dimy == dst[1]
        ):
            logger.info(
                "Using matching weights file %s for regridding.",
                weightsfile,
            )
            weightsfile_ok = True
        weights.close()

    # if no suitable weights file, generate new weights for regridding

    if not weightsfile_ok:
        logger.info(
            "Generating regridding weights. This will take"
            " about 5-10 minutes (or more)...",
        )
        # check if path for weight files exists, if not create folder
        path = os.path.split(weightsfile)[0]
        if not os.path.exists(path):
            os.makedirs(path)
        # generate weights
        cdo.genbil(
            f"{target_grid} -setgrid,{esagrid_file}",
            input=infile,
            output=weightsfile,
            options="-f nc",
        )

    # now regrid data to 0.5 deg x 0.5 deg
    cdo.remap(
        f"{target_grid},{weightsfile} -setgrid,{esagrid_file}",
        input=infile,
        output=outfile,
        options="-f nc",
    )


def _extract_variable(in_file, var, cfg, out_dir, year):
    logger.info(
        "CMORizing variable '%s' from input file '%s'",
        var["short_name"],
        in_file,
    )
    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip"]
    attributes["raw"] = var["raw"]
    cmor_table = CMOR_TABLES[attributes["project_id"]]
    definition = cmor_table.get_variable(var["mip"], var["short_name"])

    if "weights_dir" in var.keys():
        weights_dir = var["weights_dir"]
    else:
        weights_dir = "."

    # regrid input file using cdo
    # (using the preprocessor (ESMF) is too slow)

    regridded_file = f"./{year}_{var['short_name']}.nc"
    weights_file = f"{weights_dir}/{year}_{var['short_name']}_weights.nc"
    _regrid_infile(in_file, regridded_file, weights_file)

    # load input file
    cubes = iris.load(regridded_file)

    if len(cubes) > 1:
        # variable gtd contains the vertical levels as separate variables
        # (depth level can only be recognized by the variable names)
        # --> combine all depth levels into 1 cube
        for cube in cubes:
            if cube.var_name == "GST":
                sdepth = 0.0
            elif cube.var_name == "T1m":
                sdepth = 1.0
            elif cube.var_name == "T2m":
                sdepth = 2.0
            elif cube.var_name == "T5m":
                sdepth = 5.0
            elif cube.var_name == "T10m":
                sdepth = 10.0
            else:
                sdepth = 999.0
                logger.info("Could not determin depth. Check results.")
            cube.add_aux_coord(
                iris.coords.AuxCoord(
                    sdepth,
                    standard_name="depth",
                    long_name="depth",
                    units="m",
                ),
            )
            cube.var_name = "gst"
            cube.standard_name = "soil_temperature"  # "valid" standard name
            cube.attributes.pop("actual_min")
            cube.attributes.pop("actual_max")
        tmp_cube = cubes.merge_cube()
        # setting the attribute 'positive' is needed for Iris to recognize
        # this coordinate as 'Z' axis
        tmp_cube.coord("depth").attributes["positive"] = "down"
        # swap coordinates 'depth' and 'time':
        #     (depth, time, lat, lon) --> (time, depth, lat, lon)
        flipped_data = np.swapaxes(tmp_cube.core_data(), 1, 0)
        coord_spec = [
            (tmp_cube.coord("time"), 0),
            (tmp_cube.coord("depth"), 1),
            (tmp_cube.coord("latitude"), 2),
            (tmp_cube.coord("longitude"), 3),
        ]
        cube = iris.cube.Cube(flipped_data, dim_coords_and_dims=coord_spec)
        cube.metadata = tmp_cube.metadata
        # change units string so unit conversion from deg C --> K will work
        cube.units = "celsius"
        # convert units from degC to K
        cube.convert_units("K")
    else:
        cube = cubes[0]

    # --> drop attributes that differ among input files for different years
    # global attributes to remove
    drop_attrs = [
        "source",
        "date_created",
        "history",
        "tracking_id",
        "id",
        "time_coverage_start",
        "time_coverage_end",
        "platform",
        "sensor",
        "keywords",
    ]
    # variable attributes to remove
    drop_var_attrs = [
        "flag_meanings",
        "flag_values",
        "grid_mapping",
        "actual_range",
        "ancillary_variables",
    ]
    for attr in drop_attrs:
        if attr in cube.attributes.keys():
            cube.attributes.pop(attr)
    for attr in drop_var_attrs:
        if attr in cube.attributes.keys():
            cube.attributes.pop(attr)

    set_global_atts(cube, attributes)

    cube.coord("time").points = (
        cube.coord("time").core_points().astype("float64")
    )

    # Set correct names
    cube.var_name = definition.short_name
    cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    # Fix units
    # input variable for pfr reports 'percent' --> rename to '%'
    # input variable for alt reports 'metres' --> rename to 'm'
    # input variable for gtd has been converted to 'K' --> nothing to do
    cube.units = definition.units

    # Fix data type
    cube.data = cube.core_data().astype("float32")

    # Fix coordinates
    cube = _fix_coordinates(cube, definition)

    # Save results
    logger.debug("Saving cube\n%s", cube)
    logger.debug("Setting time dimension to UNLIMITED while saving!")
    save_variable(
        cube,
        cube.var_name,
        out_dir,
        attributes,
        unlimited_dimensions=["time"],
    )
    os.remove(regridded_file)  # delete temporary file
    logger.info("Finished CMORizing %s", in_file)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """CMORize ESACCI-PERMAFROST dataset."""
    glob_attrs = cfg["attributes"]

    logger.info(
        "Starting CMORization for tier%s OBS files: %s",
        glob_attrs["tier"],
        glob_attrs["dataset_id"],
    )
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    logger.info(
        "CMORizing ESACCI-PERMAFROST version %s",
        glob_attrs["version"],
    )

    if start_date is None:
        start_date = datetime(1997, 1, 1)
    if end_date is None:
        end_date = datetime(2019, 12, 31)

    loop_date = start_date
    while loop_date <= end_date:
        for short_name, var in cfg["variables"].items():
            if "short_name" not in var:
                var["short_name"] = short_name
            in_file = list(
                Path(in_dir).glob(var["file"].format(year=loop_date.year))
            )
            if not in_file:
                logger.info(
                    "%d: no data not found for variable %s",
                    loop_date.year,
                    short_name,
                )
            else:
                _extract_variable(
                    str(in_file[0]), var, cfg, out_dir, loop_date.year
                )

        loop_date += relativedelta.relativedelta(years=1)
