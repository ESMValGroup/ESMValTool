"""ESMValTool CMORizer for GLORYS12V1 data.

Tier
   Tier 2: other freely-available dataset.

Source
   https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030

Last access
   20260601
"""

import datetime
import logging

import iris
import numpy as np

from esmvaltool.cmorizers.data.utilities import (
    save_variable,
)

logger = logging.getLogger(__name__)


def _fix_latlon_coordinates(cube):
    """Fix latitude/longitude coordinates."""
    for c in ["latitude", "longitude"]:
        cube.coord(c).attributes = None
        cube.coord(c).var_name = c[:3]


def _fix_depth_coordinate(cube):
    """Fix depth coordinate."""
    cube.coord("depth").attributes = {"positive": "down"}
    cube.coord("depth").var_name = "lev"
    cube.coord("depth").long_name = "Ocean depth coordinate"


def _setup_cell_thickness(cubelist, cfg, out_dir):
    """Set up the cell thickness cube of the thkcello variable."""
    logger.info("Creating the cell thickness variable Ofx.thkcello...")
    definition = cfg["cmor_table"].get_variable("Ofx", "thkcello")
    # Extract a 3D variable to get the depth coordinate
    cube_base = cubelist.extract(iris.NameConstraint(var_name="so"))[0]
    # Extract a (depth, lat, lon) cube
    cube_depth = cube_base[0]
    cube_depth.remove_coord("time")
    # Creat cell thickness values from the depth coordinate
    depth_points = cube_depth.coord("depth").points
    thickness_points = np.concatenate(
        [[depth_points[0]], np.diff(depth_points)],
    )
    # Extract a mask from the data following the bathymetry of the original data
    mask = np.isnan(cube_depth.data)
    # Duplicate the thickness coordinate values along the lat, lon axes
    cube_depth.data = np.tile(
        thickness_points[..., np.newaxis, np.newaxis],
        (1, cube_depth.shape[1], cube_depth.shape[2]),
    )
    # Mask the cell thickness following the bathymetry
    cube_depth.data[mask] = np.nan
    # Set up names and attributes
    cube_depth.var_name = definition.short_name
    cube_depth.standard_name = definition.standard_name
    cube_depth.long_name = definition.long_name
    cube_depth.units = definition.units
    cube_depth.attributes["frequency"] = definition.frequency
    cube_depth.attributes["modeling_realm"] = definition.modeling_realm
    _fix_latlon_coordinates(cube_depth)
    _fix_depth_coordinate(cube_depth)
    drop_attrs = ["valid_min", "valid_max", "unit_long"]
    for d in drop_attrs:
        if d in cube_depth.attributes:
            cube_depth.attributes.pop(d)
    attributes = cfg["attributes"].copy()
    attributes["mip"] = "Ofx"
    # Save cube
    save_variable(
        cube_depth,
        cube_depth.var_name,
        out_dir,
        attributes,
    )


def _extract_variable(
    short_name,
    var,
    cubelist,
    start_date,
    end_date,
    cfg,
    out_dir,
):
    """Set up the cube for the variable."""
    logger.info("Cmorising variable %s.%s...", var["mip"], short_name)
    variable = var.get("raw_name", short_name)
    cube_var = cubelist.extract(iris.NameConstraint(var_name=variable))
    # Equalise attributes
    iris.util.equalise_attributes(cube_var)
    # Concatenate cube
    cube_var = cube_var.concatenate_cube()
    # Extract time range if necesary
    time_constraint = iris.Constraint(
        time=lambda t: (t >= start_date) and (t <= end_date),
    )
    cube_var = cube_var.extract(time_constraint)
    # Simple fixes
    definition = cfg["cmor_table"].get_variable(
        var["mip"],
        short_name,
    )
    cube_var.convert_units(definition.units)
    cube_var.attributes.pop("unit_long")
    cube_var.var_name = short_name
    _fix_latlon_coordinates(cube_var)
    _fix_depth_coordinate(cube_var)
    attributes = cfg["attributes"].copy()
    attributes["mip"] = var["mip"]
    # Save cube
    save_variable(
        cube_var,
        cube_var.var_name,
        out_dir,
        attributes,
        unlimited_dimensions=["time"],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    if start_date is None:
        start_date = datetime.datetime(1993, 1, 1, tzinfo=datetime.UTC)
    if end_date is None:
        end_date = datetime.datetime(2026, 4, 30, tzinfo=datetime.UTC)

    # Load input in CubeList object
    cubelist = iris.load(in_dir / "*.nc")

    for short_name, var in cfg["variables"].items():
        # If variable is Ofx.thkcello
        if short_name == "thkcello":
            _setup_cell_thickness(cubelist=cubelist, cfg=cfg, out_dir=out_dir)
        # Other variables from Omon or SImon
        else:
            # Extract variable
            _extract_variable(
                short_name=short_name,
                var=var,
                cubelist=cubelist,
                start_date=start_date,
                end_date=end_date,
                cfg=cfg,
                out_dir=out_dir,
            )
