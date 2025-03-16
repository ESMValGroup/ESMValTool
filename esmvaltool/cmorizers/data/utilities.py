"""Utils module for Python cmorizers."""

import datetime
import gzip
import json
import logging
import os
import re
import shutil
import uuid
from collections.abc import Callable
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path

import esmvalcore.cmor
import iris
import numpy as np
import yaml
from cf_units import Unit
from dask import array as da
from esmvalcore.cmor.check import CheckLevels, CMORCheckError, cmor_check
from esmvalcore.cmor.table import CMOR_TABLES
from esmvalcore.config import CFG
from iris.cube import Cube

import esmvaltool
from esmvaltool import __file__ as esmvaltool_file
from esmvaltool import __version__ as version

logger = logging.getLogger(__name__)

REFERENCES_PATH = Path(esmvaltool_file).absolute().parent / "references"


def add_height2m(cube: Cube) -> None:
    """Add scalar coordinate 'height' with value of 2m to cube in-place.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube which will get the 2m-height coordinate in-place.

    """
    add_scalar_height_coord(cube, height=2.0)


def add_height10m(cube: Cube) -> None:
    """Add scalar coordinate 'height' with value of 10m to cube in-place.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube which will get the 10m-height coordinate in-place.

    """
    add_scalar_height_coord(cube, height=10.0)


def add_scalar_depth_coord(cube: Cube, depth: float = 0.0) -> None:
    """Add scalar coordinate 'depth' to cube in-place.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube which will get the depth coordinate in-place.
    depth: float, optional (default: 0.0)
        Value for the depth in meters.

    """
    logger.debug("Adding depth coordinate (%sm)", depth)
    depth_coord = iris.coords.AuxCoord(
        depth,
        var_name="depth",
        standard_name="depth",
        long_name="depth",
        units=Unit("m"),
        attributes={"positive": "down"},
    )
    try:
        cube.coord("depth")
    except iris.exceptions.CoordinateNotFoundError:
        cube.add_aux_coord(depth_coord, ())
    return cube


def add_scalar_height_coord(cube: Cube, height: float = 2.0) -> None:
    """Add scalar coordinate 'height' to cube in-place.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube which will get the height coordinate in-place.
    height: float, optional (default: 2.0)
        Value for the height in meters.

    """
    logger.debug("Adding height coordinate (%sm)", height)
    height_coord = iris.coords.AuxCoord(
        height,
        var_name="height",
        standard_name="height",
        long_name="height",
        units=Unit("m"),
        attributes={"positive": "up"},
    )
    cube.add_aux_coord(height_coord, ())


def add_typebare(cube, value="bare_ground"):
    """Add scalar coordinate 'typebare' with value of `value`."""
    logger.debug("Adding typebare coordinate (%s)", value)
    typebare_coord = iris.coords.AuxCoord(
        value,
        var_name="typebare",
        standard_name="area_type",
        long_name="surface type",
        units=Unit("no unit"),
    )
    try:
        cube.coord("area_type")
    except iris.exceptions.CoordinateNotFoundError:
        cube.add_aux_coord(typebare_coord, ())
    return cube


@contextmanager
def constant_metadata(cube):
    """Do cube math without modifying units, attributes etc.

    Context manager that should be used when operating on a data cube
    that keeps its metadata constant (units, variable names, attributes etc.).
    Use as with any other context managers: `with constant_metadata(cube):`


    Parameters
    ----------
    cube: iris.cube.Cube
        data cube to be operated on, keeping its
        metadata constant.

    Returns
    -------
    iris.cube.Cube
        Returns the iris cube that was operated on.
    """
    metadata = cube.metadata
    yield metadata
    cube.metadata = metadata


def convert_timeunits(cube, start_year):
    """Convert time axis from malformed Year 0.

    Changes time coordinate with CMOR-like units of
    e.g. `months since START_YEAR-01-01`.

    Parameters
    ----------
    cube: iris.cube.Cube
        data cube to have its time coordinate changed.

    start_year: int
        integer start year as origin of time coordinate

    Returns
    -------
    iris.cube.Cube
        Returns the original iris cube with time coordinate reformatted.
    """
    if cube.coord("time").units == "months since 0000-01-01 00:00:00":
        real_unit = f"months since {str(start_year)}-01-01 00:00:00"
    elif cube.coord("time").units == "days since 0000-01-01 00:00:00":
        real_unit = f"days since {str(start_year)}-01-01 00:00:00"
    elif cube.coord("time").units == "days since 1950-1-1":
        real_unit = "days since 1950-1-1 00:00:00"
    else:
        real_unit = cube.coord("time").units
    cube.coord("time").units = real_unit
    return cube


def fix_coords(
    cube,
    overwrite_time_bounds=True,
    overwrite_lon_bounds=True,
    overwrite_lat_bounds=True,
    overwrite_lev_bounds=True,
    overwrite_airpres_bounds=True,
):
    """Fix coordinates to CMOR standards.

    Fixes coordinates eg time to have correct units, bounds etc;
    longitude to be CMOR-compliant 0-360deg; fixes some attributes
    and bounds - the user can avert bounds fixing by using supplied
    arguments; if bounds are None they will be fixed regardless.

    Parameters
    ----------
    cube: iris.cube.Cube
        data cube with coordinates to be fixed.

    overwrite_time_bounds: bool (optional)
        set to False not to overwrite time bounds.

    overwrite_lon_bounds: bool (optional)
        set to False not to overwrite longitude bounds.

    overwrite_lat_bounds: bool (optional)
        set to False not to overwrite latitude bounds.

    overwrite_lev_bounds: bool (optional)
        set to False not to overwrite depth bounds.

    overwrite_airpres_bounds: bool (optional)
        set to False not to overwrite air pressure bounds.

    Returns
    -------
    cube: iris.cube.Cube
        data cube with fixed coordinates.
    """
    # first fix any completely missing coord var names
    fix_dim_coordnames(cube)
    # fix individual coords
    for cube_coord in cube.coords():
        # fix time
        if cube_coord.var_name == "time":
            logger.info("Fixing time...")
            cube.coord("time").convert_units(
                Unit("days since 1950-1-1 00:00:00", calendar="gregorian")
            )
            if overwrite_time_bounds or not cube.coord("time").has_bounds():
                fix_bounds(cube, cube.coord("time"))

        # fix longitude
        if cube_coord.var_name == "lon":
            logger.info("Fixing longitude...")
            if cube_coord.ndim == 1:
                cube = cube.intersection(longitude=(0.0, 360.0))
            if overwrite_lon_bounds or not cube_coord.has_bounds():
                fix_bounds(cube, cube_coord)

        # fix latitude
        if cube_coord.var_name == "lat":
            logger.info("Fixing latitude...")
            if overwrite_lat_bounds or not cube.coord("latitude").has_bounds():
                fix_bounds(cube, cube.coord("latitude"))
            if cube_coord.core_points()[0] > cube_coord.core_points()[-1]:
                cube = iris.util.reverse(cube, cube_coord)

        # fix depth
        if cube_coord.var_name == "lev":
            logger.info("Fixing depth...")
            if overwrite_lev_bounds or not cube.coord("depth").has_bounds():
                fix_bounds(cube, cube.coord("depth"))

        # fix air_pressure
        if cube_coord.var_name == "air_pressure":
            logger.info("Fixing air pressure...")
            if (
                overwrite_airpres_bounds
                or not cube.coord("air_pressure").has_bounds()
            ):
                fix_bounds(cube, cube.coord("air_pressure"))

    # remove CS
    cube.coord("latitude").coord_system = None
    cube.coord("longitude").coord_system = None

    return cube


def fix_var_metadata(cube, var_info):
    """Fix var metadata from CMOR table.

    Sets var_name, long_name, standard_name and units
    in accordance with CMOR standards from specific CMOR table.

    Parameters
    ----------
    cube: iris.cube.Cube
        data cube to have its metadata changed.

    var_info: class
        CMOR table object holding the information to be changed in the cube.
        Attributes like standard_name, var_name, long_name are used to
        set the new metadata in the input cube.

    Returns
    -------
    iris.cube.Cube
        Returns the masked iris cube.
    """
    if var_info.standard_name == "":
        cube.standard_name = None
    else:
        cube.standard_name = var_info.standard_name
    cube.var_name = var_info.short_name
    cube.long_name = var_info.long_name
    set_units(cube, var_info.units)
    return cube


def flip_dim_coord(cube, coord_name):
    """Flip (reverse) dimensional coordinate of cube."""
    logger.info("Flipping dimensional coordinate %s...", coord_name)
    coord = cube.coord(coord_name, dim_coords=True)
    coord_idx = cube.coord_dims(coord)[0]
    coord.points = np.flip(coord.points)
    if coord.bounds is not None:
        coord.bounds = np.flip(coord.bounds, axis=0)
    cube.data = da.flip(cube.core_data(), axis=coord_idx)


def read_cmor_config(dataset):
    """Read the associated dataset-specific config file."""
    reg_path = os.path.join(
        os.path.dirname(__file__), "cmor_config", dataset + ".yml"
    )
    with open(reg_path, encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    attributes = cfg["attributes"]
    if attributes.get("activity_id", "") == "obs4MIPs":
        # Fill in various attributes automatically.
        attributes["project_id"] = "obs4MIPs"
        attributes["tier"] = "1"
        attributes["source_id"] = dataset
        cv = load_controlled_vocabulary("obs4MIPs")
        for key, value in cv["source_id"][dataset].items():
            attributes[key] = value
        attributes["institution"] = cv["institution_id"][
            attributes["institution_id"]
        ]
    elif "comment" not in attributes:
        attributes["comment"] = ""

    cfg["cmor_table"] = CMOR_TABLES[attributes["project_id"]]

    return cfg


# See https://zenodo.org/records/11500474 for the obs4MIPs specification
# See https://github.com/PCMDI/obs4MIPs-cmor-tables for the obs4MIPs CMOR tables

DRS_ATTRIBUTE = "^[a-zA-Z0-9-]+$"
FREE_FORM_ATTRIBUTE = ".*"


@lru_cache
def load_controlled_vocabulary(project: str) -> dict:
    project_config = yaml.safe_load(
        CFG["config_developer_file"].read_text(encoding="utf-8")
    )[project]
    install_dir = os.path.dirname(os.path.realpath(esmvalcore.cmor.__file__))
    cmor_type = project_config.get("cmor_type", "CMIP5")
    default_path = os.path.join(install_dir, "tables", cmor_type.lower())
    tables_path = project_config.get("cmor_path", default_path)
    tables_path = os.path.expandvars(os.path.expanduser(tables_path))
    if not os.path.exists(tables_path):
        tables_path = os.path.join(install_dir, "tables", tables_path)
    cv_paths = list(Path(tables_path, "Tables").glob("*_CV.json"))
    if not cv_paths:
        return {}
    cv_path = cv_paths[0]
    cv = json.loads(cv_path.read_text(encoding="utf-8"))
    return cv["CV"]


def check_with_controlled_vocabulary(
    project: str,
    attribute: str,
    value: str,
    attributes: dict[str, str],
) -> bool:
    cv = load_controlled_vocabulary(project)
    values = cv[attribute]
    if attribute in values:
        values = values[attribute]
    return value in values


def check_institutution(
    project: str,
    attribute: str,
    value: str,
    attributes: dict[str, str],
) -> bool:
    cv = load_controlled_vocabulary(project)
    institution_id = attributes["institution_id"]
    return value == cv["institution_id"][institution_id]


def create_source_checker(
    source_attr: str,
) -> Callable:
    def checker(
        project: str, attribute: str, value: str, attributes: dict[str, str]
    ) -> bool:
        cv = load_controlled_vocabulary(project)
        source_id = attributes["source_id"]
        return value == cv["source_id"][source_id][source_attr]

    return checker


def check_source_label(
    project: str,
    attribute: str,
    value: str,
    attributes: dict[str, str],
) -> bool:
    source_id = attributes["source_id"]
    return source_id.startswith(f"{value}-")


def check_variant_label(
    project: str,
    attribute: str,
    value: str,
    attributes: dict[str, str],
) -> bool:
    institution_id = attributes["institution_id"]
    return (value == institution_id) or re.match(
        f"^{institution_id}-[a-zA-Z0-9-]+$", value
    )


def check_datetime(
    project: str,
    attribute: str,
    value: str,
    attributes: dict[str, str],
) -> bool:
    try:
        datetime.datetime.fromisoformat(value)
    except ValueError as exc:
        logger.error("Invalid datetime format '%s'", exc)
        return False
    return True


REQUIRED_GLOBAL_ATTRIBUTES = {
    "obs4MIPs": {
        "activity_id": "obs4MIPs",
        "contact": FREE_FORM_ATTRIBUTE,
        "creation_date": check_datetime,
        "dataset_contributor": FREE_FORM_ATTRIBUTE,
        "data_specs_version": FREE_FORM_ATTRIBUTE,  # TODO: automate, this should be the GH release of obs4MIPs CMOR tables
        "frequency": check_with_controlled_vocabulary,
        "grid": FREE_FORM_ATTRIBUTE,
        "grid_label": check_with_controlled_vocabulary,
        "institution": check_institutution,
        "institution_id": check_with_controlled_vocabulary,
        "license": FREE_FORM_ATTRIBUTE,
        "nominal_resolution": check_with_controlled_vocabulary,
        "processing_code_location": FREE_FORM_ATTRIBUTE,  # TODO: automate and add check
        "product": check_with_controlled_vocabulary,
        "realm": check_with_controlled_vocabulary,
        "references": FREE_FORM_ATTRIBUTE,
        "region": create_source_checker("region"),
        "source": create_source_checker("source"),
        "source_id": check_with_controlled_vocabulary,
        "source_label": check_source_label,
        "source_type": create_source_checker("source_type"),
        "source_version_number": create_source_checker(
            "source_version_number"
        ),
        "tracking_id": FREE_FORM_ATTRIBUTE,  # TODO: improve check
        "variable_id": DRS_ATTRIBUTE,
        "variant_label": check_variant_label,
    }
}

OPTIONAL_GLOBAL_ATTRIBUTES = {
    "obs4MIPs": {
        "comment": FREE_FORM_ATTRIBUTE,
        "external_variables": check_with_controlled_vocabulary,
        "history": FREE_FORM_ATTRIBUTE,
        "source_data_notes": FREE_FORM_ATTRIBUTE,
        # TODO: Maybe we can add the two attributes below based on info from
        # the automatic download.
        "source_data_retrieval_date": check_datetime,
        "source_data_url": FREE_FORM_ATTRIBUTE,
        "title": FREE_FORM_ATTRIBUTE,
        "variant_info": FREE_FORM_ATTRIBUTE,
    }
}


def check_global_attributes(project: str, attributes: dict[str, str]) -> bool:
    """Check if the required attributes are available for the project."""
    success = True
    # Check that required attributes are present.
    for attr in REQUIRED_GLOBAL_ATTRIBUTES.get(project, {}):
        if attr not in attributes:
            logger.error("Missing global attribute '%s'", attr)
            success = False

    # Check attribute values.
    attr_definitions = REQUIRED_GLOBAL_ATTRIBUTES.get(
        project, {}
    ) | OPTIONAL_GLOBAL_ATTRIBUTES.get(project, {})

    for attr, checker in attr_definitions.items():
        if attr in attributes:
            value = attributes[attr]
            if not isinstance(checker, str):
                result = checker(project, attr, value, attributes)
            else:
                result = re.match(checker, value)

            if not result:
                logger.error(
                    "Invalid value '%s' for attribute '%s', expected a value "
                    "matching '%s'",
                    value,
                    attr,
                    checker,
                )
                success = False

    return success


def _get_attr_from_field_coord(ncfield, coord_name, attr):
    if coord_name is not None:
        attrs = ncfield.cf_group[coord_name].cf_attrs()
        attr_val = [value for (key, value) in attrs if key == attr]
        if attr_val:
            return attr_val[0]
    return None


def _load_callback(raw_cube, field, _):
    """Use this callback to fix anything Iris tries to break."""
    for coord in raw_cube.coords():
        # Iris chooses to change longitude and latitude units to degrees
        # regardless of value in file, so reinstating file value
        if coord.standard_name in ["longitude", "latitude"]:
            units = _get_attr_from_field_coord(field, coord.var_name, "units")
            if units is not None:
                coord.units = units


def _check_formatting(filename: str, attributes: dict[str, str]) -> None:
    """Run final cmorization checks."""
    project = attributes["project_id"]
    logger.info("Checking compliance with '%s' project standards", project)
    cube = iris.load_cube(filename, callback=_load_callback)

    attribute_success = check_global_attributes(
        project, cube.attributes.globals
    )

    try:
        cmor_check(
            cube=cube,
            cmor_table=project,
            mip=attributes["mip"],
            short_name=cube.var_name,
            frequency=cube.attributes.globals.get("frequency"),
            check_level=CheckLevels.STRICT,
        )
    except CMORCheckError as exc:
        logger.error("%s", exc)
        cmor_check_success = False
    else:
        cmor_check_success = True

    success = attribute_success and cmor_check_success
    msg = (
        f"Data in file {filename} is {'' if success else 'not '}"
        f"compliant with '{project}' project standards"
    )
    if success:
        logger.info(msg)
    else:
        raise ValueError(msg)
    # TODO: add concatenate test
    # TODO: add time coverage test


FILENAME_TEMPLATE = {
    "obs4MIPs": "{variable_id}_{frequency}_{source_id}_{variant_label}_{grid_label}",
    "OBS6": "{project_id}_{dataset_id}_{modeling_realm}_{version}_{mip}_{variable_id}",
    "OBS": "{project_id}_{dataset_id}_{modeling_realm}_{version}_{mip}_{variable_id}",
}


def get_output_filename(attrs: dict[str, str], time_range: str | None) -> str:
    """Get the output filename."""
    project = attrs["project_id"]
    filename = FILENAME_TEMPLATE[project].format(**attrs)
    if time_range is not None:
        filename = f"{filename}_{time_range}"
    filename = f"{filename}.nc"
    return filename


def save_variable(cube, var, outdir, attrs, **kwargs):
    """Saver function.

    Saves iris cubes (data variables) in CMOR-standard named files.

    Parameters
    ----------
    cube: iris.cube.Cube
        data cube to be saved.

    var: str
        Variable short_name e.g. ts or tas.

    outdir: str
        root directory where the file will be saved.

    attrs: dict
        dictionary holding cube metadata attributes like
        project_id, version etc.

    **kwargs: kwargs
        Keyword arguments to be passed to `iris.save`
    """
    if var != cube.var_name:
        msg = (
            f"Attempted to save cube with var_name '{cube.var_name}' as "
            f"variable '{var}'"
        )
        raise ValueError(msg)

    # Set global attributes.
    attrs["variable_id"] = cube.var_name
    set_global_atts(cube, attrs)

    # Ensure correct dtypes.
    fix_dtype(cube)

    # Determine the output filename.
    try:
        time = cube.coord("time")
    except iris.exceptions.CoordinateNotFoundError:
        time_suffix = None
    else:
        if (
            len(time.points) == 1
            and "mon" not in cube.attributes.get("mip", "")
        ) or cube.attributes.get("frequency") == "yr":
            year = str(time.cell(0).point.year)
            time_suffix = "-".join([year + "01", year + "12"])
        else:
            date1 = (
                f"{time.cell(0).point.year:d}{time.cell(0).point.month:02d}"
            )
            date2 = (
                f"{time.cell(-1).point.year:d}{time.cell(-1).point.month:02d}"
            )
            time_suffix = "-".join([date1, date2])

    file_name = get_output_filename(attrs, time_suffix)
    file_path = os.path.join(outdir, file_name)
    logger.info("Saving: %s", file_path)

    # Save the cube.
    status = "lazy" if cube.has_lazy_data() else "realized"
    logger.info("Cube has %s data [lazy is preferred]", status)
    iris.save(cube, file_path, fill_value=1e20, **kwargs)

    # Check that the cube complies with the CMOR tables for the project.
    _check_formatting(file_path, attrs)


def extract_doi_value(tags):
    """Extract doi(s) from a bibtex entry."""
    reference_doi = []
    pattern = r"doi\s*=\s*{([^}]+)}"

    if not isinstance(tags, list):
        tags = [tags]

    for tag in tags:
        bibtex_file = REFERENCES_PATH / f"{tag}.bibtex"
        if bibtex_file.is_file():
            reference_entry = bibtex_file.read_text()
            dois = re.findall(pattern, reference_entry)
            if dois:
                for doi in dois:
                    reference_doi.append(f"doi:{doi}")
            else:
                reference_doi.append("doi not found")
                logger.warning(
                    "The reference file %s does not have a doi.", bibtex_file
                )
        else:
            reference_doi.append("doi not found")
            logger.warning(
                "The reference file %s does not exist.", bibtex_file
            )
    return ", ".join(reference_doi)


def _get_processing_code_location() -> str:
    # TODO: make sure current working dir is not dirty and replace version
    # by commit that is available online
    version = ".".join(esmvaltool.__version__.split(".", 3)[:3])
    return f"https://github.com/ESMValGroup/ESMValTool/tree/{version}"


def set_global_atts(cube, attrs):
    """Complete the cmorized file with global metadata."""
    logger.debug("Setting global metadata...")
    attrs = dict(attrs)
    cube.attributes.clear()
    timestamp = datetime.datetime.now(datetime.timezone.utc)
    timestamp_format = "%Y-%m-%dT%H:%M:%SZ"
    now_time = timestamp.strftime(timestamp_format)

    # Necessary attributes
    if attrs["project_id"] == "obs4MIPs":
        glob_dict = {
            "creation_date": now_time,
            "tracking_id": f"hdl:21.14102/{uuid.uuid4()}",
            "processing_code_location": _get_processing_code_location(),
            "variable_id": cube.var_name,
        }
        required_keys = set(REQUIRED_GLOBAL_ATTRIBUTES["obs4MIPs"])
        for key in required_keys | set(OPTIONAL_GLOBAL_ATTRIBUTES["obs4MIPs"]):
            if key in attrs:
                glob_dict[key] = attrs[key]
        missing = required_keys - set(glob_dict)
        if missing:
            msg = (
                "The following required keys are missing from the "
                f"configuration file: {', '.join(sorted(missing))}"
            )
            raise KeyError(msg)
    else:
        try:
            glob_dict = {
                "title": (
                    f"{attrs.pop('dataset_id')} data reformatted for "
                    f"ESMValTool v{version}"
                ),
                "version": attrs.pop("version"),
                "tier": str(attrs.pop("tier")),
                "source": attrs.pop("source"),
                "reference": extract_doi_value(attrs.pop("reference")),
                "comment": attrs.pop("comment"),
                "user": os.environ.get("USER", "unknown user"),
                "host": os.environ.get("HOSTNAME", "unknown host"),
                "history": f"Created on {now_time}",
                "project_id": attrs.pop("project_id"),
            }
        except KeyError as original_error:
            msg = (
                "All CMORized datasets need the global attributes "
                "'dataset_id', 'version', 'tier', 'source', 'reference', "
                "'comment' and 'project_id' "
                "specified in the configuration file"
            )
            raise KeyError(msg) from original_error
        # Additional attributes
        glob_dict.update(attrs)

    cube.attributes.globals = glob_dict


def fix_bounds(cube, dim_coord):
    """Reset and fix all bounds."""
    if len(cube.coord(dim_coord).points) > 1:
        if cube.coord(dim_coord).has_bounds():
            cube.coord(dim_coord).bounds = None
        cube.coord(dim_coord).guess_bounds()

    if cube.coord(dim_coord).has_bounds():
        cube.coord(dim_coord).bounds = da.array(
            cube.coord(dim_coord).core_bounds(), dtype="float64"
        )
    return cube


def fix_dim_coordnames(cube):
    """Perform a check on dim coordinate names."""
    # first check for CMOR standard coord;
    for coord in cube.coords():
        # guess the CMOR-standard x, y, z and t axes if not there
        coord_type = iris.util.guess_coord_axis(coord)
        try:
            coord = cube.coord(axis=coord_type)
        except iris.exceptions.CoordinateNotFoundError:
            logger.warning(
                "Multiple coordinates for axis %s. "
                "This may be an error, specially for regular grids",
                coord_type,
            )
            continue

        if coord_type == "T":
            coord.var_name = "time"
            coord.attributes = {}

        if coord_type == "X":
            coord.var_name = "lon"
            coord.standard_name = "longitude"
            coord.long_name = "longitude coordinate"
            coord.units = Unit("degrees")
            coord.attributes = {}

        if coord_type == "Y":
            coord.var_name = "lat"
            coord.standard_name = "latitude"
            coord.long_name = "latitude coordinate"
            coord.units = Unit("degrees")
            coord.attributes = {}

        if coord_type == "Z":
            if coord.var_name == "depth":
                coord.standard_name = "depth"
                coord.long_name = "ocean depth coordinate"
                coord.var_name = "lev"
                coord.attributes["positive"] = "down"
            if coord.var_name == "pressure":
                coord.standard_name = "air_pressure"
                coord.long_name = "pressure"
                coord.var_name = "air_pressure"
                coord.attributes["positive"] = "up"
    return cube


def fix_dtype(cube):
    """Fix `dtype` of a cube and its coordinates."""
    if cube.dtype != np.float32:
        logger.info(
            "Converting data type of data from '%s' to 'float32'", cube.dtype
        )
        cube.data = cube.core_data().astype(np.float32, casting="same_kind")
    for coord in cube.coords():
        if coord.dtype.kind != "U" and coord.dtype != np.float64:
            logger.info(
                "Converting data type of coordinate points of '%s' from '%s' "
                "to 'float64'",
                coord.name(),
                coord.dtype,
            )
            coord.points = coord.core_points().astype(
                np.float64, casting="same_kind"
            )
        if coord.has_bounds() and coord.bounds_dtype != np.float64:
            logger.info(
                "Converting data type of coordinate bounds of '%s' from '%s' "
                "to 'float64'",
                coord.name(),
                coord.bounds_dtype,
            )
            coord.bounds = coord.core_bounds().astype(
                np.float64, casting="same_kind"
            )


def roll_cube_data(cube, shift, axis):
    """Roll a cube data on specified axis."""
    cube.data = da.roll(cube.core_data(), shift, axis=axis)
    return cube


def set_units(cube, units):
    """Set units in compliance with cf_unit."""
    special = {"psu": 1, "Sv": "1e6 m3 s-1"}
    if units in special:
        cube.units = special[units]
    else:
        cube.units = Unit(units)
    return cube


def unpack_files_in_folder(folder):
    """Unpack all compressed and tarred files in a given folder.

    This function flattens the folder hierarchy, both outside
    and inside the given folder. It also unpack nested files

    Parameters
    ----------
    folder : str
        Path to the folder to unpack
    """
    decompress = True
    while decompress:
        decompress = False
        files = os.listdir(folder)
        files.sort()
        for filename in files:
            full_path = os.path.join(folder, filename)
            if os.path.isdir(full_path):
                logger.info("Moving files from folder %s", filename)
                folder_files = os.listdir(full_path)
                for file_path in folder_files:
                    shutil.move(os.path.join(full_path, file_path), folder)
                os.rmdir(full_path)
                decompress = True
                continue
            if filename.startswith("."):
                continue
            if not filename.endswith((".gz", ".tgz", ".tar")):
                continue
            logger.info("Unpacking %s", filename)
            shutil.unpack_archive(full_path, folder)
            os.remove(full_path)
            decompress = True


def _gunzip(file_name, work_dir):
    filename = os.path.split(file_name)[-1]
    filename = re.sub(r"\.gz$", "", filename, flags=re.IGNORECASE)

    with gzip.open(file_name, "rb") as f_in:
        with open(os.path.join(work_dir, filename), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


try:
    shutil.register_unpack_format(
        "gz",
        [
            ".gz",
        ],
        _gunzip,
    )
except shutil.RegistryError:
    logger.debug("Format gz already registered. Skipping...")
