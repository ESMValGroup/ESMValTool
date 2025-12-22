"""ESMValTool CMORizer for NOAA GML surface flask data.

Tier
    Tier 2: freely available dataset.

Source
    https://gml.noaa.gov/

Last access
    20240730

Download and processing instructions
    Download one of the following file:
    https://gml.noaa.gov/aftp/data/trace_gases/ch4/flask/surface/ch4_surface-flask_ccgg_text.tar.gz
    https://gml.noaa.gov/aftp/data/trace_gases/co2/flask/surface/co2_surface-flask_ccgg_text.tar.gz
    https://gml.noaa.gov/aftp/data/trace_gases/n2o/flask/surface/n2o_surface-flask_ccgg_text.tar.gz
"""

import logging
import os
from datetime import datetime
from typing import NamedTuple

import cf_units
import dask.array as da
import iris
import iris.coords
import iris.cube
import numpy as np
import pandas as pd
import requests
from fsspec.implementations.tar import TarFileSystem
from pys2index import S2PointIndex

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

FLASK_COLUMNS = ["site", "year", "month", "value"]
DTYPE_FLASK_COLUMNS = {"site": str, "year": int, "month": int, "value": float}
TRACE_GAS_UNITS = {"ch4s": "1e-09", "co2s": "1e-06", "n2os": "1e-09"}


class FlaskStation(NamedTuple):
    """NOAA GML surface flask station data."""

    site_code: str
    site_name: str
    site_country: str
    site_latitude: float
    site_longitude: float
    site_elevation: float
    site_utc2lst: str
    data_frame: pd.DataFrame


class FlaskStations(NamedTuple):
    """NOAA GML surface flask station data."""

    site_code: list[str]
    site_name: list[str]
    site_country: list[str]
    site_latitude: list[float]
    site_longitude: list[float]
    site_elevation: list[float]
    site_utc2lst: list[str]
    data_frame: list[pd.DataFrame]


def get_station_dict():
    """Get station information from online table.

    'Code', 'Name', 'Country', 'Latitude', 'Longitude',
    'Elevation (meters)', 'Time from GMT', 'Project'
    """
    url = "https://www.esrl.noaa.gov/gmd/dv/site/?program=ccgg"
    station_dict = None
    try:
        req = requests.get(url, timeout=30)
    except requests.exceptions.Timeout:
        logger.debug("Request timed out for URL %s", url)
        req = None
    if req is not None:
        stat_list = pd.read_html(req.content)
        stats = stat_list[-1]
        # Remove asterisk from station names (flags inactive stations)
        stats["Code"] = stats["Code"].str.replace("*", "")
        stats.set_index("Code", drop=False, inplace=True)
        station_dict = stats.to_dict(orient="index")
    return station_dict


def load_file(filesystem, filepath, station_dict):
    """Load NOAA GML surface flask station data from the text file."""
    # Determine how many lines to skip in the header
    skiprows = 0
    with filesystem.open(filepath, mode="rt") as file:
        for line in file:
            if line.startswith("#"):
                skiprows = skiprows + 1
    # Read file as CSV
    with filesystem.open(filepath, mode="rt") as file:
        data_frame = pd.read_csv(
            file,
            delimiter=r"[\s]{1,20}",
            skiprows=skiprows,
            header=None,
            names=FLASK_COLUMNS,
            dtype=DTYPE_FLASK_COLUMNS,
            engine="python",
        )
    # Fetch data from station dictionary if available:
    # code, full_name, country, latitude, longitude, elevation, timezone
    site_code = filepath.split("/")[-1].split("_")[1].upper()
    site_name = "N/A"
    site_country = "N/A"
    site_latitude = np.nan
    site_longitude = np.nan
    site_elevation = np.nan
    site_utc2lst = "N/A"
    if isinstance(station_dict, dict):
        if site_code in station_dict.keys():
            site_name = station_dict[site_code]["Name"]
            site_country = station_dict[site_code]["Country"]
            site_latitude = station_dict[site_code]["Latitude"]
            site_longitude = station_dict[site_code]["Longitude"]
            site_elevation = station_dict[site_code]["Elevation (meters)"]
            site_utc2lst = station_dict[site_code]["Time from GMT"]
    # Check if site location is available otherwise return None
    station = None
    if np.any(~np.isnan([site_latitude, site_longitude])):
        # Datetime index
        data_frame.index = pd.to_datetime(
            data_frame["year"].astype(str)
            + "-"
            + data_frame["month"].astype(str),
        )
        # Create FlaskCO2Station object
        station = FlaskStation(
            site_code,
            site_name,
            site_country,
            site_latitude,
            site_longitude,
            site_elevation,
            site_utc2lst,
            data_frame,
        )
    return station


def merge_stations(stations):
    """Collect and merge station data into a FlaskStations instance."""
    columns = {}
    for name, dtype in (
        ("site_code", str),
        ("site_name", str),
        ("site_country", str),
        ("site_latitude", np.float64),
        ("site_longitude", np.float64),
        ("site_elevation", np.float64),
        ("site_utc2lst", str),
        ("data_frame", object),
    ):
        columns[name] = np.array(
            [getattr(station, name) for station in stations],
            dtype=dtype,
        )
    return FlaskStations(**columns)


def assemble_cube(stations, idx, var_attrs):
    """Assemble Iris cube with station data.

    Parameters
    ----------
    stations : FlaskStations
        Station data
    idx : int
        Unique ids of all stations
    var_attrs : dictionnary
        Contains attributes related to the trace gas

    Returns
    -------
    Iris cube
        Iris cube with station data.

    Raises
    ------
    ValueError
        If station data has inconsistent variable names.
    """
    min_time = np.array([df.index.min() for df in stations.data_frame]).min()
    max_time = np.array([df.index.max() for df in stations.data_frame]).max()
    date_index = pd.date_range(min_time, max_time, freq="MS")
    data_frames = [df.reindex(index=date_index) for df in stations.data_frame]
    all_data_columns = np.unique(
        np.array([df.columns for df in data_frames], dtype=str),
        axis=0,
    )
    if len(all_data_columns) != 1:
        raise ValueError(
            "Station data frames has different sets of column names.",
        )

    trace_gas = da.stack([df["value"].values for df in data_frames], axis=-1)[
        ...,
        idx,
    ]

    times = date_index.to_pydatetime()
    time_points = np.array(
        [datetime(year=t.year, month=t.month, day=15) for t in times],
    )
    time_bounds_lower = times
    time_bounds_upper = np.array(
        [
            datetime(
                year=t.year + (t.month == 12),
                month=t.month + 1 - (t.month == 12) * 12,
                day=1,
            )
            for t in times
        ],
    )
    time_bounds = np.stack([time_bounds_lower, time_bounds_upper], axis=-1)
    time_units = cf_units.Unit("days since 1850-01-01", calendar="standard")
    time_coord = iris.coords.DimCoord(
        points=time_units.date2num(time_points),
        standard_name="time",
        long_name="time",
        var_name="time",
        units=time_units,
        bounds=time_units.date2num(time_bounds),
    )
    index_coord = iris.coords.DimCoord(
        points=da.arange(trace_gas.shape[1]),
        standard_name=None,
        long_name="Station index (arbitrary)",
        var_name="station_index",
        units="1",
    )
    code_coord = iris.coords.AuxCoord(
        points=stations.site_code[idx],
        standard_name="platform_name",
        long_name="NOAA GML CCGG Site Name",
        var_name="site_code",
    )
    elevation_coord = iris.coords.AuxCoord(
        points=stations.site_elevation[idx],
        standard_name="height_above_mean_sea_level",
        long_name="Elevation",
        var_name="elev",
        units="m",
    )
    latitude_coord = iris.coords.AuxCoord(
        points=stations.site_latitude[idx],
        standard_name="latitude",
        long_name="Latitude",
        var_name="lat",
        units="degrees_north",
    )
    longitude_coord = iris.coords.AuxCoord(
        points=stations.site_longitude[idx],
        standard_name="longitude",
        long_name="Longitude",
        var_name="lon",
        units="degrees_east",
    )
    cube = iris.cube.Cube(
        data=da.ma.masked_array(
            trace_gas,
            da.isnan(trace_gas),
            fill_value=-999.999,
        ),
        standard_name=(var_attrs["standard_name"]),
        long_name=var_attrs["long_name"],
        var_name=var_attrs["raw_name"],
        units=TRACE_GAS_UNITS[var_attrs["raw_name"]],
        dim_coords_and_dims=[
            (time_coord, 0),
            (index_coord, 1),
        ],
        aux_coords_and_dims=[
            (latitude_coord, 1),
            (longitude_coord, 1),
            (elevation_coord, 1),
            (code_coord, 1),
        ],
    )
    return cube


def build_cube(filesystem, paths, var_attrs):
    """Build station data cube."""
    stations_dict = get_station_dict()
    individual_stations = [
        load_file(filesystem, file_path, stations_dict) for file_path in paths
    ]
    individual_stations = [s for s in individual_stations if s is not None]
    stations = merge_stations(individual_stations)
    latlon_points = np.stack(
        [stations.site_latitude, stations.site_longitude],
        axis=-1,
    )
    index = S2PointIndex(latlon_points)
    cell_ids = index.get_cell_ids()
    idx = np.argsort(cell_ids)
    cube = assemble_cube(stations, idx, var_attrs)
    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    raw_filename = cfg["filename"]

    tar_file_system = TarFileSystem(f"{in_dir}/{raw_filename}")
    paths = tar_file_system.glob(
        f"{cfg['trace_gas']}_surface-flask_ccgg_text/"
        f"{cfg['trace_gas']}_*_month.txt",
    )

    versions = np.unique(
        np.array(
            [os.path.basename(p).split("_")[-3] for p in paths],
            dtype=str,
        ),
    )
    if len(versions) != 1:
        raise ValueError(
            "All station datasets in tar file must have same version.",
        )
    version = versions[0]

    var_attrs = cfg["variables"][f"{cfg['trace_gas']}s"]
    cube = build_cube(tar_file_system, paths, var_attrs)

    attrs = cfg["attributes"].copy()
    attrs["version"] = version
    attrs["source"] = attrs["source"]

    # Run the cmorization
    for short_name, var in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", short_name)

        attrs["mip"] = var["mip"]

        # Fix metadata
        utils.set_global_atts(cube, attrs)

        # Save variable
        utils.save_variable(
            cube,
            short_name,
            out_dir,
            attrs,
            unlimited_dimensions=["time"],
        )
