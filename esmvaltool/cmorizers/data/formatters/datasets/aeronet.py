"""ESMValTool CMORizer for Aeronet data.

Tier
    Tier 3: restricted dataset.

Source
    https://aeronet.gsfc.nasa.gov/

Last access
    20230613

Download and processing instructions

Issues:

Refs:
"""

import logging
import os.path
import re
from datetime import datetime
from typing import NamedTuple

import cf_units
import dask.array as da
import iris
import iris.coords
import iris.cube
import numpy as np
import pandas as pd
from fsspec.implementations.tar import TarFileSystem
from pys2index import S2PointIndex

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

AERONET_HEADER = "AERONET Version 3;"
LEVEL_HEADER = "Version 3: AOD Level 2.0"
LEVEL_DESCRIPTION = (
    "The following data are automatically cloud cleared and quality assured "
    "with pre-field and post-field calibration applied.")
UNITS_HEADER = (
    "UNITS can be found at,,, https://aeronet.gsfc.nasa.gov/new_web/units.html"
)
DATA_QUALITY_LEVEL = "lev20"

CONTACT_PATTERN = re.compile(
    "Contact: PI=(?P<names>[^;]*); PI Email=(?P<emails>.*)")


def compress_column(df, name):
    compressed = df.pop(name).unique()
    assert len(compressed) == 1
    return compressed[0]


# @dataclass(frozen=True)
class AeronetStation(NamedTuple):
    station_name: str
    latitude: float
    longitude: float
    elevation: float
    contacts: str
    df: pd.DataFrame


class AeronetStations(NamedTuple):
    station_name: list[str]
    latitude: list[float]
    longitude: list[float]
    elevation: list[float]
    contacts: list[str]
    df: list[pd.DataFrame]


def parse_contact(contact):
    match = CONTACT_PATTERN.fullmatch(contact)
    if match is None:
        raise RuntimeError("Could not parse contact line %s", contact)
    names = match.group("names").replace("_", " ").split(" and ")
    emails = match.group("emails").split("_and_")
    mailboxes = ", ".join([
        '"{}" <{}>'.format(name, email) for name, email in zip(names, emails)
    ])
    return mailboxes


def load_file(filesystem, path_like):
    with filesystem.open(path_like, mode="rt", encoding="iso-8859-1") as f:
        aeronet_header = f.readline().strip()
        assert aeronet_header == AERONET_HEADER
        station_name = f.readline().strip()
        level_header = f.readline().strip()
        assert level_header == LEVEL_HEADER
        level_description = f.readline().strip()
        assert level_description == LEVEL_DESCRIPTION
        contact_string = f.readline().strip()
        units_header = f.readline().strip()
        assert units_header == UNITS_HEADER
        df = pd.read_csv(
            f,
            index_col=0,
            na_values=-999.0,
            date_format="%Y-%b",
            parse_dates=[0],
            usecols=lambda x: "AOD_Empty" not in x and "NUM_" not in x,
        )
    contacts = parse_contact(contact_string)
    elevation = compress_column(df, "Elevation(meters)")
    latitude = compress_column(df, "Latitude(degrees)")
    longitude = compress_column(df, "Longitude(degrees)")
    data_quality_level = compress_column(df, "Data_Quality_Level")
    assert data_quality_level == DATA_QUALITY_LEVEL
    station = AeronetStation(
        station_name,
        latitude,
        longitude,
        elevation,
        contacts,
        df,
    )
    return station


def sort_data_columns(columns):
    data_columns = [c for c in columns if "NUM_" not in c]
    # assert len(df.columns) == 3*len(data_columns)
    assert len(columns) == len(data_columns)
    aod_columns = [c for c in data_columns if c.startswith("AOD_")]
    precipitable_water_columns = [
        c for c in data_columns if c == "Precipitable_Water(cm)"
    ]
    angstrom_exponent_columns = [
        c for c in data_columns if "_Angstrom_Exponent" in c
    ]
    assert len(data_columns) == (len(aod_columns) +
                                 len(precipitable_water_columns) +
                                 len(angstrom_exponent_columns))
    return (aod_columns, precipitable_water_columns, angstrom_exponent_columns)


def merge_stations(stations):
    columns = {}
    for name, dtype in (
        ("station_name", str),
        ("latitude", np.float64),
        ("longitude", np.float64),
        ("elevation", np.float64),
        ("contacts", str),
        ("df", object),
    ):
        columns[name] = np.array(
            [getattr(station, name) for station in stations],
            dtype=dtype,
        )
    return AeronetStations(**columns)


def assemble_cube(stations, idx, wavelengths=None):
    min_time = np.array([df.index.min() for df in stations.df]).min()
    max_time = np.array([df.index.max() for df in stations.df]).max()
    date_index = pd.date_range(min_time, max_time, freq="MS")
    dfs = [df.reindex(index=date_index) for df in stations.df]
    all_data_columns = np.unique(
        np.array([df.columns for df in dfs], dtype=str),
        axis=0,
    )
    assert len(all_data_columns) == 1
    data_columns = all_data_columns[0]
    (aod_columns, precipitable_water_columns,
     angstrom_exponent_columns) = sort_data_columns(data_columns)
    if wavelengths is None:
        wavelengths = sorted([int(c[4:-2]) for c in aod_columns])

    aod = da.stack([
        np.stack([df[f"AOD_{wl}nm"].values for wl in wavelengths], axis=-1)
        for df in dfs
    ],
                   axis=-1)[..., idx]
    wavelength_points = np.array(wavelengths, dtype=np.float64)
    wavelength_coord = iris.coords.DimCoord(
        points=wavelength_points,
        standard_name="radiation_wavelength",
        long_name="Wavelength",
        var_name="wl",
        units="nm",
    )
    times = date_index.to_pydatetime()
    time_points = np.array(
        [datetime(year=t.year, month=t.month, day=15) for t in times])
    time_bounds_lower = times
    time_bounds_upper = np.array([
        datetime(year=t.year + (t.month == 12),
                 month=t.month + 1 - (t.month == 12) * 12,
                 day=1) for t in times
    ])
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
    name_coord = iris.coords.AuxCoord(
        points=stations.station_name[idx],
        standard_name="platform_name",
        long_name="Aeronet Station Name",
        var_name="station_name",
    )
    elevation_coord = iris.coords.AuxCoord(
        points=stations.elevation[idx],
        standard_name="height_above_mean_sea_level",
        long_name="Elevation",
        var_name="elev",
        units="m",
    )
    latitude_coord = iris.coords.AuxCoord(
        points=stations.latitude[idx],
        standard_name="latitude",
        long_name="Latitude",
        var_name="lat",
        units="degrees_north",
    )
    longitude_coord = iris.coords.AuxCoord(
        points=stations.longitude[idx],
        standard_name="longitude",
        long_name="Longitude",
        var_name="lon",
        units="degrees_east",
    )
    cube = iris.cube.Cube(
        data=da.ma.masked_array(aod, da.isnan(aod), fill_value=1.e20),
        standard_name=(
            "atmosphere_optical_thickness_due_to_ambient_aerosol_particles"),
        long_name="Aerosol Optical Thickness",
        var_name="aod",
        units="1",
        dim_coords_and_dims=[(time_coord, 0), (wavelength_coord, 1)],
        aux_coords_and_dims=[
            (latitude_coord, 2),
            (longitude_coord, 2),
            (elevation_coord, 2),
            (name_coord, 2),
        ],
    )
    return cube


def build_cube(filesystem, paths, wavelengths=None):
    individual_stations = [
        load_file(filesystem, file_path) for file_path in paths
    ]
    stations = merge_stations(individual_stations)
    latlon_points = np.stack([stations.latitude, stations.longitude], axis=-1)
    index = S2PointIndex(latlon_points)
    cell_ids = index.get_cell_ids()
    idx = np.argsort(cell_ids)
    cube = assemble_cube(stations, idx, wavelengths)
    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    raw_filename = cfg['filename']

    tf = TarFileSystem(f"{in_dir}/{raw_filename}")
    paths = tf.glob("AOD/AOD20/MONTHLY/*.lev20")
    versions = np.unique(
        np.array([os.path.basename(p).split("_")[1] for p in paths],
                 dtype=str))
    assert len(versions) == 1
    version = versions[0]
    wavelengths = sorted(
        [var["wavelength"] for var in cfg['variables'].values()])
    cube = build_cube(tf, paths, wavelengths)

    attrs = cfg['attributes'].copy()
    attrs['version'] = version
    attrs['source'] = attrs['source']

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)

        idx = wavelengths.index(var["wavelength"])
        sub_cube = cube[:, idx]

        attrs['mip'] = var['mip']
        # attrs['reference'] = var['reference']
        # Fix metadata
        utils.set_global_atts(sub_cube, attrs)

        # Save variable
        utils.save_variable(
            sub_cube,
            short_name,
            out_dir,
            attrs,
            unlimited_dimensions=['time'],
        )
