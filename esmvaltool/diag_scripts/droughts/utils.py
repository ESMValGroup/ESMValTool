"""Utility functions for drought diagnostics.

Functions and variables shared by several diagnostics in the drought folder.
Added functions should have a meaningfull name and docstring.
"""

from __future__ import annotations

import datetime as dt
import itertools as it
import logging
from calendar import monthrange
from pathlib import Path
from typing import Any

import cartopy as ct
import iris
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import yaml
from cartopy.mpl.geoaxes import GeoAxes
from cf_units import Unit
from esmvalcore import preprocessor as pp
from iris.cube import Cube, CubeList
from iris.util import equalise_attributes

import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughts.constants import (
    CMIP6_FNAME,
    OBS_FNAME,
)
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_cfg,
    get_diagnostic_filename,
    group_metadata,
    select_metadata,
)
from esmvaltool.diag_scripts.shared._base import _get_input_data_files

log = logging.getLogger(Path(__file__).name)


# GENERAL HELPER


def mkplotdir(cfg: dict, dname: str | Path) -> None:
    """Create a sub directory for plots if it does not exist."""
    new_dir = Path(cfg["plot_dir"] / dname)
    if not new_dir.is_dir():
        Path.mkdir(new_dir)


def quick_save(cube: Cube, name: str, cfg: dict) -> None:
    """Simply save cube to netcdf file without additional information."""
    if cfg.get("write_netcdf", True):
        diag_file = get_diagnostic_filename(name, cfg)
        log.info("quick save %s", diag_file)
        iris.save(cube, target=diag_file)


def quick_load(cfg: dict, context: dict, *, strict=True) -> Cube:
    """Load input files from config wich matches the selection.

    Select, load and return the first match.
    raises an error (if strict) or a warning for multiple matches.
    """
    meta = cfg["input_data"].values()
    var_meta = select_single_meta(meta, strict=strict, **context)
    return iris.load_cube(var_meta["filename"])


def get_plot_fname(
    cfg: dict,
    basename,
    meta: dict | None = None,
    replace: dict | None = None,
) -> str:
    """Get a valid path for saving a diagnostic plot.

    This is an alternative to shared.get_diagnostic_filename.
    It uses cfg as first argument and accept metadata to format the basename.

    Parameters
    ----------
    cfg : dict
        Dictionary with diagnostic configuration.
    basename : str
        The basename of the file.
    meta : dict, optional
        Metadata to format the basename. If None, empty dict is used.
    replace : dict
        Dictionary with strings to replace in the basename.
        If None, empty dict is used.

    Returns
    -------
    str:
        A valid path for saving a diagnostic plot.
    """
    meta = {} if meta is None else meta
    replace = {} if replace is None else replace
    basename = basename.format(**meta)
    for key, value in replace.items():
        basename = basename.replace(key, value)
    fpath = Path(cfg["plot_dir"]) / basename
    return str(fpath.with_suffix(cfg["output_file_type"]))


def add_ancestor_input(cfg: dict) -> None:
    """Read ancestors settings.yml and add it's input_data to this config."""
    log.info("add ancestors for %s", cfg[n.INPUT_FILES])
    for input_file in cfg[n.INPUT_FILES]:
        cfg_anc_file = (
            "/run/".join(input_file.rsplit("/work/", 1)) + "/settings.yml"
        )
        cfg_anc = get_cfg(cfg_anc_file)
        cfg["input_data"].update(_get_input_data_files(cfg_anc))


def abs_auxilary_path(cfg: dict, path: str | Path) -> str:
    """Return absolut path of an auxilary file."""
    if Path(path).is_absolute():
        return str(path)
    return str(Path(cfg["auxiliary_data_dir"]) / path)


def save_metadata(cfg: dict, metadata: dict) -> None:
    """Save dict as metadata.yml in work folder."""
    with (Path(cfg["work_dir"]) / "metadata.yml").open("w") as wom:
        yaml.dump(metadata, wom)


def fix_interval(interval: dict) -> dict:
    """Ensure that an interval has a label and a range.

    TODO: replace "/" with "_" in diagnostics who use this.
    """
    if "range" not in interval:
        interval["range"] = f"{interval['start']}/{interval['end']}"
    if "label" not in interval:
        interval["label"] = interval["range"]
    return interval


def log_provenance(cfg: dict, fname: str | Path, record: dict) -> None:
    """Add provenance information to the Provenancelog."""
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(fname, record)


# DATA PROCESSING


def merge_list_cube(
    cube_list: list,
    aux_name: str = "dataset",
    *,
    equalize: bool = True,
) -> Cube:
    """Merge a list of cubes into a single one with an auxiliary variable.

    Useful for applying statistics along multiple cubes. The time coordinate is
    removed by this function.

    Parameters
    ----------
    cube_list: list
        List or iterable of cubes with the same coordinates.
    aux_name: str, optional
        Name of the new auxiliary coordinate. Defaults to "dataset".
    equalize: bool, optional
        Drops differences in attributes, otherwise raises an error.
        Defaults to True.

    Returns
    -------
    iris.cube.Cube
        Merged cube with an enumerated auxiliary variable.
    """
    for ds_index, ds_cube in enumerate(cube_list):
        coord = iris.coords.AuxCoord(ds_index, long_name=aux_name)
        ds_cube.add_aux_coord(coord)
    cubes = CubeList(cube_list)
    log.info("formed %s: %s", type(cubes), cubes)
    if equalize:
        removed = equalise_attributes(cubes)
        log.info("removed different attributes: %s", removed)
    for cube in cubes:
        cube.remove_coord("time")
    return cubes.merge_cube()


def sort_cube(cube: Cube, coord: str = "longitude") -> Cube:
    """Sort data along a one-dimensional numerical coordinate.

    Parameters
    ----------
    cube : iris.cube.Cube
        The iris cube that should be sorted.
    coord : str
        The name of the one-dimensional coordinate to sort the cube by.

    Returns
    -------
    iris.cube.Cube
        The sorted cube.

    Source
    ------
    https://gist.github.com/pelson/9763057
    """
    coord_to_sort = cube.coord(coord)
    if coord_to_sort.ndim == 1:
        msg = "Only dim coords are supported."
        raise NotImplementedError(msg)
    (dim,) = cube.coord_dims(coord_to_sort)
    index = [slice(None)] * cube.ndim
    index[dim] = np.argsort(coord_to_sort.points)
    return cube[tuple(index)]


def fix_longitude(cube: Cube, coord="longitude") -> Cube:
    """Return a cube with 0 centered longitude coords.

    updating the longitude coord and sorting the data accordingly
    """
    # make sure coords are -180 to 180
    fixed_lons = [
        lon if lon < 180 else lon - 360 for lon in cube.coord(coord).points
    ]
    lon_dim = cube.coord_dims(cube.coord(coord))[0]
    cube.add_aux_coord(
        iris.coords.AuxCoord(fixed_lons, long_name="fixed_lon"),
        lon_dim,
    )
    cube = sort_cube(cube, coord="fixed_lon")
    # set new coordinates as dimcoords, add new dim and remove old and aux
    new_lon = cube.coord("fixed_lon")
    new_lon_dims = cube.coord_dims(new_lon)
    longitude = iris.coords.DimCoord.from_coord(new_lon)
    longitude.rename(coord)
    cube.remove_coord(new_lon)
    cube.remove_coord(cube.coord(coord))
    cube.add_dim_coord(longitude, new_lon_dims)
    return cube


def smooth(dat, window=32, mode="same") -> np.ndarray:
    """Smooth a 1D array with a window size.

    from scipy.ndimage.filters import uniform_filter1d as unifilter
    TODO: iris can also directly filter on cubes:
    https://scitools-iris.readthedocs.io/en/latest/generated/gallery/general/
    plot_SOI_filtering.html#sphx-glr-generated-gallery-general-plot-soi-
    filtering-py
    smoothed = unifilter(dat, 32)  # smooth over 4 year window
    """
    # using numpy convol
    np_filter = np.ones(window)
    return np.convolve(dat, np_filter, mode) / window


def date_to_months(date: str, start_year: int) -> int:
    """Translate date YYYY-MM to number of months since start_year."""
    years, months = [int(x) for x in date.split("-")]
    return int(12 * (years - start_year) + months)


def remove_attributes(
    cube: Cube | iris.Coord,
    ignore: list | None = None,
) -> None:
    """Remove most attributes of cube or coords in place.

    Used to clean up differences in cubes coordinates before merging

    Parameters
    ----------
    cube
        iris.Cube or iris.Coord
    ignore
        Optional: List of Strings of attributes, that are not removed
        By default: []
    """
    if ignore is None:
        ignore = []
    remove = [attr for attr in cube.attributes if attr not in ignore]
    for attr in remove:
        del cube.attributes[attr]


def get_time_range(cube: Cube) -> dict:
    """Guess the period of a cube based on the time coordinate."""
    if not isinstance(cube, Cube):
        cube = iris.load_cube(cube)
    time = cube.coord("time")
    start = time.units.num2date(time.points[0])
    end = time.units.num2date(time.points[-1])
    return {"start_year": start.year, "end_year": end.year}


def monthly_to_daily(
    cube: Cube,
    units: str = "mm day-1",
    *,
    leap_years: bool = True,
) -> None:
    """Convert monthly data to daily data inplace ignoring leap years.

    With leap_years=False this is similar to the same named function in utils.R
    and compatible with the pet.R diagnostic.
    """
    months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    months = months * int((cube.shape[0] / 12) + 1)
    for idx, sli in enumerate(cube.slices_over(["time"])):
        if not leap_years:
            days = months[idx]
            cube.data[idx] = cube.data[idx] / days
            continue
        # consider leap days
        time = sli.coord("time")
        date = time.units.num2date(time.points[0])
        days = monthrange(date.year, date.month)[1]
        cube.data[idx] = cube.data[idx] / days
    cube.units = units


def daily_to_monthly(
    cube: Cube,
    units: str = "mm month-1",
    *,
    leap_years: bool = True,
) -> None:
    """Convert daily data to monthly data inplace.

    With leap_years=False this is similar to the same named function in utils.R
    and compatible with the pet.R diagnostic.
    """
    months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    months = months * int((cube.shape[0] / 12) + 1)
    for idx, sli in enumerate(cube.slices_over(["time"])):
        if not leap_years:
            days = months[idx]
            cube.data[idx] = cube.data[idx] * days
            continue
        # consider leap days
        time = sli.coord("time")
        date = time.units.num2date(time.points[0])
        days = monthrange(date.year, date.month)[1]
        cube.data[idx] = cube.data[idx] * days
    cube.units = units


def _get_data_hlp(axis, data, ilat, ilon):
    """Get data_help dependend on axis."""
    if axis == 0:
        data_help = (data[:, ilat, ilon])[:, 0]
    elif axis == 1:
        data_help = (data[ilat, :, ilon])[:, 0]
    elif axis == 2:
        data_help = data[ilat, ilon, :]
    else:
        data_help = None
    return data_help


def create_cube_from_data(var, data_dict) -> Cube:
    """Create cube to prepare plotted data for saving to netCDF.

    Renamed function from cube_to_save_ploted_data

    Parameters
    ----------
    var : np.ndarray
        Data to be plotted.
    data_dict : dict
        Dictionary containing metadata for the data. It must contain: "var",
        "drought_char", "unit", "latitude", and "longitude".
    """
    cube = Cube(
        var,
        var_name=data_dict["var"],
        long_name=data_dict["drought_char"],
        units=data_dict["unit"],
    )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            data_dict["latitude"],
            var_name="lat",
            long_name="latitude",
            units="degrees_north",
        ),
        0,
    )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            data_dict["longitude"],
            var_name="lon",
            long_name="longitude",
            units="degrees_east",
        ),
        1,
    )
    return cube


def cube_to_save_ploted_ts(data_dict: dict) -> Cube:
    """Create cube to prepare plotted time series for saving to netCDF."""
    cube = Cube(
        data_dict["data"],
        var_name=data_dict["var"],
        long_name=data_dict["var"],
        units=data_dict["unit"],
    )
    coord = iris.coords.DimCoord(
        data_dict["time"],
        var_name="time",
        long_name="Time",
        units="month",
    )
    cube.add_dim_coord(coord, 0)
    return cube


def runs_of_ones_array_spei(bits, spei) -> list:
    """Set 1 at beginning ond -1 at the end of events."""
    # make sure all runs of ones are well-bounded
    bounded = np.hstack(([0], bits, [0]))
    # get 1 at run starts and -1 at run ends
    difs = np.diff(bounded)
    (run_starts,) = np.where(difs > 0)
    (run_ends,) = np.where(difs < 0)
    spei_sum = np.full(len(run_starts), 0.5)
    for iii, indexs in enumerate(run_starts):
        spei_sum[iii] = np.sum(spei[indexs : run_ends[iii]])
    return [run_ends - run_starts, spei_sum]


def count_spells(data, threshold, axis) -> np.ndarray:
    """Functions for Iris Aggregator to count spells."""
    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim
    data = data[:, :, 0, :]
    if axis > 2:
        axis = axis - 1
    listshape = []
    inoax = []
    for iii, ishape in enumerate(data.shape):
        if iii != axis:
            listshape.append(ishape)
            inoax.append(iii)
    listshape.append(4)
    return_var = np.zeros(tuple(listshape))
    for ilat in range(listshape[0]):
        for ilon in range(listshape[1]):
            data_help = _get_data_hlp(axis, data, ilat, ilon)
            if data_help.count() == 0:
                return_var[ilat, ilon, 0] = data_help[0]
                return_var[ilat, ilon, 1] = data_help[0]
                return_var[ilat, ilon, 2] = data_help[0]
                return_var[ilat, ilon, 3] = data_help[0]
            else:
                data_hits = data_help < threshold
                [events, spei_sum] = runs_of_ones_array_spei(
                    data_hits,
                    data_help,
                )
                return_var[ilat, ilon, 0] = np.count_nonzero(events)
                return_var[ilat, ilon, 1] = np.mean(events)
                return_var[ilat, ilon, 2] = np.mean(
                    (spei_sum * events)
                    / (np.mean(data_help[data_hits]) * np.mean(events)),
                )
                return_var[ilat, ilon, 3] = np.mean(spei_sum / events)
    return return_var


def get_latlon_index(coords, lim1, lim2) -> np.ndarray:
    """Get index for given values between two limits (1D), e.g. lats, lons."""
    return (
        np.where(
            np.absolute(coords - (lim2 + lim1) / 2.0) <= (lim2 - lim1) / 2.0,
        )
    )[0]


def slice_cube_interval(cube: Cube, interval: list) -> Cube:
    """Return cube slice for given interval.

    which is a list of strings (YYYY-mm) or int (index of cube)
    For 3D cubes time needs to be first dim.
    """
    if isinstance(interval[0], int) and isinstance(interval[1], int):
        return cube[interval[0] : interval[1], :, :]
    dt_start = dt.datetime.strptime(interval[0], "%Y-%m")
    dt_end = dt.datetime.strptime(interval[1], "%Y-%m")
    time = cube.coord("time")
    t_start = time.nearest_neighbour_index(time.units.date2num(dt_start))
    t_end = time.nearest_neighbour_index(time.units.date2num(dt_end))
    return cube[t_start:t_end, :, :]


def find_first(nparr: np.ndarray) -> int:
    """Return index of first nonzero element of numpy array or -1.

    Its faster than looping or using numpy.where.
    nparr requires numerical data without negative zeroes.
    """
    idx = nparr.view(bool).argmax() // nparr.itemsize
    return int(idx if np.arr[idx] else -1)


def fix_calendar(cube: Cube) -> Cube:
    """Convert cubes calendar to gregorian.

    TODO: use pp.regrid_time when available in esmvalcore.
    """
    time = cube.coord("time")
    # if time.units.name == "days since 1850-1-1 00:00:00":
    log.info("renaming unit")
    time.units = Unit("days since 1850-01-01", calendar=time.units.calendar)
    if time.units.calendar == "proleptic_gregorian":
        time.units = Unit(time.units.name, calendar="gregorian")
        log.info("renamed calendar: %s", time.units.calendar)
    if time.units.calendar != "gregorian":
        time.convert_units(Unit("days since 1850-01-01", calendar="gregorian"))
        log.info("converted time to: %s", time.units)
    return cube


def latlon_coords(cube: Cube) -> None:
    """Rename latitude, longitude coords to lat, lon inplace."""
    if "latitude" in cube.coords():
        cube.coord("latitude").rename("lat")
    if "longitude" in cube.coords():
        cube.coord("longitude").rename("lon")


def guess_lat_lon_bounds(cube: Cube) -> None:
    """Guess bounds for latitude and longitude if missing."""
    if not cube.coord("latitude").has_bounds():
        cube.coord("latitude").guess_bounds()
    if not cube.coord("longitude").has_bounds():
        cube.coord("longitude").guess_bounds()


def mmm(
    cube_list: list | CubeList,
    mdtol: float = 0,
    dropcoords: list | None = None,
    *,
    dropmethods=False,
) -> tuple:
    """Calculate mean and stdev along a cube list over all cubes.

    Return two (mean and stdev) of same shape

    Parameters
    ----------
    cube_list : list|CubeList
        List of iris cubes to be merged by mean.
    mdtol : float, optional
        Tolerance for mean calculation, by default 0
    dropcoords : list|None, optional
        Coordinates to be dropped from the cubes. If None, only time will be
        dropped. To keep all coords pass an empty list.
    dropmethods: bool, optional
        Drop cell_methods from the cubes, by default False
    """
    if dropcoords is None:
        dropcoords = ["time"]
    for idx, cube in enumerate(cube_list):
        for coord_name in dropcoords:
            if cube.coords(coord_name):
                cube.remove_coord(coord_name)
        if dropmethods:
            cube.cell_methods = None
        cube.add_aux_coord(iris.coords.AuxCoord(idx, long_name="dataset"))
    cube_list = CubeList(cube_list)
    equalise_attributes(cube_list)
    try:
        merged = cube_list.merge_cube()
    except iris.exceptions.MergeError as err:
        iris.util.describe_diff(cube_list[0], cube_list[1])
        raise iris.exceptions.MergeError from err
    if mdtol > 0:
        log.info("performing MMM with tolerance: %s", mdtol)
    mean = merged.collapsed("dataset", iris.analysis.MEAN, mdtol=mdtol)
    sdev = merged.collapsed("dataset", iris.analysis.STD_DEV)
    return mean, sdev


def regional_stats(cfg, cube, operator="mean") -> dict:
    """Calculate statistic over AR6 IPCC reference regions."""
    _ = cfg  # we might need this in the future. dont tell codacy!
    guess_lat_lon_bounds(cube)
    extracted = pp.extract_shape(cube, "ar6", decomposed=True)
    return pp.area_statistics(extracted, operator)


def transpose_by_names(cube: Cube, names: list) -> None:
    """Transpose a cube by dim-coords or their names."""
    new_dims = [cube.coord_dims(name)[0] for name in names]
    cube.transpose(new_dims)


# META DATA


def fold_meta(
    cfg: dict,
    meta: dict,
    cfg_keys: list | None = None,
    meta_keys: list | None = None,
    variables: list | None = None,
) -> tuple:
    """Create combinations of meta data and data constraints.

    cfg["variables"] overwrites meta["short_names"].

    Parameters
    ----------
    cfg: dict
        Plot specific configuration with cfg_keys on root level.
    meta: list
        Full meta data including ancestor files.
    cfg_keys: list, optional
        Data constraints as config entries used for product.
        Defaults to ["locations", "intervals"].
    meta_keys: list, optional
        Keys for each meta used for product, short_name added automatically.
        Defaults to ["dataset", "exp"].
    variables: list, optional
        Variables to be used. Defaults to None.

    Returns
    -------
    combinations : itertools.product
        All combinations of the metadata and constraints.
    groups : dict
        Grouped metadata.
    meta_keys : list
        List of metadata keys.
    """
    if variables is None:
        variables = cfg.get("variables", ["pdsi", "spi"])
    if cfg_keys is None:
        cfg_keys = ["locations", "intervals"]
    if meta_keys is None:
        meta_keys = ["dataset", "exp"]

    groups = {
        gk: list(
            group_metadata(
                select_metadata(meta, short_name=variables[0]),
                gk,
            ).keys(),
        )
        for gk in meta_keys
    }
    meta_keys.append("short_name")
    groups["short_name"] = variables
    g_map = {"locations": "location", "intervals": "interval"}
    for ckey in cfg_keys:
        try:
            groups[g_map.get(ckey, ckey)] = cfg[ckey]
        except KeyError:
            log.warning("No '%s' found in plot config", ckey)
    combinations = it.product(*groups.values())
    return combinations, groups, meta_keys


def select_meta_from_combi(meta: list, combi: dict, groups: dict) -> tuple:
    """Select one meta data from list (filter valid keys).

    Parameters
    ----------
    meta : list
        List of metadata dictionaries.
    combi : dict
        Dictionary containing the combination of metadata values.
    groups : dict
        Dictionary containing the groups of metadata.

    Returns
    -------
    tuple
        A tuple containing the selected metadata and the configuration
        dictionary.
    """
    this_cfg = dict(zip(groups.keys(), combi, strict=True))
    filter_cfg = clean_meta(this_cfg)  # remove non meta keys
    this_meta = select_metadata(meta, **filter_cfg)[0]
    return this_meta, this_cfg


def list_meta_keys(meta: list, group: dict) -> list:
    """Return a list of all keys found for a group in the meta data."""
    return list(group_metadata(meta, group).keys())


def get_meta_list(meta: dict, group: str, select: dict | None = None) -> list:
    """
    List all entries found for the group key as a list.

    With a given selection, the meta data will be filtered first.

    Parameters
    ----------
    meta : dict
        Full meta data.
    group : str
        Key to search for. Defaults to "alias".
    select : dict, optional
        Dictionary like {'short_name': 'pdsi'} that is passed to a selection.

    Returns
    -------
    list
        Collected values for the group key.
    """
    if select is not None:
        meta = select_metadata(meta, **select)
    return list(group_metadata(meta, group).keys())


def get_datasets(cfg: dict) -> dict:
    """Return a dictionary of datasets and their metadata."""
    metadata = cfg["input_data"].values()
    return group_metadata(metadata, "dataset").keys()


def get_dataset_scenarios(cfg: dict) -> list:
    """Combine datasets and scenarios to a list of pairs of strings."""
    metadata = cfg["input_data"].values()
    input_datasets = group_metadata(metadata, "dataset").keys()
    input_scenarios = group_metadata(metadata, "alias").keys()
    return list(it.product(input_datasets, input_scenarios))


def get_scenarios(meta, **kwargs) -> list:
    """Return a list of alias values for scenario names."""
    selected = select_metadata(meta, **kwargs)
    return list(group_metadata(selected, "alias").keys())


def get_basename(cfg, meta, prefix=None, suffix=None) -> str:
    """Return a formatted basename for a diagnostic file."""
    _ = cfg  # we might need this in the future. dont tell codacy!
    formats = {
        "CMIP6": CMIP6_FNAME,
        "OBS": OBS_FNAME,
    }
    basename = formats[meta["project"]].format(**meta)
    if suffix:
        basename += f"_{suffix}"
    if prefix:
        basename = f"{prefix}_{basename}"
    return basename


def clean_meta(meta) -> dict:
    """Return a copy of meta data with only selected keys.

    Keys are: short_name, dataset, alias, exp
    """
    valid_keys = ["short_name", "dataset", "alias", "exp"]
    return {key: val for key, val in meta.items() if key in valid_keys}


def select_single_metadata(
    meta: list,
    *,
    strict: bool = True,
    **kwargs: dict[str, Any],
) -> dict | None:
    """Filter meta data by arbitrary keys and return one matching result.

    For more/less then one match the first/none is returned or an error is
    raised with strict=True.

    Parameters
    ----------
    meta
        esmvaltool meta data dict
    strict, optional
        Raise error if not exactly one match exists, by default True

    Returns
    -------
        Dict: One value from the input meta

    Raises
    ------
    ValueError
        Too many matching entries
    ValueError
        No matching entry
    """
    selected_meta = select_metadata(meta, **kwargs)
    if len(selected_meta) > 1:
        log.warning("Multiple entries found for Metadata: %s", selected_meta)
        if strict:
            raise ValueError("Too many matching entries")
    elif len(selected_meta) == 0:
        log.warning("No Metadata found! For: %s", kwargs)
        if strict:
            raise ValueError("No matching entry")
        return None
    return selected_meta[0]


select_single_meta = select_single_metadata


def sub_cfg(cfg: dict, plot: str, key: str) -> dict:
    """Get get merged general and plot type specific kwargs."""
    if isinstance(cfg.get(key, {}), dict):
        general = cfg.get(key, {}).copy()
        specific = cfg.get(plot, {}).get(key, {})
        general.update(specific)
        return general
    try:
        return cfg[plot][key]
    except KeyError:
        return cfg[key]


def guess_experiment(meta: dict) -> None:
    """Guess missing 'exp' in metadata from filename."""
    exps = ["historical", "ssp126", "ssp245", "ssp370", "ssp585"]
    for exp in exps:
        if exp in meta["filename"]:
            meta["exp"] = exp


# PLOT HELPER


def date_tick_layout(
    fig,
    axes,
    dates: list | None = None,
    label: str = "Time",
    years: int | None = 1,
) -> None:
    """Update a time series figure to use date/year ticks and grid."""
    axes.set_xlabel(label)
    if dates is not None:
        datemin = np.datetime64(dates[0], "Y")
        datemax = np.datetime64(dates[-1], "Y") + np.timedelta64(1, "Y")
        axes.set_xlim(datemin, datemax)
    if years is None:
        locator = mdates.AutoDateLocator()
        min_locator = mdates.YearLocator(1)
    else:
        locator = mdates.YearLocator(years)  # type: ignore[assignment]
        min_locator = mdates.YearLocator(1)
    year_formatter = mdates.DateFormatter("%Y")
    axes.grid(True)
    axes.xaxis.set_major_locator(locator)
    axes.xaxis.set_major_formatter(year_formatter)
    axes.xaxis.set_minor_locator(min_locator)
    fig.autofmt_xdate()  # align, rotate and space for tick labels


def auto_tick_layout(fig, axes, dates=None) -> None:
    """Update a time series figure to use auto date ticks and grid.

    NOTE: can this be merged with date_tick_layout?
    """
    axes.set_xlabel("Time")
    if dates is not None:
        datemin = np.datetime64(dates[0], "Y")
        datemax = np.datetime64(dates[-1], "Y") + np.timedelta64(1, "Y")
        axes.set_xlim(datemin, datemax)
    year_locator = mdates.YearLocator(1)
    months_locator = mdates.MonthLocator()
    year_formatter = mdates.DateFormatter("%Y")
    axes.grid(True)
    axes.xaxis.set_major_locator(year_locator)
    axes.xaxis.set_major_formatter(year_formatter)
    axes.xaxis.set_minor_locator(months_locator)
    fig.autofmt_xdate()  # align, rotate and space for tick labels


def map_land_layout(axes: GeoAxes, plot, bounds, *, cbar: bool = True) -> None:
    """Plot style for rectangular drought maps with land overlay.

    Mask the ocean by overlay, add gridlines, set left/bottom tick labels.
    """
    axes.coastlines()
    axes.add_feature(
        ct.feature.OCEAN,
        edgecolor="black",
        facecolor="white",
        zorder=1,
    )
    glines = axes.gridlines(
        crs=ct.crs.PlateCarree(),
        linewidth=1,
        color="black",
        alpha=0.6,
        linestyle="--",
        draw_labels=True,
        zorder=2,
    )
    glines.xlabels_top = False
    glines.ylabels_right = False
    if bounds is not None and cbar:
        plt.colorbar(
            plot,
            ax=axes,
            ticks=bounds,
            extend="both",
            fraction=0.022,
        )
    elif cbar:
        plt.colorbar(plot, ax=axes, extend="both", fraction=0.022)


def font_color(background: str | tuple | float) -> str:
    """Return black or white depending on backgrounds greyscale.

    Parameters
    ----------
    background : str, tuple, or float
        color as string (grayscale value, name, hex) or tuple (rgb, rgba)
    """
    if sum(mpl.colors.to_rgb(background)) > 1.5:
        return "black"
    return "white"
