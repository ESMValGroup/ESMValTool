from __future__ import annotations
import datetime as dt
import itertools as it
import logging
import os
from contextlib import suppress
from os.path import dirname as par_dir
from pathlib import Path
from pprint import pformat
from tkinter import W

import cartopy as ct
import cartopy.crs as cart
import iris
import iris.plot as iplt
import matplotlib as mpl
import matplotlib.dates as mda
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from cf_units import Unit
from esmvalcore import preprocessor as pp
from iris.analysis import Aggregator
from iris.coords import AuxCoord
from iris.util import equalise_attributes

import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts import shared
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_cfg,
    get_diagnostic_filename,
    get_plot_filename,
    select_metadata,
    group_metadata,
)
from esmvaltool.diag_scripts.shared._base import _get_input_data_files

logger = logging.getLogger(os.path.basename(__file__))

# fmt: off
DENSITY = AuxCoord(
    1000,
    long_name="density",
    units="kg m-3")

# fmt: off
FNAME_FORMAT = "{project}_{reference_dataset}_{mip}_{exp}_{ensemble}_{short_name}_{start_year}-{end_year}"

CONTINENTAL_REGIONS = {
    "Global": ["GLO"],  # global
    "North America": ["GIC", "NWN", "NEN", "WNA", "CNA", "ENA"],
    "Central America": ["NCA", "SCA", "CAR"],
    "Southern America": ["NWS", "NSA", "NES", "SAM", "SWS", "SES", "SSA"],
    "Europe": ["NEU", "WCE", "EEU", "MED"],
    "Africa": ["SAH", "WAF", "CAF", "NEAF", "SEAF", "WSAF", "ESAF", "MDG"],
    "Asia": ["RAR", "WSB", "ESB", "RFE", "WCA", "ECA", "TIB", "EAS", "ARP", "SAS", "SEA"],
    "Australia": ["NAU", "CAU", "EAU", "SAU", "NZ", "WAN", "EAN"],
}
# fmt: on

# REGION_NAMES = {
# 'Arabian-Peninsula',
# 'Arabian-Sea',
# 'Arctic-Ocean',
# 'Bay-of-Bengal'
# 'C.Australia',
# 'C.North-America',
# 'Caribbean',
# 'Central-Africa'
# 'E.Antarctica',
# 'E.Asia',
# 'E.Australia',
# 'E.C.Asia',
# 'E.Europe'
# 'E.North-America',
# 'E.Siberia',
# 'E.Southern-Africa'
# 'Equatorial.Atlantic-Ocean',
# 'Equatorial.Indic-Ocean'
# 'Equatorial.Pacific-Ocean',
# 'Greenland/Iceland',
# 'Madagascar',
# 'Mediterranean',
# 'N.Atlantic-Ocean',
# 'N.Australia',
# 'N.Central-America'
# 'N.E.North-America',
# 'N.E.South-America',
# 'N.Eastern-Africa',
# 'N.Europe',
# 'N.Pacific-Ocean',
# 'N.South-America',
# 'N.W.North-America'
# 'N.W.South-America',
# 'New-Zealand',
# 'Russian-Arctic',
# 'Russian-Far-East',
# 'S.Asia',
# 'S.Atlantic-Ocean',
# 'S.Australia',
# 'S.Central-America',
# 'S.E.Asia',
# 'S.E.South-America',
# 'S.Eastern-Africa',
# 'S.Indic-Ocean',
# 'S.Pacific-Ocean',
# 'S.South-America',
# 'S.W.South-America',
# 'Sahara',
# 'South-American-Monsoon',
# 'Southern-Ocean',
# 'Tibetan-Plateau',
# 'W.Antarctica',
# 'W.C.Asia',
# 'W.North-America',
# 'W.Siberia',
# 'W.Southern-Africa',
# 'West&Central-Europe',
# 'Western-Africa'
# }


def merge_list_cube(cube_list, aux_name="dataset", points=None, equalize=True):
    """Merge a list of cubes into a single one with an auxiliary variable.

    Useful for applying statistics along multiple cubes. The time coordinate is
    removed by this function.

    Parameters
    ----------
    cube_list : list
        List or iterable of cubes with the same coordinates.
    aux_name : str, optional
        Name of the new auxiliary coordinate. Defaults to "dataset".
    equalize : bool, optional
        Drops differences in attributes, otherwise raises an error.
        Defaults to True.
    points : list, optional
        Set values or labels as new coordinate points.

    Returns
    -------
    iris.cube.Cube
        Merged cube with an enumerated auxiliary variable.
    """
    for ds_index, ds_cube in enumerate(cube_list):
        coord = iris.coords.AuxCoord(ds_index, long_name=aux_name)
        ds_cube.add_aux_coord(coord)
    cubes = iris.cube.CubeList(cube_list)
    logger.info("formed %s: %s", type(cubes), cubes)
    if equalize:
        removed = equalise_attributes(cubes)
        logger.info("removed different attributes: %s", removed)
    for cube in cubes:
        cube.remove_coord("time")
    merged = cubes.merge_cube()
    return merged


def fold_meta(
    cfg: dict,
    meta: dict,
    cfg_keys: list|None=None,
    meta_keys: list|None=None,
    variables: list|None=None,
) -> tuple:
    """Create combinations of meta data and data constraints.
    
    cfg["variables"] overwrites meta["short_names"].

    Parameters
    ----------
    cfg : dict
        Plot specific configuration with cfg_keys on root level.
    meta : list
        Full meta data including ancestor files.
    cfg_keys : list, optional
        Data constraints as config entries used for product.
        Defaults to ["locations", "intervals"].
    meta_keys : list, optional
        Keys for each meta used for product, short_name added automatically.
        Defaults to ["dataset", "exp"].
    variables : list, optional
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
                select_metadata(meta, short_name=variables[0]), gk
            ).keys()
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
            logger.warning("No '%s' found in plot config", ckey)
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
        A tuple containing the selected metadata and the configuration dictionary.
    """
    this_cfg = dict(zip(groups.keys(), combi))
    filter_cfg = clean_meta(this_cfg)  # remove non meta keys
    this_meta = select_metadata(meta, **filter_cfg)[0]
    return this_meta, this_cfg


def list_meta_keys(meta:list, group:dict) -> list:
    """Return a list of all keys found for a group in the meta data."""
    return list(group_metadata(meta, group).keys())


def mkplotdir(cfg: dict, dname: str|Path)-> None:
    """Create a sub directory for plots if it does not exist."""
    new_dir = Path(cfg["plot_dir"] / dname)
    if not new_dir.is_dir():
        Path.mkdir(new_dir)


def sort_cube(cube:iris.cube, coord:str="longitude")-> iris.cube:
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


def fix_longitude(cube: iris.cube) -> iris.cube.Cube:
    """Return a cube with 0 centered longitude coords.

    updating the longitude coord and sorting the data accordingly
    """
    # make sure coords are -180 to 180
    fixed_lons = [
        l if l < 180 else l - 360 for l in cube.coord("longitude").points
    ]
    try:
        cube.add_aux_coord(
            iris.coords.AuxCoord(fixed_lons, long_name="fixed_lon"), 2,
        )
    except Exception as e:
        logger.warning("TODO: hardcoded dimensions in ut.fix_longitude")
        cube.add_aux_coord(
            iris.coords.AuxCoord(fixed_lons, long_name="fixed_lon"), 1
        )
    # sort data and fixed coordinates
    cube = sort_cube(cube, coord="fixed_lon")
    # set new coordinates as dimcoords
    new_lon = cube.coord("fixed_lon")
    new_lon_dims = cube.coord_dims(new_lon)
    # Create a new coordinate which is a DimCoord.
    # NOTE: The var name becomes the dim name
    longitude = iris.coords.DimCoord.from_coord(new_lon)
    longitude.rename("longitude")
    # Remove the AuxCoord, old longitude and add the new DimCoord.
    cube.remove_coord(new_lon)
    cube.remove_coord(cube.coord("longitude"))
    cube.add_dim_coord(longitude, new_lon_dims)
    return cube


def get_datetime_coords(cube, coord="time"):
    """returns the time coordinate of the cube converted
    from cf_date to mpl compatible datetime objects.
    TODO: this seems not to work with current Iris version?
    but cf_date.dateime contains the year aswell
    TODO: the difference behaviour may be caused by different time coords in datasets
    nc.num2date default behaviour changed to use only_cf_times, change back with
    only_use_python_datetimes=True
    """
    tc = cube.coord(coord)
    # logger.info(tc.units)
    fixed = iplt._fixup_dates(tc, tc.points)
    # logger.info(tc.points[0:12])
    # fixed = num2date(tc.points, str(tc.units), only_use_python_datetimes=True)
    # logger.info(fixed)
    # try:
    #     fixed = [f.datetime for f in fixed]
    # except Exception as e:
    #     logger.info("probably already in datetime format?")
    #     logger.info(e)
    return fixed


def get_start_year(cube, coord="time"):
    """returns the year of the first time coord point.
    Works for datetime and cf_calendar types
    """
    tc = cube.coord(coord)
    first = iplt._fixup_dates(tc, tc.points)[0]
    try:
        return first.datetime.year
    except:
        return first.year


def get_meta_list(meta, group, select=None):
    """list all entries found for the group key as a list.
    with a given selection, the meta data will be filtered first.

    Args:
        meta (dict): full meta data
        group (str): key to search for. Defaults to "alias"
        select (dict, optional): dict like {'short_name': 'pdsi'} that is passed to a selection.

    Returns:
        list: of collected values for group key.
    """
    if select is not None:
        meta = select_metadata(meta, **select)
    return list(group_metadata(meta, group).keys())


def get_datasets(cfg):
    metadata = cfg["input_data"].values()
    return group_metadata(metadata, "dataset").keys()


def get_dataset_scenarios(cfg):
    """
    returns iterable of meta data for all dataset/scenario combinations
    :param cfg:
    :return:
    """
    metadata = cfg["input_data"].values()
    input_datasets = group_metadata(metadata, "dataset").keys()
    input_scenarios = group_metadata(metadata, "alias").keys()
    return it.product(input_datasets, input_scenarios)


def date_tick_layout(fig, ax, dates=None, label="Time", auto=True, years=1):
    """update a figure (timeline) to use
    date/year ticks and grid.
    :param fig: figure that will be updated
    :param ax: ax to set ticks/labels/limits on
    :param dates: optional, to set limits
    :param label: ax label
    :param auto: if true auto format instead of year
    :param years: if not auto, tick every x years
    :return: nothing, updates figure in place
    """
    ax.set_xlabel(label)
    if dates is not None:
        datemin = np.datetime64(dates[0], "Y")
        datemax = np.datetime64(dates[-1], "Y") + np.timedelta64(1, "Y")
        ax.set_xlim(datemin, datemax)
    if auto:
        locator = mdates.AutoDateLocator()
        min_locator = mdates.YearLocator(1)
    else:
        locator = mdates.YearLocator(years)
        min_locator = mdates.YearLocator(1)
    year_formatter = mdates.DateFormatter("%Y")
    ax.grid(True)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(year_formatter)
    ax.xaxis.set_minor_locator(min_locator)
    fig.autofmt_xdate()  # align, rotate and space for tick labels


def auto_tick_layout(fig, ax, dates=None):
    ax.set_xlabel("Time")
    if dates is not None:
        datemin = np.datetime64(dates[0], "Y")
        datemax = np.datetime64(dates[-1], "Y") + np.timedelta64(1, "Y")
        ax.set_xlim(datemin, datemax)
    year_locator = mdates.YearLocator(1)
    months_locator = mdates.MonthLocator()
    year_formatter = mdates.DateFormatter("%Y")
    ax.grid(True)
    ax.xaxis.set_major_locator(year_locator)
    ax.xaxis.set_major_formatter(year_formatter)
    ax.xaxis.set_minor_locator(months_locator)
    fig.autofmt_xdate()  # align, rotate and space for tick labels


def map_land_layout(fig, ax, plot, bounds, var, cbar=True):
    """plot style for rectangular drought maps
    masks the ocean by overlay,
    add gridlines,
    set left/bottom tick labels
    TODO: make this a method of MapPlot class
    """
    ax.coastlines()
    ax.add_feature(
        ct.feature.OCEAN, edgecolor="black", facecolor="white", zorder=1
    )
    gl = ax.gridlines(
        crs=ct.crs.PlateCarree(),
        linewidth=1,
        color="black",
        alpha=0.6,
        linestyle="--",
        draw_labels=True,
        zorder=2,
    )
    gl.xlabels_top = False
    gl.ylabels_right = False
    if bounds is not None and cbar:
        cb = plt.colorbar(
            plot, ax=ax, ticks=bounds, extend="both", fraction=0.022
        )
        # cb = plt.colorbar(plot, ax=ax, ticks=bounds[0:-1:]+[bounds[-1]], extend="both", fraction=0.022)
    elif cbar:
        cb = plt.colorbar(plot, ax=ax, extend="both", fraction=0.022)
    # cb.set_label(f"{var.upper()}")
    # fig.tight_layout()


def get_cubes_dataset_alias(cfg):
    pass


def get_scenarios(meta, **kwargs):
    selected = select_metadata(meta, **kwargs)
    return list(group_metadata(selected, "alias").keys())


def add_preprocessor_input(cfg):
    """
    NOTE: One can add preprocessor output/variables as ancestor too, no need
    to do it in the diagnostic...
    """
    logger.warning(
        "Please add variables to the ancestor list instead. This function will be removed in the future."
    )
    run_dir = os.path.dirname(cfg["run_dir"])
    pp_dir = "/preproc/".join(run_dir.rsplit("/run/", 1))
    fake_cfg = {"input_files": [pp_dir]}
    cfg["input_data"].update(_get_input_data_files(fake_cfg))


def add_ancestor_input(cfg):
    """Read ancestors settings.yml and
    add it's input_data to this config.
    TODO: make sure it don't break for non ancestor scripts
    TODO: recursive? (optional)
    """
    logger.info(f"add ancestors for {cfg[n.INPUT_FILES]}")
    for input_file in cfg[n.INPUT_FILES]:
        cfg_anc_file = (
            "/run/".join(input_file.rsplit("/work/", 1)) + "/settings.yml"
        )
        cfg_anc = get_cfg(cfg_anc_file)
        cfg["input_data"].update(_get_input_data_files(cfg_anc))


def add_meta_files(cfg):
    pass


def _fixup_dates(coord, values):
    """copy from iris plot.py source code"""
    if coord.units.calendar is not None and values.ndim == 1:
        # Convert coordinate values into tuples of
        # (year, month, day, hour, min, sec)
        dates = [coord.units.num2date(val).timetuple()[0:6] for val in values]
        if coord.units.calendar == "gregorian":
            r = [dt.datetime(*date) for date in dates]
        else:
            try:
                import cftime
                import nc_time_axis
            except ImportError:
                msg = (
                    "Cannot plot against time in a non-gregorian "
                    'calendar, because "nc_time_axis" is not available :  '
                    "Install the package from "
                    "https://github.com/SciTools/nc-time-axis to enable "
                    "this usage."
                )
                raise iris.IrisError(msg)

            r = [
                nc_time_axis.CalendarDateTime(
                    cftime.datetime(*date), coord.units.calendar
                )
                for date in dates
            ]
        values = np.empty(len(r), dtype=object)
        values[:] = r
    return values


def quick_save(cube, name, cfg):
    """Simply save cube to netcdf file without additional information

    Args:
        cube: iris.cube object to save
        name: basename for the created netcdf file
        cfg: user configuration containing file output path
    """
    if cfg.get("write_netcdf", True):
        diag_file = get_diagnostic_filename(name, cfg)
        logger.info(f"quick save {diag_file}")
        iris.save(cube, target=diag_file)


def quick_load(cfg, context, strict=True):
    """
    select input files from config wich matches the selection. loads and returns the first match.
    raises an error (if strict) or a warning for multiple matches.
    """
    meta = cfg["input_data"].values()
    var_meta = select_metadata(meta, **context)
    if len(var_meta) != 1:
        if strict:
            raise Exception()
        # logger.warning(f"Unexpected amount of matching meta data: {len(var_meta)}")
        print("warning meta data missmatch")
    return iris.load_cube(var_meta[0]["filename"])


def smooth(dat, window=32, mode="same"):
    # smooth
    # from scipy.ndimage.filters import uniform_filter1d as unifilter
    # TODO: iris can also directly filter on cubes:
    # https://scitools-iris.readthedocs.io/en/latest/generated/gallery/general/plot_SOI_filtering.html#sphx-glr-generated-gallery-general-plot-soi-filtering-py
    # smoothed = unifilter(dat, 32)  # smooth over 4 year window
    # using numpy convol
    filter = np.ones(window)
    return np.convolve(dat, filter, mode) / window


def get_basename(cfg, meta, prefix=None, suffix=None):
    formats = {  # TODO: load this from config-developer.yml?
        "CMIP6": "{project}_{dataset}_{mip}_{exp}_{ensemble}_{short_name}_{grid}_{start_year}-{end_year}",
        "OBS": "{project}_{dataset}_{type}_{version}_{mip}_{short_name}_{start_year}-{end_year}",
    }
    basename = formats[meta["project"]].format(**meta)
    if suffix:
        basename += f"_{suffix}"
    if prefix:
        basename = f"{prefix}_{basename}"
    return basename


def get_custom_basename(meta, folder="plots", prefix=None, suffix=None):
    """manually create basenames for output files

    NOTE: for basename in configured naming format use get_basename
    """
    defaults = {
        "variable": "variable",
        "alias": "alias",
    }
    defaults.update(meta)
    base = f"{folder}/"
    if prefix:
        base += f"{prefix}_"
    base += "{alias}_{variable}".format(**defaults)
    if suffix:
        base += f"_{suffix}"
    return base


def clean_meta(meta, **kwargs):
    # TODO: incomplete keylist
    valid_keys = ["short_name", "dataset", "alias", "exp"]
    return {key: val for key, val in meta.items() if key in valid_keys}


def select_single_metadata(meta, strict=True, **kwargs):
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
    # logger.info("dataset from alias: " + dataset)
    if len(selected_meta) > 1:
        logger.warning(f"Multiple entries found for Metadata: {selected_meta}")
        if strict:
            raise ValueError("Too many matching entries")
    elif len(selected_meta) == 0:
        logger.warning(f"No Metadata found! For: {kwargs}")
        logger.debug(f"{meta}")
        if strict:
            raise ValueError("No matching entry")
        return None
    return selected_meta[0]


def select_single_meta(*args, **kwargs):
    return select_single_metadata(*args, **kwargs)


def date_to_months(date, start_year):
    """translates datestings YYYY-MM to
    total months since begin of start_year
    """
    years, months = [int(x) for x in date.split("-")]
    return 12 * (years - start_year) + months


def fix_interval(interval):
    """Ensure that an interval has a label and a range.
    TODO: replace "/" with "_" in diagnostics who use this.
    """
    if "range" not in interval:
        interval["range"] = f"{interval['start']}/{interval['end']}"
    if "label" not in interval:
        interval["label"] = interval["range"]
    return interval


def get_plot_filename(cfg, basename, meta=dict(), replace=dict()):
    """Get a valid path for saving a diagnostic plot.
    This is an alternative to shared.get_diagnostic_filename.
    It uses cfg as first argument and accept metadata to format the basename.

    Parameters
    ----------
    cfg: dict
        Dictionary with diagnostic configuration.
    basename: str
        The basename of the file.
    meta: dict
        Metadata to format the basename.
    replace: dict
        Dictionary with strings to replace in the basename.

    Returns
    -------
    str:
        A valid path for saving a diagnostic plot.
    """
    basename = basename.format(**meta)
    for key, value in replace.items():
        basename = basename.replace(key, value)
    return os.path.join(
        cfg["plot_dir"],
        f"{basename}.{cfg['output_file_type']}",
    )


def slice_cube_interval(cube, interval):
    """returns a cube slice for given interval
    which is a list of strings (YYYY-mm) or int (index of cube)
    NOTE: For 3D cubes (time first)
    """
    if isinstance(interval[0], int) and isinstance(interval[1], int):
        return cube[interval[0] : interval[1], :, :]
    dt_start = dt.datetime.strptime(interval[0], "%Y-%m")
    dt_end = dt.datetime.strptime(interval[1], "%Y-%m")
    time = cube.coord("time")
    t_start = time.nearest_neighbour_index(time.units.date2num(dt_start))
    t_end = time.nearest_neighbour_index(time.units.date2num(dt_end))
    return cube[t_start:t_end, :, :]


def find_first(nparr):
    """finds first nonzero element of numpy array and returns index or -1.
    Its faster than looping or using numpy.where.
    https://stackoverflow.com/a/61117770/7268121
    nparr: (numpy.array) requires numerical data without negative zeroes.
    """
    idx = nparr.view(bool).argmax() // nparr.itemsize
    return idx if np.arr[idx] else -1


def add_spei_meta(cfg, name="spei", pos=0):
    """NOTE: workaround to add missing meta for specific ancestor script"""
    logger.info("adding meta file for save_spei output")
    spei_fname = (
        cfg["tmp_meta"]["filename"].split("/")[-1].replace("_pr_", f"_{name}_")
    )
    # FIXME: hardcoded index in input files.. needs to match recipe
    spei_file = os.path.join(cfg["input_files"][pos], spei_fname)
    logger.info(f"spei file path: {spei_file}")
    logger.info("generating missing meta info")
    meta = cfg["tmp_meta"].copy()
    meta["filename"] = spei_file
    meta["short_name"] = name
    meta["long_name"] = "Standardised Precipitation-Evapotranspiration Index"
    # meta['standard_name'] = 'standardised_precipitation_evapotranspiration_index'
    if name.lower() == "spi":
        meta["long_name"] = "Standardised Precipitation Index"
        # meta['standard_name'] = 'standardised_precipitation_index'
    cfg["input_data"][spei_file] = meta


def fix_calendar(cube):
    """Fix Calendar.
    Fixes wrong calendars to 'gregorian' instead of 'proleptic_gregorian' or '365_day' or any other.
    TODO: use pp.regrid_time when available in esmvalcore.

    Args:
        cube (iris.cube.Cube): Input cubes which need to be fixed.

    Returns:
        iris.cube.Cube: with fixed calendars

    """
    time = cube.coord("time")
    # if time.units.name == "days since 1850-1-1 00:00:00":
    logger.info("renaming unit")
    time.units = Unit("days since 1850-01-01", calendar=time.units.calendar)
    if time.units.calendar == "proleptic_gregorian":
        time.units = Unit(time.units.name, calendar="gregorian")
        logger.info(f"renamed calendar: {time.units.calendar}")
    if time.units.calendar != "gregorian":
        time.convert_units(Unit("days since 1850-01-01", calendar="gregorian"))
        logger.info(f"converted time to: {time.units}")
    return cube


def latlon_coords(cube):
    """replace latitude and longitude
    with lat and lon inplace
    TODO: make this good!
    """
    try:
        cube.coord("longitude").rename("lon")
    except Exception as e:
        logging.info("no coord named longitude")
        print(e)
    try:
        cube.coord("latitude").rename("lat")
    except:
        logging.info("no coord named latitude")


def standard_time(cubes):
    """Make sure all cubes' share the standard time coordinate.
    This function extracts the date information from the cube and
    reconstructs the time coordinate, resetting the actual dates to the
    15th of the month or 1st of july for yearly data (consistent with
    `regrid_time`), so that there are no mismatches in the time arrays.
    It will use reset the calendar to
    a default gregorian calendar with unit "days since 1850-01-01".
    Might not work for (sub)daily data, because different calendars may have
    different number of days in the year.
    NOTE: this might be replaced by preprocessor
    """
    from datetime import datetime

    from esmvalcore.iris_helpers import date2num

    t_unit = Unit("days since 1850-01-01", calendar="standard")

    for cube in cubes:
        # Extract date info from cube
        coord = cube.coord("time")
        years = [p.year for p in coord.units.num2date(coord.points)]
        months = [p.month for p in coord.units.num2date(coord.points)]
        days = [p.day for p in coord.units.num2date(coord.points)]

        # Reconstruct default calendar
        if 0 not in np.diff(years):
            # yearly data
            dates = [datetime(year, 7, 1, 0, 0, 0) for year in years]
        elif 0 not in np.diff(months):
            # monthly data
            dates = [
                datetime(year, month, 15, 0, 0, 0)
                for year, month in zip(years, months)
            ]
        elif 0 not in np.diff(days):
            # daily data
            dates = [
                datetime(year, month, day, 0, 0, 0)
                for year, month, day in zip(years, months, days)
            ]
            if coord.units != t_unit:
                logger.warning(
                    "Multimodel encountered (sub)daily data and inconsistent "
                    "time units or calendars. Attempting to continue, but "
                    "might produce unexpected results."
                )
        else:
            raise ValueError(
                "Multimodel statistics preprocessor currently does not "
                "support sub-daily data."
            )

        # Update the cubes' time coordinate (both point values and the units!)
        cube.coord("time").points = date2num(dates, t_unit, coord.dtype)
        cube.coord("time").units = t_unit
        cube.coord("time").bounds = None
        cube.coord("time").guess_bounds()


def guess_lat_lon_bounds(cube):
    """guesses bounds for latitude and longitude if not existent."""
    if not cube.coord("latitude").has_bounds():
        cube.coord("latitude").guess_bounds()
    if not cube.coord("longitude").has_bounds():
        cube.coord("longitude").guess_bounds()


def mmm(cube_list, mdtol=0, dropcoords=["time"], dropmethods=False):
    """calculates mean and stdev along a cube list over all cubes returns two (mean and stdev) of same shape
    TODO: merge alreadey exist, use that one, mean and std is trivial than.
    Args:
        cube_list ([type]): [description]
    """
    for i, cube in enumerate(cube_list):
        # code.interact(local=locals())
        for coord_name in dropcoords:
            if cube.coords(coord_name):
                cube.remove_coord(coord_name)
        if dropmethods:
            cube.cell_methods = None
        cube.add_aux_coord(iris.coords.AuxCoord(i, long_name="dataset"))
    # equalise_attributes(cube_list)
    # cube.remove_coord("season_number")  # add as drop_coord
    cube_list = iris.cube.CubeList(cube_list)
    equalise_attributes(cube_list)
    # common_depth =
    try:
        # Note: just for testing:
        merged = cube_list.merge_cube()
    except iris.exceptions.MergeError as err:
        # unique_seasons = set()
        # for c in cube_list:
        #     unique_seasons.add(c.coord("season_number"))
        # print(unique_seasons)
        # for c in cube_list:
        #     print(c.coord("depth"))
        iris.util.describe_diff(cube_list[0], cube_list[1])
        raise err
    if mdtol > 0:
        logger.info(f"performing MMM with tolerance: {mdtol}")
    mean = merged.collapsed("dataset", iris.analysis.MEAN, mdtol=mdtol)
    sdev = merged.collapsed("dataset", iris.analysis.STD_DEV)
    return mean, sdev


def get_hex_positions():
    return {
        "NWN": [2, 0],
        "NEN": [4, 0],
        "GIC": [6.5, -0.5],
        "NEU": [14, 0],
        "RAR": [20, 0],
        "WNA": [1, 1],
        "CNA": [3, 1],
        "ENA": [5, 1],
        "WCE": [13, 1],
        "EEU": [15, 1],
        "WSB": [17, 1],
        "ESB": [19, 1],
        "RFE": [21, 1],
        "NCA": [2, 2],
        "MED": [14, 2],
        "WCA": [16, 2],
        "ECA": [18, 2],
        "TIB": [20, 2],
        "EAS": [22, 2],
        "SCA": [3, 3],
        # "CAR": [5, 3],
        "SAH": [13, 3],
        "ARP": [15, 3],
        "SAS": [19, 3],
        "SEA": [23, 3],
        # "PAC": [27.5, 3.3],
        "NWS": [6, 4],
        "NSA": [8, 4],
        "WAF": [12, 4],
        "CAF": [14, 4],
        "NEAF": [16, 4],
        "NAU": [24.5, 4.3],
        "SAM": [7, 5],
        "NES": [9, 5],
        "WSAF": [13, 5],
        "SEAF": [15, 5],
        "MDG": [17.5, 5.3],
        "CAU": [23.5, 5.3],
        "EAU": [25.5, 5.3],
        "SWS": [6, 6],
        "SES": [8, 6],
        "ESAF": [14, 6],
        "SAU": [24.5, 6.3],
        "NZ": [27, 6.5],
        "SSA": [7, 7],
    }


def get_region_data():
    """reads shapes.txt as csv file and returns a list of region names"""
    fname = "/work/bd0854/b309169/ESMValTool-private/esmvaltool/diag_scripts/droughtindex/shapes.txt"
    data = pd.read_csv(fname, sep=",", header=0, skiprows=0)
    print(data)
    return data


def get_region_abbrs():
    # abbr_positions = get_hex_positions()
    data = get_region_data()
    data.set_index("Name", inplace=True)
    return {i: data.loc[i]["Abbr"] for i in data.index.tolist()}


def add_aux_regions(cube, mask=None):
    """Add an auxilary coordinate with region numbers.
    Dimension names: 'lat' and 'lon'.
    TODO: add option/fallback to use pp with shapefile instead
    """
    import regionmask

    region = regionmask.defined_regions.ar6.land
    xarr = xr.DataArray.from_iris(cube)
    mask2d = mask if mask else region.mask(xarr).to_iris()
    # newcube = cube.copy()
    region_coord = iris.coords.AuxCoord(
        mask2d.data, long_name="AR6 reference region", var_name="region"
    )
    cube.add_aux_coord(region_coord, [1, 2])
    return cube


def regional_mean(
    cube,
    mask=None,
    minmax=False,
    stddev=False,
):
    """Returns a cube with one dim_coord 'region'
    which contains weighted means over 'lat', 'lon'.
    TODO: allow minmax and stddev at the same time. Different returns or
    return a dict? maybe add an metrics array: [mean, min, max, sum ...]
    TODO: provide an alternative based on shape file and preprocessor in case
    regionmask is not installed?
    TODO: pass config and/or aux file path?
    TODO: DEPRECATED use more general regiona_stats instead (includes maean)
    """
    logger.warning("Please use regional_stats instead of regional_mean")
    if minmax:
        res = regional_stats_xarr(
            {}, cube, operators=["min", "mean", "max"]
        ).values()
        return res["min", "mean", "max"]
    if stddev:
        res = regional_stats_xarr({}, cube, operators["mean", "std_dev"])
        return res["mean"], res["std_dev"]
    res = regional_stats_xarr(
        {}, cube, operators=["min", "mean", "max"]
    ).values()
    return res["mean"]


def regional_stats_xarr(cfg, cube, operators=["mean"], shapefile=None):
    results = {}
    if shapefile:
        return regional_stats(cfg, cube, operators, shapefile)
    try:
        import regionmask
    except:
        err = "No Module regionmask. Install it or provide a shapefile."
        raise ModuleError(err)
    latlon_coords(cube)
    region = regionmask.defined_regions.ar6.land
    xarr = xr.DataArray.from_iris(cube)
    mask3d = mask if mask else region.mask_3D(xarr)
    weights = np.cos(np.deg2rad(xarr.lat))  # TODO: use cell area
    weights3d = mask3d * weights
    weighted = xarr.weighted(weights3d)
    if "mean" in operators:
        mean = weighted.mean(dim=("lat", "lon"))
        result["mean"] = (mean.to_iris(),)
    # stddev = weighted.sum_of_squares() / weights3d.sum()
    if "std_dev" in operators:
        stddev = weighted.std(dim=("lat", "lon"))
        result["std_dev"] = stddev.to_iris()
    if "min" in operators:
        mini = xarr.where(mask3d).min(dim=("lat", "lon"))
        result["min"] = mini.to_iris()
    if "max" in operators:
        maxi = xarr.where(mask3d).max(dim=("lat", "lon"))
        result["max"] = maxi.to_iris()
    return result


def regional_stats(cfg, cube, operator="mean"):
    """Mean over IPCC Reference regions using shape file.
    The shapefile (string) must be contained in cft (given by recipe).
    It can be an absolute path or relative to the folder for auxilary data
    configured in esmvaltool config.
    """
    # if "shapefile" not in cfg:
    #     raise ValueError("A shapefile must be given or \
    #         utils.regional_stats_xarr() be used with module regionmask.")
    # shapefile = Path(cfg['auxiliary_data_dir']) / cfg['shapefile']
    guess_lat_lon_bounds(cube)
    extracted = pp.extract_shape(cube, "ar6", decomposed=True)

    return pp.area_statistics(extracted, operator)


def transpose_by_names(cube, names):
    """transposes an iris cube by dim-coords or their names"""
    new_dims = [cube.coord_dims(name)[0] for name in names]
    print(new_dims)
    cube.transpose(new_dims)


def generate_metadata(work_folder, diag=None):
    """create metadata from ancestor config and nc files.

    Can be used as workaround, to handle output of ancestors, that don't create
    metadata.yml, similar to preprocessed variables. If metadata.yml exists in
    the work folder (subfolders ignored) its content is returned instead.
    """
    try:
        return yaml.load(os.path.join(work_folder), "metadata.yml")
    except:
        logger.warning(f"no metadata.yml found in {work_folder}.")
    raise NotImplementedError("TODO: generate metadata from cubes and cfg")
    # TODO: provide fixes for some diagnostics or directly implement it.
    # meta = {}
    # cfg = yaml.read()
    # nc_files =


def save_metadata(cfg, metadata):
    """save dict as metadata.yml in workfolder."""
    with open(os.path.join(cfg["work_dir"], "metadata.yml"), "w") as wom:
        yaml.dump(metadata, wom)


def get_meta(index):
    index = index.lower()
    if index == "pdsi":
        return dict(
            variable_group="index",
            standard_name="palmer_drought_severity_index",
            short_name="pdsi",
            long_name="Palmer Drought Severity Index",
            units="1",
        )
    elif index in ["scpdsi", "sc-pdsi"]:
        return dict(
            variable_group="index",
            standard_name="self_calibrated_palmer_drought_severity_index",
            short_name="scpdsi",
            long_name="Self Calibrated Palmer Drought Severity Index",
            units="1",
        )
    else:
        logger.error(f"No default meta data for Index: {index}")
        raise NotImplementedError


def set_defaults(target, defaults):
    """Applies set_default on target for each entry of defaults

    This checks if a key exists in target, and only if not the keys are set with
    values. It does not checks recursively for nested entries.

    NOTE: this might be obsolete since python 3.9, as there are direct
    fallback assignments like: `target = defaults | target`

    Parameters
    ----------
    target
        dictionary to set the defaults on
    defaults
        dictionary containtaining defaults to be set on target
    """
    for key in defaults.keys():
        target.setdefault(key, defaults[key])


def sub_cfg(cfg, plot, key):
    """Get get merged general and plot type specific kwargs."""
    if isinstance(cfg.get(key, {}), dict):
        general = cfg.get(key, {}).copy()
        specific = cfg.get(plot, {}).get(key, {})
        general.update(specific)
        return general
    else:
        try:
            return cfg[plot][key]
        except KeyError:
            return cfg[key]


def aux_path(cfg, path):
    """returns absolut path of an aux file."""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(cfg["auxiliary_data_dir"], path)


def remove_attributes(cube, ignore=[]):
    """remove attributes of cubes or coords in place
    used to clean up differences in cubes coordinates before merging

    Parameters
    ----------
    cube
        iris.Cube or iris.Coord
    ignore
        Optional: List of Strings of attributes, that are not removed
        Default: []
    """
    remove = []
    for attr in cube.attributes:
        if not attr in ignore:
            remove.append(attr)
    for attr in remove:
        del cube.attributes[attr]


def convert_to_mmday_xarray(pr):
    """convert precipitation of xarray from kg/m2/s to mm/day

    Args:
        pr: xarray.dataarrray precipitation in kgm-2s-1
    """
    pr.values = pr.values * 60 * 60 * 24
    pr.attrs["units"] = "mm day-1"
    return pr


def font_color(background):
    """black or white depending on greyscale of the background

    Parameters
    ----------
    bacgkround
        matplotlib color
    """
    if sum(mpl.colors.to_rgb(background)) > 1.5:
        return "black"
    else:
        return "white"


def log_provenance(cfg, fname, record):
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(fname, record)


def get_time_range(cube):
    """guesses the period of a cube based on the time coordinate."""
    if not isinstance(cube, iris.cube.Cube):
        cube = iris.load_cube(cube)
    time = cube.coord("time")
    print(time)
    print(time.points)
    start = time.units.num2date(time.points[0])
    end = time.units.num2date(time.points[-1])
    return {"start_year": start.year, "end_year": end.year}


def guess_experiment(meta):
    """guess experiment from filename
    TODO: this is a workaround for incomplete meta data..
    fix this in ancestor diagnostics rather than use this function.
    """
    exps = ["historical", "ssp126", "ssp245", "ssp370", "ssp585"]
    for exp in exps:
        if exp in meta["filename"]:
            meta["exp"] = exp


def monthly_to_daily(cube, units="mm day-1", leap_years=True):
    """convert monthly data to daily data inplace ignoring leap years"""
    months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    months = months * int((cube.shape[0] / 12) + 1)
    for i, s in enumerate(cube.slices_over(["time"])):
        if not leap_years:
            days = months[i]
            cube.data[i] = cube.data[i] / days
        try:
            from calendar import monthrange

            time = s.coord("time")
            date = time.units.num2date(time.points[0])
            days = monthrange(date.year, date.month)[1]
        except Exception as e:
            logger.warning("date failed, using fixed days without leap year")
            logger.warning(e)
            days = months[i]
        cube.data[i] = cube.data[i] / days
    cube.units = units


def daily_to_monthly(cube, units="mm month-1", leap_years=True):
    """convert daily data to monthly data inplace
    with leap_years=False this is similar to the same named function in utils.R
    and compatible with pet.R
    """
    months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    months = months * int((cube.shape[0] / 12) + 1)
    for i, s in enumerate(cube.slices_over(["time"])):
        if not leap_years:
            days = months[i]
            cube.data[i] = cube.data[i] * days
            continue
        try:  # consider leapday
            from calendar import monthrange

            time = s.coord("time")
            date = time.units.num2date(time.points[0])
            days = monthrange(date.year, date.month)[1]
        except Exception as e:
            logger.warning("date failed, using fixed days without leap year")
            logger.warning(e)
            days = months[i]
        cube.data[i] = cube.data[i] * days
    cube.units = units


def _get_data_hlp(axis, data, ilat, ilon):
    """Get data_help dependend on axis."""
    if axis == 0:
        data_help = (data[:, ilat, ilon])[:, 0]
    elif axis == 1:
        data_help = (data[ilat, :, ilon])[:, 0]
    elif axis == 2:
        data_help = data[ilat, ilon, :]

    return data_help


def _get_drought_data(cfg, cube):
    """Prepare data and calculate characteristics."""
    # make a new cube to increase the size of the data array
    # Make an aggregator from the user function.
    spell_no = Aggregator(
        "spell_count", count_spells, units_func=lambda units: 1
    )
    new_cube = _make_new_cube(cube)

    # calculate the number of drought events and their average duration
    drought_show = new_cube.collapsed(
        "time", spell_no, threshold=cfg["threshold"]
    )
    drought_show.rename("Drought characteristics")
    # length of time series
    time_length = len(new_cube.coord("time").points) / 12.0
    # Convert number of droughtevents to frequency (per year)
    drought_show.data[:, :, 0] = drought_show.data[:, :, 0] / time_length
    return drought_show


def _provenance_map_spei(cfg, name_dict, spei, dataset_name):
    """Set provenance for plot_map_spei."""
    caption = (
        "Global map of "
        + name_dict["drought_char"]
        + " ["
        + name_dict["unit"]
        + "] "
        + "based on "
        + cfg["indexname"]
        + "."
    )

    if cfg["indexname"].lower == "spei":
        set_refs = [
            "martin18grl",
            "vicente10jclim",
        ]
    elif cfg["indexname"].lower == "spi":
        set_refs = [
            "martin18grl",
            "mckee93proc",
        ]
    else:
        set_refs = [
            "martin18grl",
        ]

    provenance_record = get_provenance_record(
        [name_dict["input_filenames"]], caption, ["global"], set_refs
    )

    diagnostic_file = get_diagnostic_filename(
        cfg["indexname"]
        + "_map"
        + name_dict["add_to_filename"]
        + "_"
        + dataset_name,
        cfg,
    )
    plot_file = get_plot_filename(
        cfg["indexname"]
        + "_map"
        + name_dict["add_to_filename"]
        + "_"
        + dataset_name,
        cfg,
    )

    logger.info("Saving analysis results to %s", diagnostic_file)

    cubesave = cube_to_save_ploted(spei, name_dict)
    iris.save(cubesave, target=diagnostic_file)

    logger.info(
        "Recording provenance of %s:\n%s",
        diagnostic_file,
        pformat(provenance_record),
    )
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_file, provenance_record)
        provenance_logger.log(diagnostic_file, provenance_record)


def _provenance_map_spei_multi(cfg, data_dict, spei, input_filenames):
    """Set provenance for plot_map_spei_multi."""
    caption = (
        "Global map of the multi-model mean of "
        + data_dict["drought_char"]
        + " ["
        + data_dict["unit"]
        + "] "
        + "based on "
        + cfg["indexname"]
        + "."
    )

    if cfg["indexname"].lower == "spei":
        set_refs = [
            "martin18grl",
            "vicente10jclim",
        ]
    elif cfg["indexname"].lower == "spi":
        set_refs = [
            "martin18grl",
            "mckee93proc",
        ]
    else:
        set_refs = [
            "martin18grl",
        ]

    provenance_record = get_provenance_record(
        input_filenames, caption, ["global"], set_refs
    )

    diagnostic_file = get_diagnostic_filename(
        cfg["indexname"]
        + "_map"
        + data_dict["filename"]
        + "_"
        + data_dict["datasetname"],
        cfg,
    )
    plot_file = get_plot_filename(
        cfg["indexname"]
        + "_map"
        + data_dict["filename"]
        + "_"
        + data_dict["datasetname"],
        cfg,
    )

    logger.info("Saving analysis results to %s", diagnostic_file)

    iris.save(cube_to_save_ploted(spei, data_dict), target=diagnostic_file)

    logger.info(
        "Recording provenance of %s:\n%s",
        diagnostic_file,
        pformat(provenance_record),
    )
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_file, provenance_record)
        provenance_logger.log(diagnostic_file, provenance_record)


def _provenance_time_series_spei(cfg, data_dict):
    """Provenance for time series plots."""
    caption = (
        "Time series of " + data_dict["var"] + " at" + data_dict["area"] + "."
    )

    if cfg["indexname"].lower == "spei":
        set_refs = [
            "vicente10jclim",
        ]
    elif cfg["indexname"].lower == "spi":
        set_refs = [
            "mckee93proc",
        ]
    else:
        set_refs = [
            "martin18grl",
        ]

    provenance_record = get_provenance_record(
        [data_dict["filename"]], caption, ["reg"], set_refs, plot_type="times"
    )

    diagnostic_file = get_diagnostic_filename(
        cfg["indexname"]
        + "_time_series_"
        + data_dict["area"]
        + "_"
        + data_dict["dataset_name"],
        cfg,
    )
    plot_file = get_plot_filename(
        cfg["indexname"]
        + "_time_series_"
        + data_dict["area"]
        + "_"
        + data_dict["dataset_name"],
        cfg,
    )
    logger.info("Saving analysis results to %s", diagnostic_file)

    cubesave = cube_to_save_ploted_ts(data_dict)
    iris.save(cubesave, target=diagnostic_file)

    logger.info(
        "Recording provenance of %s:\n%s",
        diagnostic_file,
        pformat(provenance_record),
    )
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_file, provenance_record)
        provenance_logger.log(diagnostic_file, provenance_record)


def cube_to_save_ploted(var, data_dict):
    """Create cube to prepare plotted data for saving to netCDF."""
    plot_cube = iris.cube.Cube(
        var,
        var_name=data_dict["var"],
        long_name=data_dict["drought_char"],
        units=data_dict["unit"],
    )
    plot_cube.add_dim_coord(
        iris.coords.DimCoord(
            data_dict["latitude"],
            var_name="lat",
            long_name="latitude",
            units="degrees_north",
        ),
        0,
    )
    plot_cube.add_dim_coord(
        iris.coords.DimCoord(
            data_dict["longitude"],
            var_name="lon",
            long_name="longitude",
            units="degrees_east",
        ),
        1,
    )

    return plot_cube


def cube_to_save_ploted_ts(data_dict):
    """Create cube to prepare plotted time series for saving to netCDF."""
    plot_cube = iris.cube.Cube(
        data_dict["data"],
        var_name=data_dict["var"],
        long_name=data_dict["var"],
        units=data_dict["unit"],
    )
    plot_cube.add_dim_coord(
        iris.coords.DimCoord(
            data_dict["time"], var_name="time", long_name="Time", units="month"
        ),
        0,
    )

    return plot_cube


def get_provenance_record(
    ancestor_files, caption, domains, refs, plot_type="geo"
):
    """Get Provenance record."""
    record = {
        "caption": caption,
        "statistics": ["mean"],
        "domains": domains,
        "plot_type": plot_type,
        "themes": ["phys"],
        "authors": [
            "weigel_katja",
            "adeniyi_kemisola",
        ],
        "references": refs,
        "ancestors": ancestor_files,
    }
    return record


def _make_new_cube(cube):
    """Make a new cube with an extra dimension for result of spell count."""
    new_shape = cube.shape + (4,)
    new_data = iris.util.broadcast_to_shape(cube.data, new_shape, [0, 1, 2])
    new_cube = iris.cube.Cube(new_data)
    new_cube.add_dim_coord(
        iris.coords.DimCoord(cube.coord("time").points, long_name="time"), 0
    )
    new_cube.add_dim_coord(
        iris.coords.DimCoord(
            cube.coord("latitude").points, long_name="latitude"
        ),
        1,
    )
    new_cube.add_dim_coord(
        iris.coords.DimCoord(
            cube.coord("longitude").points, long_name="longitude"
        ),
        2,
    )
    new_cube.add_dim_coord(
        iris.coords.DimCoord([0, 1, 2, 3], long_name="z"), 3
    )
    return new_cube


def _plot_multi_model_maps(
    cfg, all_drought_mean, lats_lons, input_filenames, tstype
):
    """Prepare plots for multi-model mean."""
    data_dict = {
        "latitude": lats_lons[0],
        "longitude": lats_lons[1],
        "model_kind": tstype,
    }
    if tstype == "Difference":
        # RCP85 Percentage difference
        data_dict.update({
            "data": all_drought_mean[:, :, 0],
            "var": "diffnumber",
            "datasetname": "Percentage",
            "drought_char": "Number of drought events",
            "unit": "%",
            "filename": "Percentage_difference_of_No_of_Events",
            "drought_numbers_level": np.arange(-100, 110, 10),
        })
        plot_map_spei_multi(
            cfg, data_dict, input_filenames, colormap="rainbow"
        )

        data_dict.update({
            "data": all_drought_mean[:, :, 1],
            "var": "diffduration",
            "drought_char": "Duration of drought events",
            "filename": "Percentage_difference_of_Dur_of_Events",
            "drought_numbers_level": np.arange(-100, 110, 10),
        })
        plot_map_spei_multi(
            cfg, data_dict, input_filenames, colormap="rainbow"
        )

        data_dict.update({
            "data": all_drought_mean[:, :, 2],
            "var": "diffseverity",
            "drought_char": "Severity Index of drought events",
            "filename": "Percentage_difference_of_Sev_of_Events",
            "drought_numbers_level": np.arange(-50, 60, 10),
        })
        plot_map_spei_multi(
            cfg, data_dict, input_filenames, colormap="rainbow"
        )

        data_dict.update({
            "data": all_drought_mean[:, :, 3],
            "var": "diff" + (cfg["indexname"]).lower(),
            "drought_char": "Average "
            + cfg["indexname"]
            + " of drought events",
            "filename": "Percentage_difference_of_Avr_of_Events",
            "drought_numbers_level": np.arange(-50, 60, 10),
        })
        plot_map_spei_multi(
            cfg, data_dict, input_filenames, colormap="rainbow"
        )
    else:
        data_dict.update({
            "data": all_drought_mean[:, :, 0],
            "var": "frequency",
            "unit": "year-1",
            "drought_char": "Number of drought events per year",
            "filename": tstype + "_No_of_Events_per_year",
            "drought_numbers_level": np.arange(0, 0.4, 0.05),
        })
        if tstype == "Observations":
            data_dict["datasetname"] = "Mean"
        else:
            data_dict["datasetname"] = "MultiModelMean"
        plot_map_spei_multi(
            cfg, data_dict, input_filenames, colormap="gnuplot"
        )

        data_dict.update({
            "data": all_drought_mean[:, :, 1],
            "var": "duration",
            "unit": "month",
            "drought_char": "Duration of drought events [month]",
            "filename": tstype + "_Dur_of_Events",
            "drought_numbers_level": np.arange(0, 6, 1),
        })
        plot_map_spei_multi(
            cfg, data_dict, input_filenames, colormap="gnuplot"
        )

        data_dict.update({
            "data": all_drought_mean[:, :, 2],
            "var": "severity",
            "unit": "1",
            "drought_char": "Severity Index of drought events",
            "filename": tstype + "_Sev_index_of_Events",
            "drought_numbers_level": np.arange(0, 9, 1),
        })
        plot_map_spei_multi(
            cfg, data_dict, input_filenames, colormap="gnuplot"
        )
        namehlp = "Average " + cfg["indexname"] + " of drought events"
        namehlp2 = tstype + "_Average_" + cfg["indexname"] + "_of_Events"
        data_dict.update({
            "data": all_drought_mean[:, :, 3],
            "var": (cfg["indexname"]).lower(),
            "unit": "1",
            "drought_char": namehlp,
            "filename": namehlp2,
            "drought_numbers_level": np.arange(-2.8, -1.8, 0.2),
        })
        plot_map_spei_multi(
            cfg, data_dict, input_filenames, colormap="gnuplot"
        )


def _plot_single_maps(cfg, cube2, drought_show, tstype, input_filenames):
    """Plot map of drought characteristics for individual models and times."""
    cube2.data = drought_show.data[:, :, 0]
    name_dict = {
        "add_to_filename": tstype + "_No_of_Events_per_year",
        "name": tstype + " Number of drought events per year",
        "var": "frequency",
        "unit": "year-1",
        "drought_char": "Number of drought events per year",
        "input_filenames": input_filenames,
    }
    plot_map_spei(cfg, cube2, np.arange(0, 0.4, 0.05), name_dict)

    # plot the average duration of drought events
    cube2.data = drought_show.data[:, :, 1]
    name_dict.update({
        "add_to_filename": tstype + "_Dur_of_Events",
        "name": tstype + " Duration of drought events(month)",
        "var": "duration",
        "unit": "month",
        "drought_char": "Number of drought events per year",
        "input_filenames": input_filenames,
    })
    plot_map_spei(cfg, cube2, np.arange(0, 6, 1), name_dict)

    # plot the average severity index of drought events
    cube2.data = drought_show.data[:, :, 2]
    name_dict.update({
        "add_to_filename": tstype + "_Sev_index_of_Events",
        "name": tstype + " Severity Index of drought events",
        "var": "severity",
        "unit": "1",
        "drought_char": "Number of drought events per year",
        "input_filenames": input_filenames,
    })
    plot_map_spei(cfg, cube2, np.arange(0, 9, 1), name_dict)

    # plot the average spei of drought events
    cube2.data = drought_show.data[:, :, 3]

    namehlp = tstype + "_Avr_" + cfg["indexname"] + "_of_Events"
    namehlp2 = tstype + "_Average_" + cfg["indexname"] + "_of_Events"
    name_dict.update({
        "add_to_filename": namehlp,
        "name": namehlp2,
        "var": "severity",
        "unit": "1",
        "drought_char": "Number of drought events per year",
        "input_filenames": input_filenames,
    })
    plot_map_spei(cfg, cube2, np.arange(-2.8, -1.8, 0.2), name_dict)


def runs_of_ones_array_spei(bits, spei):
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


def count_spells(data, threshold, axis):
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
                    data_hits, data_help
                )

                return_var[ilat, ilon, 0] = np.count_nonzero(events)
                return_var[ilat, ilon, 1] = np.mean(events)
                return_var[ilat, ilon, 2] = np.mean(
                    (spei_sum * events)
                    / (np.mean(data_help[data_hits]) * np.mean(events))
                )
                return_var[ilat, ilon, 3] = np.mean(spei_sum / events)

    return return_var


def get_latlon_index(coords, lim1, lim2):
    """Get index for given values between two limits (1D), e.g. lats, lons."""
    index = (
        np.where(
            np.absolute(coords - (lim2 + lim1) / 2.0) <= (lim2 - lim1) / 2.0
        )
    )[0]
    return index


def plot_map_spei_multi(cfg, data_dict, input_filenames, colormap="jet"):
    """Plot contour maps for multi model mean."""
    spei = np.ma.array(data_dict["data"], mask=np.isnan(data_dict["data"]))

    # Get latitudes and longitudes from cube
    lons = data_dict["longitude"]
    if max(lons) > 180.0:
        lons = np.where(lons > 180, lons - 360, lons)
        # sort the array
        index = np.argsort(lons)
        lons = lons[index]
        spei = spei[np.ix_(range(data_dict["latitude"].size), index)]

    # Plot data
    # Create figure and axes instances
    subplot_kw = {"projection": cart.PlateCarree(central_longitude=0.0)}
    fig, axx = plt.subplots(figsize=(6.5, 4), subplot_kw=subplot_kw)
    axx.set_extent(
        [-180.0, 180.0, -90.0, 90.0], cart.PlateCarree(central_longitude=0.0)
    )

    # Draw filled contours
    cnplot = plt.contourf(
        lons,
        data_dict["latitude"],
        spei,
        data_dict["drought_numbers_level"],
        transform=cart.PlateCarree(central_longitude=0.0),
        cmap=colormap,
        extend="both",
        corner_mask=False,
    )
    # Draw coastlines
    axx.coastlines()

    # Add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.6, orientation="horizontal")

    # Add colorbar title string
    if data_dict["model_kind"] == "Difference":
        cbar.set_label(
            data_dict["model_kind"] + " " + data_dict["drought_char"] + " [%]"
        )
    else:
        cbar.set_label(
            data_dict["model_kind"] + " " + data_dict["drought_char"]
        )

    # Set labels and title to each plot
    axx.set_xlabel("Longitude")
    axx.set_ylabel("Latitude")
    axx.set_title(
        data_dict["datasetname"]
        + " "
        + data_dict["model_kind"]
        + " "
        + data_dict["drought_char"]
    )

    # Sets number and distance of x ticks
    axx.set_xticks(np.linspace(-180, 180, 7))
    # Sets strings for x ticks
    axx.set_xticklabels([
        "180°W",
        "120°W",
        "60°W",
        "0°",
        "60°E",
        "120°E",
        "180°E",
    ])
    # Sets number and distance of y ticks
    axx.set_yticks(np.linspace(-90, 90, 7))
    # Sets strings for y ticks
    axx.set_yticklabels(["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])

    fig.tight_layout()
    fig.savefig(
        get_plot_filename(
            cfg["indexname"]
            + "_map"
            + data_dict["filename"]
            + "_"
            + data_dict["datasetname"],
            cfg,
        ),
        dpi=300,
    )
    plt.close()

    _provenance_map_spei_multi(cfg, data_dict, spei, input_filenames)


def plot_map_spei(cfg, cube, levels, name_dict):
    """Plot contour map."""
    mask = np.isnan(cube.data)
    spei = np.ma.array(cube.data, mask=mask)
    np.ma.masked_less_equal(spei, 0)

    # Get latitudes and longitudes from cube
    name_dict.update({"latitude": cube.coord("latitude").points})
    lons = cube.coord("longitude").points
    lons = np.where(lons > 180, lons - 360, lons)
    # sort the array
    index = np.argsort(lons)
    lons = lons[index]
    name_dict.update({"longitude": lons})
    spei = spei[np.ix_(range(len(cube.coord("latitude").points)), index)]

    # Get data set name from cube
    try:
        dataset_name = cube.metadata.attributes["model_id"]
    except KeyError:
        try:
            dataset_name = cube.metadata.attributes["source_id"]
        except KeyError:
            dataset_name = "Observations"

    # Plot data
    # Create figure and axes instances
    subplot_kw = {"projection": cart.PlateCarree(central_longitude=0.0)}
    fig, axx = plt.subplots(figsize=(8, 4), subplot_kw=subplot_kw)
    axx.set_extent(
        [-180.0, 180.0, -90.0, 90.0], cart.PlateCarree(central_longitude=0.0)
    )

    # np.set_printoptions(threshold=np.nan)

    # Draw filled contours
    cnplot = plt.contourf(
        lons,
        cube.coord("latitude").points,
        spei,
        levels,
        transform=cart.PlateCarree(central_longitude=0.0),
        cmap="gnuplot",
        extend="both",
        corner_mask=False,
    )
    # Draw coastlines
    axx.coastlines()

    # Add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.6, orientation="horizontal")

    # Add colorbar title string
    cbar.set_label(name_dict["name"])

    # Set labels and title to each plot
    axx.set_xlabel("Longitude")
    axx.set_ylabel("Latitude")
    axx.set_title(dataset_name + " " + name_dict["name"])

    # Sets number and distance of x ticks
    axx.set_xticks(np.linspace(-180, 180, 7))
    # Sets strings for x ticks
    axx.set_xticklabels([
        "180°W",
        "120°W",
        "60°W",
        "0°",
        "60°E",
        "120°E",
        "180°E",
    ])
    # Sets number and distance of y ticks
    axx.set_yticks(np.linspace(-90, 90, 7))
    # Sets strings for y ticks
    axx.set_yticklabels(["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])

    fig.tight_layout()

    fig.savefig(
        get_plot_filename(
            cfg["indexname"]
            + "_map"
            + name_dict["add_to_filename"]
            + "_"
            + dataset_name,
            cfg,
        ),
        dpi=300,
    )
    plt.close()

    _provenance_map_spei(cfg, name_dict, spei, dataset_name)


def plot_time_series_spei(cfg, cube, filename, add_to_filename=""):
    """Plot time series."""
    # SPEI vector to plot
    spei = cube.data
    # Get time from cube
    time = cube.coord("time").points
    # Adjust (ncdf) time to the format matplotlib expects
    add_m_delta = mda.datestr2num("1850-01-01 00:00:00")
    time = time + add_m_delta

    # Get data set name from cube
    try:
        dataset_name = cube.metadata.attributes["model_id"]
    except KeyError:
        try:
            dataset_name = cube.metadata.attributes["source_id"]
        except KeyError:
            dataset_name = "Observations"

    data_dict = {
        "data": spei,
        "time": time,
        "var": cfg["indexname"],
        "dataset_name": dataset_name,
        "unit": "1",
        "filename": filename,
        "area": add_to_filename,
    }

    fig, axx = plt.subplots(figsize=(16, 4))
    axx.plot_date(
        time,
        spei,
        "-",
        tz=None,
        xdate=True,
        ydate=False,
        color="r",
        linewidth=4.0,
        linestyle="-",
        alpha=1.0,
        marker="x",
    )
    axx.axhline(y=-2, color="k")

    # Plot labels and title
    axx.set_xlabel("Time")
    axx.set_ylabel(cfg["indexname"])
    axx.set_title(
        "Mean "
        + cfg["indexname"]
        + " "
        + data_dict["dataset_name"]
        + " "
        + data_dict["area"]
    )

    # Set limits for y-axis
    axx.set_ylim(-4.0, 4.0)

    # Often improves the layout
    fig.tight_layout()
    # Save plot to file
    fig.savefig(
        get_plot_filename(
            cfg["indexname"]
            + "_time_series_"
            + data_dict["area"]
            + "_"
            + data_dict["dataset_name"],
            cfg,
        ),
        dpi=300,
    )
    plt.close()

    _provenance_time_series_spei(cfg, data_dict)
