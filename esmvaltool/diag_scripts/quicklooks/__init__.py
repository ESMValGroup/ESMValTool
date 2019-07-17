"""Shared scripts for quicklook diagnostics."""

import logging

import cftime
import iris
import matplotlib.pyplot as plt
from cf_units import Unit

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_plot_filename, group_metadata,
                                            io)

logger = logging.getLogger(__name__)


def _convert_units(cube, target_units):
    """Convert units of cube if possible."""
    cube_units = cube.units

    # Try conversion
    try:
        cube.convert_units(target_units)
    except ValueError:
        if cube_units == Unit('kg m-2 s-1'):
            cube.units = Unit('mm s-1')

    # Try conversion again, might work for precipitation-like units
    try:
        cube.convert_units(target_units)
    except ValueError:
        logger.warning(
            "Cannot convert cube units '%s' to desired plot units '%s'",
            cube_units, target_units)
        cube.units = cube_units


def _guess_calendar_datetime(cube):
    """Guess the cftime.datetime form to create datetimes.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.

    Returns
    -------
    cftime.datetime
        A datetime creator function from cftime, based on the cube's calendar.
    """
    time_coord = cube.coord('time')

    if time_coord.units.calendar in [
            '360_day',
    ]:
        datetime = cftime.Datetime360Day
    elif time_coord.units.calendar in ['365_day', 'noleap']:
        datetime = cftime.DatetimeNoLeap
    elif time_coord.units.calendar in [
            'julian',
    ]:
        datetime = cftime.DatetimeJulian
    elif time_coord.units.calendar in [
            'gregorian',
    ]:
        datetime = cftime.DatetimeGregorian
    elif time_coord.units.calendar in [
            'proleptic_gregorian',
    ]:
        datetime = cftime.DatetimeProlepticGregorian
    else:
        logger.warning('Calendar set to Gregorian, instead of %s',
                       time_coord.units.calendar)
        datetime = cftime.DatetimeGregorian
    return datetime


def cube_time_to_float(cube):
    """Convert time coordinate into decimal time.

    Takes an iris time coordinate and returns a list of floats.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube.

    Returns
    -------
    list
        List of floats showing the time coordinate in decimal time.

    """
    times = cube.coord('time')
    datetime = _guess_calendar_datetime(cube)

    dtimes = times.units.num2date(times.points)
    floattimes = []
    for dtime in dtimes:
        # TODO: it would be better to have a calendar dependent value
        # for daysperyear, as this is not accurate for 360 day calendars.
        daysperyear = 365.25

        try:
            dayofyr = dtime.dayofyr
        except AttributeError:
            time = datetime(dtime.year, dtime.month, dtime.day)
            time0 = datetime(dtime.year, 1, 1, 0, 0)
            dayofyr = (time - time0).days

        floattime = dtime.year + dayofyr / daysperyear + dtime.hour / (
            24.0 * daysperyear)
        if dtime.hour:
            floattime += dtime.hour / (24.0 * daysperyear)
        if dtime.minute:
            floattime += dtime.minute / (24.0 * 60.0 * daysperyear)
        floattimes.append(floattime)
    return floattimes


def get_grouped_input_data(cfg):
    """Get input data."""
    input_data = []
    if cfg['quicklook']['active']:
        quicklook_dir = cfg['quicklook']['preproc_dir']
        logger.info("Reading data from quicklook directory %s", quicklook_dir)
        raw_data = io.netcdf_to_metadata(cfg, root=quicklook_dir)
        for dataset in raw_data:
            if 'derive_input' in dataset['filename']:
                logger.debug("Ignoring derivation input variable '%s'",
                             dataset['filename'])
                continue
            input_data.append(dataset)
    else:
        logger.info("Reading data from regular ESMValTool preproc directory")
        input_data = cfg['input_data'].values()

    return group_metadata(input_data, 'short_name')


def get_provenance_record(caption, ancestors):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'authors': ['bock_ls'],
        'references': ['acknow_project'],
        'ancestors': ancestors,
    }
    return record


def load_cube(dataset, required_coords):
    """Load cube and compute average over all excessive diagnostics."""
    filename = dataset['filename']
    logger.debug("Loading '%s'", filename)
    cube = iris.load_cube(filename)
    coords = [coord.name() for coord in cube.coords(dim_coords=True)]

    # Check if cubes has desired coordinates
    for coord_name in required_coords:
        if coord_name not in coords:
            logger.warning(
                "File '%s' does not contain necessary coordinate '%s', "
                "skipping", filename, coord_name)
            return None
        coords.remove(coord_name)

    # Calculate mean
    if cube.ndim > len(required_coords):
        if 'latitude' in coords:
            grid_areas = iris.analysis.cartography.area_weights(cube)
        else:
            grid_areas = None
        cube = cube.collapsed(coords, iris.analysis.MEAN, weights=grid_areas)

    # If desired, convert units of cube
    if dataset.get('plot_units'):
        _convert_units(cube, dataset['plot_units'])

    return cube


def save_plot(basename, plot_type, cfg, **save_kwargs):
    """Save matplotlib plot."""
    provenance_record = {}
    if cfg['write_plots']:
        save_kwargs.setdefault('bbox_inches', 'tight')
        save_kwargs.setdefault('orientation', 'landscape')
        plot_path = get_plot_filename(basename, cfg)
        plt.savefig(plot_path, **save_kwargs)
        logger.info("Wrote %s", plot_path)
        provenance_record['plot_file'] = plot_path
        provenance_record['plot_type'] = plot_type
    plt.close()
    return provenance_record


def set_plot_title(title, y_shift=1.08):
    """Set title of plot."""
    plt.title(title, y=y_shift)


def write_provenance(netcdf_path, provenance_record, cfg):
    """Write provenance information."""
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)
